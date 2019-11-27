// Include statements
#include "kdmig.h"

namespace rockseis {

    // Elastic: Derivate V Source gradients inside scaleGrad
// =============== ABSTRACT KDMIG CLASS =============== //
template<typename T>
Kdmig<T>::Kdmig() {
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
    freqinc=1;
}

template<typename T>
bool Kdmig<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return KDMIG_ERR;
	}else{
		Flog.close();
		return KDMIG_OK;
	}
}

template<typename T>
void Kdmig<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Kdmig<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Kdmig<T>::writeProgressbar(int x, int n, int r, int w){
	if(!logfile.empty()){
		if ( n <= 0 ) n = 1;
		if ( r > n ) r = n;
		if ( w > 48 ) w = 48;
		// Only update r times.
		if (  x % (n/r) != 0 ) return;

		// Calculuate the ratio of complete-to-incomplete.
		float ratio = x/(float)n;
		int   c     = ratio * w;
		prog.current = clock();
		float time_spent  = (float) ( prog.current - prog.previous ) / CLOCKS_PER_SEC;
		prog.previous = prog.current;
		prog.persec = 0.0;
		if(time_spent > 0) prog.persec = (n/r) / time_spent;
		snprintf(prog.speed, 48, "%.2f its/s", prog.persec); 

		// Show the percentage complete.
		sprintf(prog.progress,"%3d%% [", (int)(ratio*100) );

		// Show the load bar.
		for (x=0; x<c; x++)
			strcat(prog.progress,"=");

		for (x=c; x<w; x++)
			strcat(prog.progress," ");

		strcat(prog.progress,"]"); 
		strcat(prog.progress,"  "); 
		strcat(prog.progress, prog.speed); 
		writeLog(prog.progress); // Write to file
	}
}

template<typename T>
void Kdmig<T>::writeProgress(int x, int n, int r, int w){
	if(!logfile.empty()){
		if ( n <= 0 ) n = 1;
		if ( r > n ) r = n;
		if ( w > 48 ) w = 48;
		// Only update r times.
		if (  x % (n/r) != 0 ) return;

		// Calculuate the ratio of complete-to-incomplete.
		float ratio = x/(float)n;
		prog.current = clock();
		float time_spent  = (float) ( prog.current - prog.previous ) / CLOCKS_PER_SEC;
		prog.previous = prog.current;
		prog.persec = 0.0;
		if(time_spent > 0) prog.persec = (n/r) / time_spent;
		snprintf(prog.speed, 48, "%.2f its/s", prog.persec); 

		// Show the percentage complete.
		sprintf(prog.progress,"%3d%%", (int)(ratio*100) );
		strcat(prog.progress,"  "); 
		strcat(prog.progress, prog.speed); 
		writeLog(prog.progress); // Write to file
	}
}

template<typename T>
Kdmig<T>::~Kdmig() {
    // Nothing here
}

// =============== ACOUSTIC 2D KDMIG CLASS =============== //

template<typename T>
KdmigAcoustic2D<T>::KdmigAcoustic2D(){
    dataset = false;
    modelset = false;
    pimageset = false;
}

template<typename T>
KdmigAcoustic2D<T>::KdmigAcoustic2D(std::shared_ptr<ModelEikonal2D<T>> _model, std::shared_ptr<Ttable<T>> _ttable, std::shared_ptr<Data2D<T>> _data, std::shared_ptr<Image2D<T>> _pimage):Kdmig<T>(){
    data = _data;
    ttable = _ttable;
    model = _model;
    pimage = _pimage;
    dataset = true;
    ttableset = true;
    modelset = true;
    pimageset = true;
}

template<typename T>
int KdmigAcoustic2D<T>::run()
{
     int result = KDMIG_ERR;
     int nt;
     T dt;
	 T ot;

     nt = data->getNt();
     dt = data->getDt();
     ot = data->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<Ttable<T>> ttable_sou (new Ttable<T>(model, 1));
     std::shared_ptr<Ttable<T>> ttable_rec (new Ttable<T>(model, 1));

     /* Prepare FFT of data */
     int ntr = data->getNtrace();
     Index Idata(nt,ntr);
     T *rdata_array = data->getData();
     std::shared_ptr<Fft<T>> fft1d (new Fft<T>(2*nt));

     // Compute size of complex array
     unsigned long nf,nfs;
     nf=fft1d->getNfft();
     nfs = nf/2;
     T df=2.0*PI/(nf*dt);
     T *cdata_array;
     cdata_array = fft1d->getData();

     // Create image
     pimage->allocateImage();

     // Create ttable arrays
     ttable_sou->allocTtable();
     ttable_rec->allocTtable();

     this->writeLog("Running 2D Acoustic first arrival tomography gradient with fast sweeping method.");

     this->writeLog("Doing forward Loop.");
     // Inserting source point
     ttable_sou->insertSource(data, SMAP, 0);

     // Solving Eikonal equation for source traveltime
     ttable->interpTtable(ttable_sou);

     //Loop over data traces
     int i,j;
     for (i=0; i<ntr; i++){
         // Inserting new receiver point
         ttable_rec->insertSource(data, GMAP, i);

         // Solving Eikonal equation for receiver traveltime
         ttable->interpTtable(ttable_rec);

         /* Applying forward fourier transform over data trace */
         for(j=0; j<nt; j++){
             cdata_array[2*j] = rdata_array[Idata(j,i)];
             cdata_array[2*j+1] = 0.0;
         }
         for(j=nt; j < nf; j++){
             cdata_array[2*j] = 0.0;
             cdata_array[2*j+1] = 0.0;
         }
         fft1d->fft1d(1);

         // Build image contribution
         this->crossCorr(ttable_sou, ttable_rec, cdata_array, nfs, df, ot);
     }
        
    result=KDMIG_OK;
    return result;
}


template<typename T>
void KdmigAcoustic2D<T>::crossCorr(std::shared_ptr<Ttable<T>> ttable_sou, std::shared_ptr<Ttable<T>> ttable_rec, T *cdata, unsigned long nfs, T df, T ot) {
    /* Build image */
    if(!pimage->getAllocated()) pimage->allocateImage();
	int ix, iz, ihx, ihz, iw;
    T *TT_sou = ttable_sou->getData();
    T *TT_rec = ttable_rec->getData();
	T *imagedata = pimage->getImagedata();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nz = pimage->getNz();
	int hx, hz;
    T wsr,wsi;
    T TTsum=0;
    T omega=0;
    for (ihx=0; ihx<nhx; ihx++){
        hx= -(nhx-1)/2 + ihx;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            for (ix=0; ix<nx; ix++){
                if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
                {
                    for (iz=0; iz<nz; iz++){
                        if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz)){
                            TTsum = TT_sou[kt2D(ix-hx, iz-hz)] + TT_rec[kt2D(ix-hx, iz-hz)] - ot;

                            for (iw=0; iw<nfs; iw += this->getFreqinc()){
                                omega = iw*df;
                                wsr = cos(-omega*TTsum);
                                wsi = sin(-omega*TTsum);
                                imagedata[ki2D(ix,iz,ihx,ihz)] -= cdata[2*iw]*wsr - cdata[2*iw+1]*wsi;
                            }
                        }
                    }	
                }
            }
        }
    }
}

template<typename T>
KdmigAcoustic2D<T>::~KdmigAcoustic2D() {
    // Nothing here
}

// =============== ELASTIC 2D KDMIG CLASS =============== //

template<typename T>
KdmigElastic2D<T>::KdmigElastic2D(){
    dataset = false;
    vpmodelset = false;
    vsmodelset = false;
    simageset = false;
}

template<typename T>
KdmigElastic2D<T>::KdmigElastic2D(std::shared_ptr<ModelEikonal2D<T>> _vpmodel, std::shared_ptr<ModelEikonal2D<T>> _vsmodel, std::shared_ptr<Data2D<T>> _data, std::shared_ptr<Image2D<T>> _simage):Kdmig<T>(){
    data = _data;
    vpmodel = _vpmodel;
    vsmodel = _vsmodel;
    simage = _simage;
    dataset = true;
    vpmodelset = true;
    vsmodelset = true;
    simageset = true;
}

template<typename T>
int KdmigElastic2D<T>::run()
{
     int result = KDMIG_ERR;
     int nt;
     T dt;
	 T ot;

     nt = data->getNt();
     dt = data->getDt();
     ot = data->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<RaysAcoustic2D<T>> rays_sou (new RaysAcoustic2D<T>(vpmodel));
     std::shared_ptr<RaysAcoustic2D<T>> rays_rec (new RaysAcoustic2D<T>(vsmodel));

     /* Prepare FFT of data */
     int ntr = data->getNtrace();
     Index Idata(nt,ntr);
     T *rdata_array = data->getData();
     std::shared_ptr<Fft<T>> fft1d (new Fft<T>(2*nt));

     // Compute size of complex array
     unsigned long nf,nfs;
     nf=fft1d->getNfft();
     nfs = nf/2;
     T df=2.0*PI/(nf*dt);
     T *cdata_array;
     cdata_array = fft1d->getData();


     // Create image
     simage->allocateImage();

     this->writeLog("Running 2D Acoustic first arrival tomography gradient with fast sweeping method.");

     this->writeLog("Doing forward Loop.");
     // Inserting source point
     rays_sou->insertSource(data, SMAP, 0);

     // Solving Eikonal equation for source traveltime
     rays_sou->solve();

     //Loop over data traces
     int i,j;
     for (i=0; i<ntr; i++){
         // Reset traveltime for receiver
         rays_rec->clearTT();

         // Inserting new receiver point
         rays_rec->insertSource(data, GMAP, i);

         // Solving Eikonal equation for receiver traveltime
         rays_rec->solve();

         /* Applying forward fourier transform over data trace */
         for(j=0; j<nt; j++){
             cdata_array[2*j] = rdata_array[Idata(j,i)];
             cdata_array[2*j+1] = 0.0;
         }
         for(j=nt; j < nf; j++){
             cdata_array[2*j] = 0.0;
             cdata_array[2*j+1] = 0.0;
         }
         fft1d->fft1d(1);

         // Build image contribution
         this->crossCorr(rays_sou, rays_rec, cdata_array, nfs, df, ot, vpmodel->getLpml());
     }

        
    result=KDMIG_OK;
    return result;
}


template<typename T>
void KdmigElastic2D<T>::crossCorr(std::shared_ptr<RaysAcoustic2D<T>> rays_sou, std::shared_ptr<RaysAcoustic2D<T>> rays_rec, T *cdata, unsigned long nfs, T df, T ot, int pad) {
    /* Build image */
    if(!simage->getAllocated()) simage->allocateImage();
	int ix, iz, ihx, ihz, iw;
    T *TT_sou = rays_sou->getTT();
    T *TT_rec = rays_rec->getTT();
	T *imagedata = simage->getImagedata();
	int nhx = simage->getNhx();
	int nhz = simage->getNhz();
	int nx = simage->getNx();
	int nxt = nx+2*pad;
	int nz = simage->getNz();
	int hx, hz;
    T dx0 = simage->getDx();
    T dz0 = simage->getDz();
    T dx,dz;

    T vxxsr;
    T vxxsi;

    T vzzsr;
    T vzzsi;

    T vxzsr;
    T vxzsi;

    T vxxrr;
    T vxxri;
    T vxzrr;
    T vxzri;

    T wsr,wsi;
    //T TTsum=0;
    T omega=0;
    for (ihx=0; ihx<nhx; ihx++){
        hx= -(nhx-1)/2 + ihx;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            for (ix=0; ix<nx; ix++){
                if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
                {
                    for (iz=0; iz<nz; iz++){
                        if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz)){
                            //TTsum = TT_sou[kt2D(ix-hx+pad, iz-hz+pad)] + TT_rec[kt2D(ix-hx+pad, iz-hz+pad)] - ot;
                            
                            for (iw=1; iw<nfs; iw += this->getFreqinc()){
                                omega = iw*df;
                                //wsr = cos(-omega*TTsum);
                                //wsi = sin(-omega*TTsum);
                                //imagedata[ki2D(ix,iz,ihx,ihz)] -= cdata[2*iw]*wsr - cdata[2*iw+1]*wsi;
                                dx = dx0*omega;
                                dz = dz0*omega;
                                vzzsr =  (cos(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad+1)]) - 2*cos(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad)]) + cos(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad-1)]))/(dz*dz);
                                vxzsr =  ((cos(-omega*TT_sou[kt2D(ix-hx+pad+1, iz-hz+pad+1)]) - cos(-omega*TT_sou[kt2D(ix-hx+pad-1, iz-hz+pad+1)]))/(2*dx) - (cos(-omega*TT_sou[kt2D(ix-hx+pad+1, iz-hz+pad-1)]) - cos(-omega*TT_sou[kt2D(ix-hx+pad-1, iz-hz+pad-1)]))/(2*dx))/(2*dz); 
                                vxxrr =  cdata[2*iw]*(cos(-omega*TT_rec[kt2D(ix-hx+pad+1, iz-hz+pad)]) - cos(-omega*TT_rec[kt2D(ix-hx+pad-1, iz-hz+pad)]))/(2*dx);
                                vxzrr = cdata[2*iw]*(cos(-omega*TT_rec[kt2D(ix-hx+pad, iz-hz+pad+1)]) - cos(-omega*TT_rec[kt2D(ix-hx+pad, iz-hz+pad-1)]))/(2*dz);

                                vzzsi =  (sin(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad+1)]) - 2*sin(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad)]) + sin(-omega*TT_sou[kt2D(ix-hx+pad, iz-hz+pad-1)]))/(dz*dz);
                                vxzsi =  ((sin(-omega*TT_sou[kt2D(ix-hx+pad+1, iz-hz+pad+1)]) - sin(-omega*TT_sou[kt2D(ix-hx+pad-1, iz-hz+pad+1)]))/(2*dx) - (sin(-omega*TT_sou[kt2D(ix-hx+pad+1, iz-hz+pad-1)]) - sin(-omega*TT_sou[kt2D(ix-hx+pad-1, iz-hz+pad-1)]))/(2*dx))/(2*dz); 
                                vxxri =  cdata[2*iw+1]*(sin(-omega*TT_rec[kt2D(ix-hx+pad+1, iz-hz+pad)]) - sin(-omega*TT_rec[kt2D(ix-hx+pad-1, iz-hz+pad)]))/(2*dx);
                                vxzri = cdata[2*iw+1]*(sin(-omega*TT_rec[kt2D(ix-hx+pad, iz-hz+pad+1)]) - sin(-omega*TT_rec[kt2D(ix-hx+pad, iz-hz+pad-1)]))/(2*dz);
                                imagedata[ki2D(ix,iz,ihx,ihz)] +=  (-2.0*CMULR(vzzsr,vzzsi,vxxrr,vxxri) + CMULR(vxzsr,vxzsi,vxzrr,vxzri));
                            }
                        }
                    }	
                }
            }
        }
    }
}

template<typename T>
KdmigElastic2D<T>::~KdmigElastic2D() {
    // Nothing here
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Kdmig<float>;
template class Kdmig<double>;
template class KdmigAcoustic2D<float>;
template class KdmigAcoustic2D<double>;
template class KdmigElastic2D<float>;
template class KdmigElastic2D<double>;
}
