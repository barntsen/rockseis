// Include statements
#include "mva.h"

namespace rockseis {

// =============== ABSTRACT MVA CLASS =============== //
template<typename T>
Mva<T>::Mva() {
	order = 4;
    snapinc=1;
    snapmethod = FULL;
    misfit_type = SI;
    ncheck = 0;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
}

template<typename T>
Mva<T>::Mva(int _order, int _snapinc) {
    if(_order > 1 && _order < 9)
    {
        order = _order;
    }else{
        order = 4;
    }
    if(_snapinc > 0)
    {
        snapinc = _snapinc;
    }else{
        snapinc=1;
    }

    snapmethod = FULL;
    misfit_type = SI;
    ncheck = 0;

	prog.previous = 0;
	prog.current = 0;
    prog.persec = 1;
}

template<typename T>
bool Mva<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return MVA_ERR;
	}else{
		Flog.close();
		return MVA_OK;
	}
}

template<typename T>
void Mva<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Mva<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Mva<T>::writeProgressbar(int x, int n, int r, int w){
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
void Mva<T>::writeProgress(int x, int n, int r, int w){
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
Mva<T>::~Mva() {
    // Nothing here
}

// =============== ACOUSTIC 2D MVA CLASS =============== //

template<typename T>
MvaAcoustic2D<T>::MvaAcoustic2D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    pimageset = false;
}

template<typename T>
MvaAcoustic2D<T>::MvaAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<Image2D<T>> _pimage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataP, int order, int snapinc):Mva<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    pimage = _pimage;
    sourceset = true;
    dataPset = true;
    modelset = true;
    pimageset = true;
}

template<typename T>
bool MvaAcoustic2D<T>::checkStability(){
    T *Vp = model->getVp();
    // Find maximum Vp
    T Vpmax;
    Vpmax=Vp[0];
    size_t n=model->getNx()*model->getNz();
    for(size_t i=1; i<n; i++){
        if(Vp[i] > Vpmax){
            Vpmax = Vp[i];
        }
    }

    T dx = model->getDx();
    T dz = model->getDz();
    T dt = source->getDt();
    T dt_stab;
    dt_stab = 2.0/(3.1415*sqrt((1.0/(dx*dx))+(1/(dz*dz)))*Vpmax); 
    if(dt < dt_stab){
        return true;
    }else{
        rs_warning("Modeling time interval exceeds maximum stable number of: ", std::to_string(dt_stab));
        return false;
    }
}

template<typename T>
void MvaAcoustic2D<T>::calcAdjointsource(T *adjsrc_fw, T* wsp, int pads, T *adjsrc_bw, T* wrp, int padr)
{
	if(!pimage->getAllocated()) rs_error("MvaAcoustic2D<T>::MvaAcoustic2D<T>::calcAdjointsource: pimage is not allocated.");
	int ix, iz, ihx, ihz;
	T *imagedata = pimage->getImagedata();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nxs = nx+2*pads;
	int nxr = nx+2*padr;
	int nz = pimage->getNz();
	int hx, hz;
    //Reset arrays
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            adjsrc_fw[km2D(ix,iz)] = 0.0;
            adjsrc_bw[km2D(ix,iz)] = 0.0;
        }
    }
    // Calculate integral
    for (ihx=0; ihx<nhx; ihx++){
        hx= -(nhx-1)/2 + ihx;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            for (ix=0; ix<nx; ix++){
                if( ((ix-2*hx) >= 0) && ((ix-2*hx) < nx) && ((ix+2*hx) >= 0) && ((ix+2*hx) < nx))
                {
                    for (iz=0; iz<nz; iz++){
                        if( ((iz-2*hz) >= 0) && ((iz-2*hz) < nz) && ((iz+2*hz) >= 0) && ((iz+2*hz) < nz))
                        {
							adjsrc_fw[km2D(ix,iz)] += imagedata[ki2D(ix+hx,iz+hz,ihx,ihz)]*wsp[kr2D(ix+2*hx+padr, iz+2*hz+padr)];
							adjsrc_bw[km2D(ix,iz)] += imagedata[ki2D(ix-hx,iz-hz,ihx,ihz)]*wrp[ks2D(ix-2*hx+pads, iz-2*hz+pads)];

                        }

                    }	
                }
            }
        }
    }
}

template<typename T>
void MvaAcoustic2D<T>::insertAdjointsource(std::shared_ptr<WavesAcoustic2D<T>> waves_fw, T* adjsrc_fw, std::shared_ptr<WavesAcoustic2D<T>> waves_bw, T* adjsrc_bw, T* L)
{
    int ix, iz;
    int nx = pimage->getNx();
    int nz = pimage->getNz();
    T dt = waves_fw->getDt();
    int padw = waves_fw->getLpml();
    int nxw = nx+2*padw;
    T* ws = waves_fw->getP2();
    T* wr = waves_bw->getP2();
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            ws[kw2D(ix+padw,iz+padw)] += dt*dt*L[kw2D(ix+padw,iz+padw)]*adjsrc_fw[km2D(ix, iz)];
            wr[kw2D(ix+padw,iz+padw)] += dt*dt*L[kw2D(ix+padw,iz+padw)]*adjsrc_bw[km2D(ix, iz)];

        }

    }	
}

template<typename T>
void MvaAcoustic2D<T>::modifyImage()
{
	if(!pimage->getAllocated()) rs_error("MvaAcoustic2D<T>::modifyImage(): pimage is not allocated.");
	int ix, iz, ihx, ihz;
	T *imagedata = pimage->getImagedata();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nz = pimage->getNz();
	int hx, hz;
    T G1,G2;
    // Allocate work array
    T *wrk = (T *) calloc(nz, sizeof(T));
	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
        G1 = GAUSS(hx, 0.25*nhx);
		for (ihz=0; ihz<nhz; ihz++){
			hz= -(nhz-1)/2 + ihz;
            G2 = G1*GAUSS(hz, 0.25*nhz);
			for (ix=0; ix<nx; ix++){
				{
					for (iz=1; iz<nz-1; iz++){
						wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
					}	
					for (iz=0; iz<nz; iz++){
						imagedata[ki2D(ix,iz,ihx,ihz)] = G2*wrk[iz]; 
					}	
				}
			}
		}
	}
    // Free work array
    free(wrk);
}

template<typename T>
void MvaAcoustic2D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T* Vp, T* Rho, T* adjsrc)
{
    if(!vpgradset) rs_error("MvaAcoustic2D:crossCorr: No gradient set in mva class");
    if(!vpgrad->getAllocated()) vpgrad->allocateImage();

    int ix, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T mrxx, mrzz;
    T L;
    T vpscale;
    int nx = vpgrad->getNx();
    int nz = vpgrad->getNz();
    int dx = vpgrad->getDx();
    int dz = vpgrad->getDz();
    int nxs = nx+2*pads;
    int nxr = nx+2*padr;
    for (ix=1; ix<nx-1; ix++){
        {
            for (iz=1; iz<nz-1; iz++){
                L = Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)];
                vpscale = -2.0/(Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)]);
                mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
                mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
                vpgraddata[km2D(ix,iz)] -= vpscale*wsp[ks2D(ix+pads, iz+pads)]*L*(mrxx + mrzz + adjsrc[km2D(ix,iz)]);
            }
        }	
    }
}

template<typename T>
int MvaAcoustic2D<T>::run(){
     int result = MVA_ERR;
     int nt;
     float dt;
	 float ot;
     T *wrp;
     T *wrx;
     T *wrz;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("MvaAcoustic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves_adj_fw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_adj_bw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_adj_fw->getNx_pml(), 1, waves_adj_fw->getNz_pml(), waves_adj_fw->getDx(), 1.0, waves_adj_fw->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Fwsnap;
     Fwsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Fwsnap->openSnap(this->getSnapfile(), 'r'); // Open a snapshot file for reading
     Fwsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Bwsnap;
     Bwsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Bwsnap->openSnap(this->getSnapfile() + "-bw", 'r'); // Open a snapshot file for reading
     Bwsnap->allocSnap(0);

    // Create gradient array
    vpgrad->allocateImage();

    // Modify local image
    //this->modifyImage();

    // Allocate memory for adjoint sources 
    T* adjsrc_fw, *adjsrc_bw;
    adjsrc_fw = (T *) calloc(waves_adj_fw->getNx()*waves_adj_fw->getNz(), sizeof(T));
    adjsrc_bw = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Acoustic RTMVA gradient with full checkpointing.");
     this->writeLog("Doing time loop.");
     int imageit=0;
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves_adj_fw->forwardstepAcceleration(model, der);
    	waves_adj_fw->forwardstepStress(model, der);

    	waves_adj_bw->forwardstepAcceleration(model, der);
    	waves_adj_bw->forwardstepStress(model, der);
    
        //Read snapshots (interpolation needed here if snapinc not 1!)
        Fwsnap->setSnapit(imageit);
        Fwsnap->readSnap(it);
        Bwsnap->setSnapit(imageit);
        Bwsnap->readSnap(it);
    	
        // Inserting adjoint sources
    	this->calcAdjointsource(adjsrc_fw, Bwsnap->getData(0), 0, adjsrc_bw, Fwsnap->getData(0), 0);
    	this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model->getL());

        //Read snapshots 
        Fwsnap->setSnapit(Fwsnap->getSnapnt() - 1 - imageit);
        Fwsnap->readSnap(it);
        Bwsnap->setSnapit(Bwsnap->getSnapnt() - 1 - imageit);
        Bwsnap->readSnap(it);
        // Do Crosscorrelation
        if((((nt - 1 - it)-Fwsnap->getEnddiff()) % Fwsnap->getSnapinc()) == 0){
            wrp = waves_adj_fw->getP1();
            wrx = waves_adj_fw->getAx(); 
            wrz = waves_adj_fw->getAz(); 
            crossCorr(Fwsnap->getData(0), 0, wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);
            wrp = waves_adj_bw->getP1();
            wrx = waves_adj_bw->getAx(); 
            wrz = waves_adj_bw->getAz(); 
            crossCorr(Bwsnap->getData(0), 0, wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);
            imageit++;
        }

    	// Roll the pointers P1 and P2
    	waves_adj_fw->roll();
    	waves_adj_bw->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of time loop
    
	//Remove snapshot file
	Fwsnap->removeSnap();
	Bwsnap->removeSnap();

    result=MVA_OK;
    return result;
}

template<typename T>
int MvaAcoustic2D<T>::run_optimal(){
     int result = MVA_ERR;
     int nt;
     float dt;
	 float ot;
     T *wsp;
     T *wrp;
     T *wrx;
     T *wrz;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("MvaAcoustic2D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves_fw1 (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_fw2 (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_bw1 (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_bw2 (new WavesAcoustic2D<T>(model, nt, dt, ot));

     std::shared_ptr<WavesAcoustic2D<T>> waves_adj_fw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_adj_bw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw1->getNx_pml(), 1, waves_fw1->getNz_pml(), waves_fw1->getDx(), 1.0, waves_fw1->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal_fw (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     std::shared_ptr<Revolve<T>> optimal_bw (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint files
     optimal_fw->openCheck(this->getSnapfile() + "-fw", waves_fw1, 'w');
     optimal_bw->openCheck(this->getSnapfile() + "-bw", waves_bw1, 'w');

     // Create image
     vpgrad->allocateImage();

    // Modify local image
    this->modifyImage();

    // Allocate memory for adjoint sources 
    T* adjsrc_fw, *adjsrc_bw;
    adjsrc_fw = (T *) calloc(waves_adj_fw->getNx()*waves_adj_fw->getNz(), sizeof(T));
    adjsrc_bw = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Acoustic RTMVA gradient with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal_fw->getCapo();
        whatodo = optimal_fw->revolve();
        whatodo = optimal_bw->revolve();
        capo = optimal_fw->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping
                waves_fw1->forwardstepAcceleration(model, der);
                waves_fw1->forwardstepStress(model, der);

                waves_bw1->forwardstepAcceleration(model, der);
                waves_bw1->forwardstepStress(model, der);

                // Inserting source 
                waves_fw1->insertSource(model, source, SMAP, it);
                waves_bw1->insertSource(model, source, GMAP, nt-1-it);

                // Roll the pointers P1 and P2
                waves_fw1->roll();
                waves_bw1->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw1->forwardstepAcceleration(model, der);
            waves_fw1->forwardstepStress(model, der);
            waves_bw1->forwardstepAcceleration(model, der);
            waves_bw1->forwardstepStress(model, der);

            waves_fw2->forwardstepAcceleration(model, der);
            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepAcceleration(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepAcceleration(model, der);
            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepAcceleration(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw1->insertSource(model, source, SMAP, capo);
            waves_bw1->insertSource(model, dataP, GMAP, nt-1-capo);
            waves_fw2->insertSource(model, source, SMAP, nt-1-capo);
            waves_bw2->insertSource(model, dataP, GMAP, capo);

            // Inserting adjoint sources
            wsp = waves_fw2->getP1();
	    wrp = waves_bw2->getP1();
	    this->calcAdjointsource(adjsrc_fw, wrp, waves_adj_fw->getLpml(), adjsrc_bw, wsp, waves_adj_bw->getLpml());
	    this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model->getL());

            /* Do Crosscorrelation */
            wsp = waves_fw1->getP1();
            wrp = waves_adj_fw->getP1();
            wrx = waves_adj_fw->getAx(); 
            wrz = waves_adj_fw->getAz(); 
            crossCorr(wsp, waves_fw1->getLpml(), wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);
            wsp = waves_bw1->getP1();
            wrp = waves_adj_bw->getP1();
            wrx = waves_adj_bw->getAx(); 
            wrz = waves_adj_bw->getAz(); 
            crossCorr(wsp, waves_bw1->getLpml(), wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);

            // Roll the pointers P1 and P2
            waves_fw1->roll();
            waves_bw1->roll();
            waves_fw2->roll();
            waves_bw2->roll();
            waves_adj_fw->roll();
            waves_adj_bw->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);

            //Close checkpoint file for w and reopen for rw
            optimal_fw->closeCheck();
            optimal_fw->openCheck(this->getSnapfile() + "-fw", waves_fw1, 'a');

            optimal_bw->closeCheck();
            optimal_bw->openCheck(this->getSnapfile() + "-bw", waves_bw1, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping
            waves_fw2->forwardstepAcceleration(model, der);
            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepAcceleration(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepAcceleration(model, der);
            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepAcceleration(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting source
            waves_fw2->insertSource(model, source, SMAP, nt-1-capo);
            waves_bw2->insertSource(model, dataP, GMAP, capo);

            // Inserting adjoint sources
            T *wsp = waves_fw2->getP1();
            T *wrp = waves_bw2->getP1();
	    this->calcAdjointsource(adjsrc_fw, wrp, waves_adj_fw->getLpml(), adjsrc_bw, wsp, waves_adj_bw->getLpml());
	    this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model->getL());

            /* Do Crosscorrelation */
            wsp = waves_fw1->getP1();
            wrp = waves_adj_fw->getP1();
            T* wrx = waves_adj_fw->getAx(); 
            T* wrz = waves_adj_fw->getAz(); 
            crossCorr(wsp, waves_fw1->getLpml(), wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);

            wsp = waves_bw1->getP1();
            wrp = waves_adj_bw->getP1();
            wrx = waves_adj_bw->getAx(); 
            wrz = waves_adj_bw->getAz(); 
            crossCorr(wsp, waves_bw1->getLpml(), wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);

            // Roll the pointers P1 and P2
            waves_fw2->roll();
            waves_bw2->roll();
    	    waves_adj_fw->roll();
    	    waves_adj_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal_fw->writeCheck(waves_fw1);
            optimal_bw->writeCheck(waves_bw1);
        }
        if (whatodo == restore)
        {
            optimal_fw->readCheck(waves_fw1);
            optimal_bw->readCheck(waves_bw1);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));


	//Remove snapshot file
	optimal_fw->removeCheck();
	optimal_bw->removeCheck();

    result=MVA_OK;
    return result;
}

template<typename T>
MvaAcoustic2D<T>::~MvaAcoustic2D() {
    // Nothing here
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Mva<float>;
template class Mva<double>;

template class MvaAcoustic2D<float>;
template class MvaAcoustic2D<double>;

}
