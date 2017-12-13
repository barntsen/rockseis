// Include statements
#include "mva.h"

namespace rockseis {

// =============== ABSTRACT MVA CLASS =============== //
template<typename T>
Mva<T>::Mva() {
	order = 4;
    snapinc=1;
    snapmethod = FULL;
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

// =============== ELASTIC 2D MVA CLASS =============== //

template<typename T>
MvaElastic2D<T>::MvaElastic2D(){
    sourceset = false;
    dataUxset = false;
    dataUzset = false;
    modelset = false;
    pimageset = false;
    simageset = false;
}

template<typename T>
MvaElastic2D<T>::MvaElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataUx, std::shared_ptr<Data2D<T>> _dataUz, int order, int snapinc):Mva<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUz = _dataUz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUzset = true;
    pimageset = false;
    simageset = false;
}

template<typename T>
void MvaElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads, T* wrx, T* wrz, int padr, T* Vp, T* Vs, T* Rho)
{
	int ix, iz, ihx, ihz;
	T *pimagedata = NULL; 
	T *simagedata = NULL;
	T msxx, mszz, msxz, mrxx, mrzz, mrxz;
	T C33_minus, C33_plus;
	T C44_minus, C44_plus;
	int nhx; 
	int nhz;
	int nx;
	T dx;
	T dz;
	int nz;
	int hx, hz;

	if(pimageset){
		if(!pimage->getAllocated()){
			pimage->allocateImage();
		}
		pimagedata = pimage->getImagedata();
	}
	if(simageset){
		if(!simage->getAllocated()){
			simage->allocateImage();
		}
		simagedata = simage->getImagedata();
	}
	// Getting sizes
	if(pimageset) {
		nhx = pimage->getNhx();
		nhz = pimage->getNhz();
		nx = pimage->getNx();
		nz = pimage->getNz();
		dx = pimage->getDx(); 
		dz = pimage->getDz(); 
	}else{
		nhx = simage->getNhx();
		nhz = simage->getNhz();
		nx = simage->getNx();
		nz = simage->getNz();
		dx = simage->getDx(); 
		dz = simage->getDz(); 
	}

	int nxs = nx+2*pads;
	int nxr = nx+2*padr;

	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihz=0; ihz<nhz; ihz++){
			hz= -(nhz-1)/2 + ihz;
			for (ix=0; ix<nx; ix++){
				if( ((ix-hx) >= 1) && ((ix-hx) < nx-1) && ((ix+hx) >= 1) && ((ix+hx) < nx-1))
				{
					for (iz=0; iz<nz; iz++){
						if( ((iz-hz) >= 1) && ((iz-hz) < nz-1) && ((iz+hz) >= 1) && ((iz+hz) < nz-1))
						{
							C33_minus = Rho[km2D(ix-hx, iz-hz)]*Vp[km2D(ix-hx, iz-hz)]*Vp[km2D(ix-hx, iz-hz)];
							C33_plus = Rho[km2D(ix+hx, iz+hz)]*Vp[km2D(ix+hx, iz+hz)]*Vp[km2D(ix+hx, iz+hz)];
							C44_minus = Rho[km2D(ix-hx, iz-hz)]*Vs[km2D(ix-hx, iz-hz)]*Vs[km2D(ix-hx, iz-hz)];
							C44_plus = Rho[km2D(ix+hx, iz+hz)]*Vs[km2D(ix+hx, iz+hz)]*Vs[km2D(ix+hx, iz+hz)];

							msxx = (wsx[ks2D(ix-hx+pads, iz-hz+pads)] - wsx[ks2D(ix-hx+pads-1, iz-hz+pads)])/dx;
							mszz = (wsz[ks2D(ix-hx+pads, iz-hz+pads)] - wsz[ks2D(ix-hx+pads, iz-hz+pads-1)])/dz;
							mrxx = (wrx[kr2D(ix+hx+padr, iz+hz+padr)] - wrx[kr2D(ix+hx+padr-1, iz+hz+padr)])/dx;
							mrzz = (wrz[kr2D(ix+hx+padr, iz+hz+padr)] - wrz[kr2D(ix+hx+padr, iz+hz+padr-1)])/dz;

							if(pimageset){
								pimagedata[ki2D(ix,iz,ihx,ihz)] += C33_minus*C33_plus*(msxx + mszz) * (mrxx + mrzz);
							}

							if(simageset){
								msxz = 0.5*(wsx[ks2D(ix-hx+pads, iz-hz+pads+1)] - wsx[ks2D(ix-hx+pads, iz-hz+pads)])/dz;
								msxz += 0.5*(wsx[ks2D(ix-hx+pads-1, iz-hz+pads)] - wsx[ks2D(ix-hx+pads-1, iz-hz+pads-1)])/dz;
								msxz += 0.5*(wsz[ks2D(ix-hx+pads+1, iz-hz+pads)] - wsz[ks2D(ix-hx+pads, iz-hz+pads)])/dx;
								msxz += 0.5*(wsz[ks2D(ix-hx+pads, iz-hz+pads-1)] - wsz[ks2D(ix-hx+pads-1, iz-hz+pads-1)])/dx;

								mrxz = 0.5*(wrx[kr2D(ix+hx+padr, iz+hz+padr+1)] - wrx[kr2D(ix+hx+padr, iz+hz+padr)])/dz;
								mrxz += 0.5*(wrx[kr2D(ix+hx+padr-1, iz+hz+padr)] - wrx[kr2D(ix+hx+padr-1, iz+hz+padr-1)])/dz;
								mrxz += 0.5*(wrz[kr2D(ix+hx+padr+1, iz+hz+padr)] - wrz[kr2D(ix+hx+padr, iz+hz+padr)])/dx;
								mrxz += 0.5*(wrz[kr2D(ix+hx+padr, iz+hz+padr-1)] - wrz[kr2D(ix+hx+padr-1, iz+hz+padr-1)])/dx;
								simagedata[ki2D(ix,iz,ihx,ihz)] += C44_minus*C44_plus*(-2.0*msxx*mrzz + -2.0*mszz*mrxx + msxz*mrxz);
							}
						}
					}	
				}
			}
		}
	}
}
/*
template<typename T>
void MvaElastic2D<T>::adjSource(T *wx, T *wz, int pads, std::shared_ptr<WavesElastic2D<T>> waves_adj, T* Vp, T* Vs, T* Rho, bool sign)
{
	int ix, iz, ihx, ihz;
	T *pimagedata = NULL; 
	T *simagedata = NULL;
	T mxx, mzz;
	T C33_minus, C33_plus;
	T C44_minus, C44_plus;
	int nhx; 
	int nhz;
	int nx;
	T dx;
	T dz;
	int nz;
	int hx, hz;

	if(pimageset){
		if(!pimage->getAllocated()){
			pimage->allocateImage();
		}
		pimagedata = pimage->getImagedata();
	}
	if(simageset){
		if(!simage->getAllocated()){
			simage->allocateImage();
		}
		simagedata = simage->getImagedata();
	}
	// Getting sizes
	if(pimageset) {
		nhx = pimage->getNhx();
		nhz = pimage->getNhz();
		nx = pimage->getNx();
		nz = pimage->getNz();
		dx = pimage->getDx(); 
		dz = pimage->getDz(); 
	}else{
		nhx = simage->getNhx();
		nhz = simage->getNhz();
		nx = simage->getNx();
		nz = simage->getNz();
		dx = simage->getDx(); 
		dz = simage->getDz(); 
	}

	int nxs = nx+2*pads;
	int nxr = nx+2*padr;

	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihz=0; ihz<nhz; ihz++){
			hz= -(nhz-1)/2 + ihz;
			for (ix=0; ix<nx; ix++){
				if( ((ix-hx) >= 1) && ((ix-hx) < nx-1) && ((ix+hx) >= 1) && ((ix+hx) < nx-1))
				{
					for (iz=0; iz<nz; iz++){
						if( ((iz-hz) >= 1) && ((iz-hz) < nz-1) && ((iz+hz) >= 1) && ((iz+hz) < nz-1))
						{
							C33_minus = Rho[km2D(ix-hx, iz-hz)]*Vp[km2D(ix-hx, iz-hz)]*Vp[km2D(ix-hx, iz-hz)];
							C33_plus = Rho[km2D(ix+hx, iz+hz)]*Vp[km2D(ix+hx, iz+hz)]*Vp[km2D(ix+hx, iz+hz)];
							C44_minus = Rho[km2D(ix-hx, iz-hz)]*Vs[km2D(ix-hx, iz-hz)]*Vs[km2D(ix-hx, iz-hz)];
							C44_plus = Rho[km2D(ix+hx, iz+hz)]*Vs[km2D(ix+hx, iz+hz)]*Vs[km2D(ix+hx, iz+hz)];

							msxx = (wsx[ks2D(ix-hx+pads, iz-hz+pads)] - wsx[ks2D(ix-hx+pads-1, iz-hz+pads)])/dx;
							mszz = (wsz[ks2D(ix-hx+pads, iz-hz+pads)] - wsz[ks2D(ix-hx+pads, iz-hz+pads-1)])/dz;
							mrxx = (wrx[kr2D(ix+hx+padr, iz+hz+padr)] - wrx[kr2D(ix+hx+padr-1, iz+hz+padr)])/dx;
							mrzz = (wrz[kr2D(ix+hx+padr, iz+hz+padr)] - wrz[kr2D(ix+hx+padr, iz+hz+padr-1)])/dz;

							if(pimageset){
								pimagedata[ki2D(ix,iz,ihx,ihz)] += C33_minus*C33_plus*(msxx + mszz) * (mrxx + mrzz);
							}

							if(simageset){
								msxz = 0.5*(wsx[ks2D(ix-hx+pads, iz-hz+pads+1)] - wsx[ks2D(ix-hx+pads, iz-hz+pads)])/dz;
								msxz += 0.5*(wsx[ks2D(ix-hx+pads-1, iz-hz+pads)] - wsx[ks2D(ix-hx+pads-1, iz-hz+pads-1)])/dz;
								msxz += 0.5*(wsz[ks2D(ix-hx+pads+1, iz-hz+pads)] - wsz[ks2D(ix-hx+pads, iz-hz+pads)])/dx;
								msxz += 0.5*(wsz[ks2D(ix-hx+pads, iz-hz+pads-1)] - wsz[ks2D(ix-hx+pads-1, iz-hz+pads-1)])/dx;

								mrxz = 0.5*(wrx[kr2D(ix+hx+padr, iz+hz+padr+1)] - wrx[kr2D(ix+hx+padr, iz+hz+padr)])/dz;
								mrxz += 0.5*(wrx[kr2D(ix+hx+padr-1, iz+hz+padr)] - wrx[kr2D(ix+hx+padr-1, iz+hz+padr-1)])/dz;
								mrxz += 0.5*(wrz[kr2D(ix+hx+padr+1, iz+hz+padr)] - wrz[kr2D(ix+hx+padr, iz+hz+padr)])/dx;
								mrxz += 0.5*(wrz[kr2D(ix+hx+padr, iz+hz+padr-1)] - wrz[kr2D(ix+hx+padr-1, iz+hz+padr-1)])/dx;
								simagedata[ki2D(ix,iz,ihx,ihz)] += C44_minus*C44_plus*(-2.0*msxx*mrzz + -2.0*mszz*mrxx + msxz*mrxz);
							}
						}
					}	
				}
			}
		}
	}
}

*/

template<typename T>
int MvaElastic2D<T>::runRtm(){
	int result = MVA_ERR;
	if(!pimageset && !simageset) {
		rs_warning("MvaElastic2D::runRtm: No image set");
		return result;
	}
	int nt;
	float dt;
	float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Uxsnap;
     Uxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Uxsnap->openSnap(this->getSnapfile() + "-ux", 'w'); // Create a new snapshot file
     Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field

     std::shared_ptr<Snapshot2D<T>> Uzsnap;
     Uzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Uzsnap->openSnap(this->getSnapfile() + "-uz", 'w'); // Create a new snapshot file
     Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field

     this->writeLog("Running 2D Elastic reverse-time migration with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field
        Uxsnap->writeSnap(it);

        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field
        Uzsnap->writeSnap(it);

    	// Time stepping
    	waves->forwardstepStress(model, der);
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Uxsnap->closeSnap();
    Uzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic2D_DS<T>>(model, nt, dt, ot);

    // Create image
    if(this->pimageset) pimage->allocateImage();
    if(this->simageset) simage->allocateImage();

    Uxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
    Uxsnap->allocSnap(0);

    Uzsnap->openSnap(this->getSnapfile() + "-uz", 'r');
    Uzsnap->allocSnap(0);

    // Get models for scaling
    T *Vp, *Vs, *Rho;
    Vp = model->getVp();
    Vs = model->getVs();
    Rho = model->getR();

    this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping 
    	waves->forwardstepStress(model, der);
    	waves->forwardstepDisplacement(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataUx, GMAP, (nt - 1 - it));
    	waves->insertSource(model, dataUz, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Uxsnap->readSnap(nt - 1 - it);
        Uzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Uxsnap->getEnddiff()) % Uxsnap->getSnapinc()) == 0){
            T *Uxr = waves->getUx1();
            T *Uzr = waves->getUz1();
            crossCorr(Uxsnap->getData(0), Uzsnap->getData(0), 0, Uxr, Uzr, waves->getLpml(), Vp, Vs, Rho);
        }
        
        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Uxsnap->removeSnap();
	Uzsnap->removeSnap();

    result=MVA_OK;
    return result;
}

template<typename T>
int MvaElastic2D<T>::runRtm_optimal(){
     int result = MVA_ERR;
     if(!pimageset && !simageset) {
         rs_warning("MvaElastic2D::runRtm: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), 1, waves_fw->getNz_pml(), waves_fw->getDx(), 1.0, waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
     pimage->allocateImage();
     simage->allocateImage();

     // Get models for scaling
     T *Vp, *Vs, *Rho;
     Vp = model->getVp();
     Vs = model->getVs();
     Rho = model->getR();

     this->writeLog("Running 2D Elastic reverse-time migration with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal->getCapo();
        whatodo = optimal->revolve();
        capo = optimal->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping
                waves_fw->forwardstepStress(model, der);
                waves_fw->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Roll pointers
                waves_fw->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }

        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepStress(model, der);
            waves_fw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertSource(model, dataUx, GMAP, capo);
            waves_bw->insertSource(model, dataUz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getUx1();
            T *wsz = waves_fw->getUz1();
            T *wrx = waves_bw->getUx1();
            T *wrz = waves_bw->getUz1();

            crossCorr(wsx, wsz, waves_fw->getLpml(), wrx, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

            // Roll pointers
            waves_fw->roll();
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);
      
            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping
            waves_bw->forwardstepStress(model, der);
            waves_bw->forwardstepDisplacement(model, der);

            // Inserting data
            waves_bw->insertSource(model, dataUx, GMAP, capo);
            waves_bw->insertSource(model, dataUz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getUx1();
            T *wsz = waves_fw->getUz1();
            T *wrx = waves_bw->getUx1();
            T *wrz = waves_bw->getUz1();
            crossCorr(wsx, wsz, waves_fw->getLpml(), wrx, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

            // Roll pointers
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal->writeCheck(waves_fw);
        }
        if (whatodo == restore)
        {
            optimal->readCheck(waves_fw);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));

	//Remove snapshot file
	optimal->removeCheck();

    result=MVA_OK;
    return result;
}

template<typename T>
MvaElastic2D<T>::~MvaElastic2D() {
    // Nothing here
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Mva<float>;
template class Mva<double>;

template class MvaElastic2D<float>;
template class MvaElastic2D<double>;

}
