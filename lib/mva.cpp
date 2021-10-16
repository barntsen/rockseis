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
    vpgradset = false;
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
    vpgradset = false;
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
    T* ws = waves_fw->getP();
    T* wr = waves_bw->getP();
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            ws[kw2D(ix+padw,iz+padw)] += dt*L[kw2D(ix+padw,iz+padw)]*adjsrc_fw[km2D(ix, iz)];
            wr[kw2D(ix+padw,iz+padw)] += dt*L[kw2D(ix+padw,iz+padw)]*adjsrc_bw[km2D(ix, iz)];

        }

    }	
}

template<typename T>
void MvaAcoustic2D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T* Vp, T* Rho, T* adjsrc)
{
    if(!vpgradset) rs_error("MvaAcoustic2D:crossCorr: No gradient set in mva class");
    if(!vpgrad->getAllocated()) vpgrad->allocateImage();

    int ix, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T mrxx, mrzz;
    T vpscale;
    int nx = vpgrad->getNx();
    int nz = vpgrad->getNz();
    T dx = vpgrad->getDx();
    T dz = vpgrad->getDz();
    int nxs = nx+2*pads;
    int nxr = nx+2*padr;
    for (ix=1; ix<nx-1; ix++){
        {
            for (iz=1; iz<nz-1; iz++){
                vpscale = -2.0/(Vp[km2D(ix, iz)]);
                mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
                mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
                vpgraddata[km2D(ix,iz)] -= vpscale*wsp[ks2D(ix+pads, iz+pads)]*(mrxx + mrzz + adjsrc[km2D(ix,iz)]);
            }
        }	
    }
}

template<typename T>
int MvaAcoustic2D<T>::run(){
     int result = MVA_ERR;
     int nt;
     T dt;
	 T ot;
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
    	waves_adj_fw->forwardstepVelocity(model, der);
    	waves_adj_fw->forwardstepStress(model, der);

    	waves_adj_bw->forwardstepVelocity(model, der);
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
            wrp = waves_adj_fw->getP();
            wrx = waves_adj_fw->getVx(); 
            wrz = waves_adj_fw->getVz(); 
            crossCorr(Fwsnap->getData(0), 0, wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);
            wrp = waves_adj_bw->getP();
            wrx = waves_adj_bw->getVx(); 
            wrz = waves_adj_bw->getVz(); 
            crossCorr(Bwsnap->getData(0), 0, wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);
            imageit++;
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of time loop
    
	//Remove snapshot file
	Fwsnap->removeSnap();
	Bwsnap->removeSnap();

    // Free arrays
    free(adjsrc_fw);
    free(adjsrc_bw);

    result=MVA_OK;
    return result;
}

template<typename T>
int MvaAcoustic2D<T>::run_optimal(){
     int result = MVA_ERR;
     int nt;
     T dt;
	 T ot;
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
                waves_fw1->forwardstepVelocity(model, der);
                waves_fw1->forwardstepStress(model, der);

                waves_bw1->forwardstepVelocity(model, der);
                waves_bw1->forwardstepStress(model, der);

                // Inserting source 
                waves_fw1->insertSource(model, source, SMAP, it);
                waves_bw1->insertSource(model, dataP, GMAP, nt-1-it);

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw1->forwardstepVelocity(model, der);
            waves_fw1->forwardstepStress(model, der);
            waves_bw1->forwardstepVelocity(model, der);
            waves_bw1->forwardstepStress(model, der);

            waves_fw2->forwardstepVelocity(model, der);
            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepVelocity(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepVelocity(model, der);
            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepVelocity(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw1->insertSource(model, source, SMAP, capo);
            waves_bw1->insertSource(model, dataP, GMAP, nt-1-capo);
            waves_fw2->insertSource(model, source, SMAP, nt-1-capo);
            waves_bw2->insertSource(model, dataP, GMAP, capo);

            // Inserting adjoint sources
            wsp = waves_fw2->getP();
            wrp = waves_bw2->getP();
            this->calcAdjointsource(adjsrc_fw, wrp, waves_adj_fw->getLpml(), adjsrc_bw, wsp, waves_adj_bw->getLpml());
            this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model->getL());

            /* Do Crosscorrelation */
            wsp = waves_fw1->getP();
            wrp = waves_adj_fw->getP();
            wrx = waves_adj_fw->getVx(); 
            wrz = waves_adj_fw->getVz(); 
            crossCorr(wsp, waves_fw1->getLpml(), wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);
            wsp = waves_bw1->getP();
            wrp = waves_adj_bw->getP();
            wrx = waves_adj_bw->getVx(); 
            wrz = waves_adj_bw->getVz(); 
            crossCorr(wsp, waves_bw1->getLpml(), wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);

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
            waves_fw2->forwardstepVelocity(model, der);
            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepVelocity(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepVelocity(model, der);
            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepVelocity(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting source
            waves_fw2->insertSource(model, source, SMAP, nt-1-capo);
            waves_bw2->insertSource(model, dataP, GMAP, capo);

            // Inserting adjoint sources
            T *wsp = waves_fw2->getP();
            T *wrp = waves_bw2->getP();
            this->calcAdjointsource(adjsrc_fw, wrp, waves_adj_fw->getLpml(), adjsrc_bw, wsp, waves_adj_bw->getLpml());
            this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model->getL());

            /* Do Crosscorrelation */
            wsp = waves_fw1->getP();
            wrp = waves_adj_fw->getP();
            T* wrx = waves_adj_fw->getVx(); 
            T* wrz = waves_adj_fw->getVz(); 
            crossCorr(wsp, waves_fw1->getLpml(), wrp, wrx, wrz, waves_adj_fw->getLpml(), model->getVp(), model->getR(), adjsrc_fw);

            wsp = waves_bw1->getP();
            wrp = waves_adj_bw->getP();
            wrx = waves_adj_bw->getVx(); 
            wrz = waves_adj_bw->getVz(); 
            crossCorr(wsp, waves_bw1->getLpml(), wrp, wrx, wrz, waves_adj_bw->getLpml(), model->getVp(), model->getR(), adjsrc_bw);

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

    // Free arrays
    free(adjsrc_fw);
    free(adjsrc_bw);

    result=MVA_OK;
    return result;
}

template<typename T>
MvaAcoustic2D<T>::~MvaAcoustic2D() {
    // Nothing here
}


// =============== ELASTIC 2D PPmva CLASS =============== //

template<typename T>
PPmvaElastic2D<T>::PPmvaElastic2D(){
    sourceset = false;
    dataUxset = false;
    dataUzset = false;
    modelset = false;
    pimageset = false;
    vpgradset = false;
}

template<typename T>
PPmvaElastic2D<T>::PPmvaElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Image2D<T>> _pimage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataUx, std::shared_ptr<Data2D<T>> _dataUz, int order, int snapinc):Mva<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUz = _dataUz;
    model = _model;
    pimage = _pimage;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUzset = true;
    pimageset = false;
    vpgradset = false;
}

template<typename T>
bool PPmvaElastic2D<T>::checkStability(){
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
void PPmvaElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads, std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, T *adjsrc)
{
    if(!vpgradset) rs_error("PPmvaElastic2D<T>::crossCorr: No gradient set for computation.");
	int ix, iz;
    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getUx1();
    T* wrz = waves_bw->getUz1();

	T *vpgraddata = NULL; 
	T msxx=0, mszz=0, mrxx=0, mrzz=0;
	int nx;
	int nz;
	T dx;
	T dz;

    T* Vp = model->getVp();
    T* Rho = model->getR();
    T vpscale;

    if(!vpgrad->getAllocated()){
        vpgrad->allocateImage();
    }
    vpgraddata = vpgrad->getImagedata();

    // Getting sizes
    nx = waves_bw->getNx();
    nz = waves_bw->getNz();
    dx = waves_bw->getDx(); 
    dz = waves_bw->getDz(); 
    T adj;

	int nxs = nx+2*pads;
    int nxr = nx+2*padr;

    for (ix=1; ix<nx-1; ix++){
        for (iz=1; iz<nz-1; iz++){
            msxx = (wsx[ks2D(ix+pads, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads)])/dx;
            mszz = (wsz[ks2D(ix+pads, iz+pads)] - wsz[ks2D(ix+pads, iz+pads-1)])/dz;
            mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
            mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
            adj  = adjsrc[km2D(ix,iz)];

            vpscale = 2.0*Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)];
            vpgraddata[km2D(ix,iz)] += vpscale*(msxx + mszz) * (mrxx + mrzz + adj);

        }
    }	
}

template<typename T>
void PPmvaElastic2D<T>::calcAdjointsource(T *adjsrc_fw, T *wsx, T *wsz, int pads, T *adjsrc_bw, T* wrx, T* wrz, int padr, std::shared_ptr<ModelElastic2D<T>> model)
{
	if(!pimage->getAllocated()) rs_error("PPmvaElastic2D<T>::calcAdjointsource: pimage is not allocated.");
	int ix, iz, ihx, ihz;
	T *imagedata = pimage->getImagedata();
	T msxx, mszz, mrxx, mrzz;
	T C33_minus, C33_plus;
    T* Vp = model->getVp();
    T* Rho = model->getR();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nxs = nx+2*pads;
	int nxr = nx+2*padr;
	int nz = pimage->getNz();
    T dx = pimage->getDx(); 
    T dz = pimage->getDz(); 
	int hx, hz;


    float C;
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
        C = 1.0;
        if(((ihx == 0) || (ihx == nhx-1)) && (nhx > 1)) C*=0.5;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            if(((ihz == 0) || (ihz == nhz-1)) && (nhz > 1)) C*=0.5;
            for (ix=1; ix<nx-1; ix++){
                if( ((ix-2*hx) >= 0) && ((ix-2*hx) < nx) && ((ix+2*hx) >= 0) && ((ix+2*hx) < nx))
                {
                    for (iz=1; iz<nz-1; iz++){
                        if( ((iz-2*hz) >= 0) && ((iz-2*hz) < nz) && ((iz+2*hz) >= 0) && ((iz+2*hz) < nz))
                        {
							C33_minus = Rho[km2D(ix-2*hx, iz-2*hz)]*Vp[km2D(ix-2*hx, iz-2*hz)]*Vp[km2D(ix-2*hx, iz-2*hz)];
							C33_plus = Rho[km2D(ix+2*hx, iz+2*hz)]*Vp[km2D(ix+2*hx, iz+2*hz)]*Vp[km2D(ix+2*hx, iz+2*hz)];

							msxx = (wsx[ks2D(ix-2*hx+pads, iz-2*hz+pads)] - wsx[ks2D(ix-2*hx+pads-1, iz-2*hz+pads)])/dx;
							mszz = (wsz[ks2D(ix-2*hx+pads, iz-2*hz+pads)] - wsz[ks2D(ix-2*hx+pads, iz-2*hz+pads-1)])/dz;
							mrxx = (wrx[kr2D(ix+2*hx+padr, iz+2*hz+padr)] - wrx[kr2D(ix+2*hx+padr-1, iz+2*hz+padr)])/dx;
							mrzz = (wrz[kr2D(ix+2*hx+padr, iz+2*hz+padr)] - wrz[kr2D(ix+2*hx+padr, iz+2*hz+padr-1)])/dz;
							adjsrc_fw[km2D(ix,iz)] += C*imagedata[ki2D(ix+hx,iz+hz,ihx,ihz)]*C33_plus*(mrxx + mrzz);
							adjsrc_bw[km2D(ix,iz)] += C*imagedata[ki2D(ix-hx,iz-hz,ihx,ihz)]*C33_minus*(msxx + mszz);
                        }

                    }	
                }
            }
        }
    }
}

template<typename T>
void PPmvaElastic2D<T>::insertAdjointsource(std::shared_ptr<WavesElastic2D_DS<T>> waves_fw, T* adjsrc_fw, std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, T* adjsrc_bw, std::shared_ptr<ModelElastic2D<T>> model)
{
    int ix, iz;
    int nx = pimage->getNx();
    int nz = pimage->getNz();
    int padw = waves_fw->getLpml();
    int nxw = nx+2*padw;
	T C33;
    T* Vp = model->getVp();
    T* Rho = model->getR();
    T* wsx = waves_fw->getSxx();
    T* wsz = waves_fw->getSzz();

    T* wrx = waves_bw->getSxx();
    T* wrz = waves_bw->getSzz();
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            C33 = Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)];
            wsx[kw2D(ix+padw,iz+padw)] += C33*adjsrc_fw[km2D(ix, iz)];
            wsz[kw2D(ix+padw,iz+padw)] += C33*adjsrc_fw[km2D(ix, iz)];
            wrx[kw2D(ix+padw,iz+padw)] += C33*adjsrc_bw[km2D(ix, iz)];
            wrz[kw2D(ix+padw,iz+padw)] += C33*adjsrc_bw[km2D(ix, iz)];
        }
    }	
}

template<typename T>
int PPmvaElastic2D<T>::run(){
	int result = RTM_ERR;
	if(!vpgradset) {
		rs_warning("PPmvaElastic2D::run: No gradient set");
		return result;
	}
	int nt;
	T dt;
	T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("PPmvaElastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_fw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_adj_fw->getNx_pml(), 1, waves_adj_fw->getNz_pml(), waves_adj_fw->getDx(), 1.0, waves_adj_fw->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Fwuxsnap;
     Fwuxsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Fwuxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
     Fwuxsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Fwuzsnap;
     Fwuzsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Fwuzsnap->openSnap(this->getSnapfile() + "-uz", 'r');
     Fwuzsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Bwuxsnap;
     Bwuxsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Bwuxsnap->openSnap(this->getSnapfile() + "-ux-bw", 'r');
     Bwuxsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Bwuzsnap;
     Bwuzsnap = std::make_shared<Snapshot2D<T>>(waves_adj_fw, this->getSnapinc());
     Bwuzsnap->openSnap(this->getSnapfile() + "-uz-bw", 'r');
     Bwuzsnap->allocSnap(0);

    // Create gradient array
    if(this->vpgradset) vpgrad->allocateImage();

    // Allocate memory for adjoint sources 
    T* adjsrc_fw, *adjsrc_bw;
    adjsrc_fw = (T *) calloc(waves_adj_fw->getNx()*waves_adj_fw->getNz(), sizeof(T));
    adjsrc_bw = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Elastic PP RTMVA gradient with full checkpointing.");
     this->writeLog("Doing time loop.");
    // Loop over reverse time
    int imageit=0;
    for(int it=0; it < nt; it++)
    {
    	// Time stepping stress
    	waves_adj_fw->forwardstepStress(model, der);
    	waves_adj_bw->forwardstepStress(model, der);
    
        //Read snapshots (interpolation needed here if snapinc not 1!)
        Fwuxsnap->setSnapit(imageit);
        Fwuxsnap->readSnap(it);
        Fwuzsnap->setSnapit(imageit);
        Fwuzsnap->readSnap(it);

        Bwuxsnap->setSnapit(imageit);
        Bwuxsnap->readSnap(it);
        Bwuzsnap->setSnapit(imageit);
        Bwuzsnap->readSnap(it);

        // Inserting adjoint sources
    	this->calcAdjointsource(adjsrc_fw, Fwuxsnap->getData(0), Fwuzsnap->getData(0), 0, adjsrc_bw, Bwuxsnap->getData(0), Bwuzsnap->getData(0), 0, model);
    	this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model);

    	// Time stepping displacement
    	waves_adj_fw->forwardstepDisplacement(model, der);
    	waves_adj_bw->forwardstepDisplacement(model, der);

        //Read forward snapshot
        Fwuxsnap->setSnapit(Fwuxsnap->getSnapnt() - 1 - imageit);
        Fwuxsnap->readSnap(it);
        Fwuzsnap->setSnapit(Fwuzsnap->getSnapnt() - 1 - imageit);
        Fwuzsnap->readSnap(it);

        Bwuxsnap->setSnapit(Bwuxsnap->getSnapnt() - 1 - imageit);
        Bwuxsnap->readSnap(it);
        Bwuzsnap->setSnapit(Bwuzsnap->getSnapnt() - 1 - imageit);
        Bwuzsnap->readSnap(it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Fwuxsnap->getEnddiff()) % Fwuxsnap->getSnapinc()) == 0){
            crossCorr(Fwuxsnap->getData(0), Fwuzsnap->getData(0), 0, waves_adj_fw, model, adjsrc_fw);
            crossCorr(Bwuxsnap->getData(0), Bwuzsnap->getData(0), 0, waves_adj_bw, model, adjsrc_bw);
            imageit++;
        }
        
        // Roll pointers
    	waves_adj_fw->roll();
    	waves_adj_bw->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Fwuxsnap->removeSnap();
	Fwuzsnap->removeSnap();

	Bwuxsnap->removeSnap();
	Bwuzsnap->removeSnap();

    // Free arrays
    free(adjsrc_fw);
    free(adjsrc_bw);

    result=RTM_OK;
    return result;
}

template<typename T>
int PPmvaElastic2D<T>::run_optimal(){
     int result = MVA_ERR;
     int nt;
     T dt;
	 T ot;
     T *wsx;
     T *wrx;
     T *wsz;
     T *wrz;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("PPmvaElastic2D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw1 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw2 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_bw1 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_bw2 (new WavesElastic2D_DS<T>(model, nt, dt, ot));

     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_fw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
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

    // Allocate memory for adjoint sources 
    T* adjsrc_fw, *adjsrc_bw;
    adjsrc_fw = (T *) calloc(waves_adj_fw->getNx()*waves_adj_fw->getNz(), sizeof(T));
    adjsrc_bw = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Elastic PP RTMVA gradient with optimal checkpointing.");
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
                waves_fw1->forwardstepStress(model, der);
                waves_bw1->forwardstepStress(model, der);

                // Inserting pressure source 
                waves_fw1->insertPressuresource(model, source, SMAP, it);

                // Inserting data
                waves_bw1->insertPressuresource(model, dataUx, GMAP, nt-1-it);
                waves_bw1->insertPressuresource(model, dataUz, GMAP, nt-1-it);

                // Time stepping displacement
                waves_fw1->forwardstepDisplacement(model, der);
                waves_bw1->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw1->insertForcesource(model, source, SMAP, it);

                // Inserting data
                waves_bw1->insertForcesource(model, dataUx, GMAP, nt-1-it);
                waves_bw1->insertForcesource(model, dataUz, GMAP, nt-1-it);

                // Roll the pointers
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
            waves_fw1->forwardstepStress(model, der);
            waves_bw1->forwardstepStress(model, der);

            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw1->insertPressuresource(model, source, SMAP, capo);
            waves_fw2->insertPressuresource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw1->insertPressuresource(model, dataUx, GMAP, nt-1-capo);
            waves_bw1->insertPressuresource(model, dataUz, GMAP, nt-1-capo);
            waves_bw2->insertPressuresource(model, dataUx, GMAP, capo);
            waves_bw2->insertPressuresource(model, dataUz, GMAP, capo);


            // Inserting adjoint sources
            wsx = waves_fw2->getUx1();
            wsz = waves_fw2->getUz1();
            wrx = waves_bw2->getUx1();
            wrz = waves_bw2->getUz1();
            this->calcAdjointsource(adjsrc_fw, wsx, wsz, waves_fw2->getLpml(), adjsrc_bw, wrx, wrz, waves_bw2->getLpml(), model);
            this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model);

            // Time stepping displacement
            waves_fw1->forwardstepDisplacement(model, der);
            waves_bw1->forwardstepDisplacement(model, der);
            waves_fw2->forwardstepDisplacement(model, der);
            waves_bw2->forwardstepDisplacement(model, der);
            waves_adj_fw->forwardstepDisplacement(model, der);
            waves_adj_bw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw1->insertForcesource(model, source, SMAP, capo);
            waves_fw2->insertForcesource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw1->insertForcesource(model, dataUx, GMAP, nt-1-capo);
            waves_bw1->insertForcesource(model, dataUz, GMAP, nt-1-capo);
            waves_bw2->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw2->insertForcesource(model, dataUz, GMAP, capo);

            /* Do Crosscorrelation */
            wsx = waves_fw1->getUx1();
            wsz = waves_fw1->getUz1();
            crossCorr(wsx, wsz, waves_fw1->getLpml(), waves_adj_fw, model, adjsrc_fw);
            wrx = waves_bw1->getUx1();
            wrz = waves_bw1->getUz1();
            crossCorr(wrx, wrz, waves_bw1->getLpml(), waves_adj_bw, model, adjsrc_bw);

            // Roll the pointers
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
            waves_fw2->forwardstepStress(model, der);
            waves_bw2->forwardstepStress(model, der);

            waves_adj_fw->forwardstepStress(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw2->insertPressuresource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw2->insertPressuresource(model, dataUx, GMAP, capo);
            waves_bw2->insertPressuresource(model, dataUz, GMAP, capo);

            // Inserting adjoint sources
            wsx = waves_fw2->getUx1();
            wsz = waves_fw2->getUz1();
            wrx = waves_bw2->getUx1();
            wrz = waves_bw2->getUz1();
            this->calcAdjointsource(adjsrc_fw, wsx, wsz, waves_fw2->getLpml(), adjsrc_bw, wrx, wrz, waves_bw2->getLpml(), model);
            this->insertAdjointsource(waves_adj_fw, adjsrc_fw, waves_adj_bw, adjsrc_bw, model);

            // Time stepping displacement
            waves_fw2->forwardstepDisplacement(model, der);
            waves_bw2->forwardstepDisplacement(model, der);
            waves_adj_fw->forwardstepDisplacement(model, der);
            waves_adj_bw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw2->insertForcesource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw2->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw2->insertForcesource(model, dataUz, GMAP, capo);

            /* Do Crosscorrelation */
            wsx = waves_fw1->getUx1();
            wsz = waves_fw1->getUz1();
            crossCorr(wsx, wsz, waves_fw1->getLpml(), waves_adj_fw, model, adjsrc_fw);
            wrx = waves_bw1->getUx1();
            wrz = waves_bw1->getUz1();
            crossCorr(wrx, wrz, waves_bw1->getLpml(), waves_adj_bw, model, adjsrc_bw);

            // Roll the pointers 
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

    // Free arrays
    free(adjsrc_fw);
    free(adjsrc_bw);

    return result;
}

template<typename T>
PPmvaElastic2D<T>::~PPmvaElastic2D() {
    // Nothing here
}

// =============== ELASTIC 2D PSmva CLASS =============== //

template<typename T>
PSmvaElastic2D<T>::PSmvaElastic2D(){
    sourceset = false;
    dataUxset = false;
    dataUzset = false;
    modelset = false;
    simageset = false;
    vsgradset = false;
}

template<typename T>
PSmvaElastic2D<T>::PSmvaElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Image2D<T>> _simage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataUx, std::shared_ptr<Data2D<T>> _dataUz, int order, int snapinc):Mva<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUz = _dataUz;
    model = _model;
    simage = _simage;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUzset = true;
    simageset = false;
    vsgradset = false;
}

template<typename T>
bool PSmvaElastic2D<T>::checkStability(){
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
void PSmvaElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads, std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, T *adjsrcxx, T* adjsrczz, T* adjsrcxz)
{
    if(!vsgradset) rs_error("PSmvaElastic2D<T>::crossCorr: No gradient set for computation.");
	int ix, iz;
    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getUx1();
    T* wrz = waves_bw->getUz1();

	T *vsgraddata = NULL; 
	T msxx=0, mszz=0, msxz=0, mrxx=0, mrzz=0, mrxz=0;
	int nx;
	int nz;
	T dx;
	T dz;

    T* Vs = model->getVp();
    T* Rho = model->getR();
    T vsscale;

    if(!vsgrad->getAllocated()){
        vsgrad->allocateImage();
    }
    vsgraddata = vsgrad->getImagedata();

    // Getting sizes
    nx = waves_bw->getNx();
    nz = waves_bw->getNz();
    dx = waves_bw->getDx(); 
    dz = waves_bw->getDz(); 

	int nxs = nx+2*pads;
    int nxr = nx+2*padr;

    for (ix=1; ix<nx-1; ix++){
        for (iz=1; iz<nz-1; iz++){
            msxx = (wsx[ks2D(ix+pads, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads)])/dx;
            mszz = (wsz[ks2D(ix+pads, iz+pads)] - wsz[ks2D(ix+pads, iz+pads-1)])/dz;
            mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
            mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
            msxz = 0.5*(wsx[ks2D(ix+pads, iz+pads+1)] - wsx[ks2D(ix+pads, iz+pads)])/dz;
            msxz += 0.5*(wsx[ks2D(ix+pads-1, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads-1)])/dz;
            msxz += 0.5*(wsz[ks2D(ix+pads+1, iz+pads)] - wsz[ks2D(ix+pads, iz+pads)])/dx;
            msxz += 0.5*(wsz[ks2D(ix+pads, iz+pads-1)] - wsz[ks2D(ix+pads-1, iz+pads-1)])/dx;

            mrxz = 0.5*(wrx[kr2D(ix+padr, iz+padr+1)] - wrx[kr2D(ix+padr, iz+padr)])/dz;
            mrxz += 0.5*(wrx[kr2D(ix+padr-1, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr-1)])/dz;
            mrxz += 0.5*(wrz[kr2D(ix+padr+1, iz+padr)] - wrz[kr2D(ix+padr, iz+padr)])/dx;
            mrxz += 0.5*(wrz[kr2D(ix+padr, iz+padr-1)] - wrz[kr2D(ix+padr-1, iz+padr-1)])/dx;


            vsscale = 2.0*Rho[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            vsgraddata[km2D(ix,iz)] += vsscale*(msxz*mrxz - 2.0*msxx*mrzz - 2.0*mszz*mrxx);
            vsgraddata[km2D(ix,iz)] += vsscale*(msxx*adjsrcxx[km2D(ix,iz)] + mszz*adjsrczz[km2D(ix,iz)] + msxz*adjsrcxz[km2D(ix,iz)]);

        }
    }	
}

template<typename T>
void PSmvaElastic2D<T>::calcAdjointsource(T *adjsrcxx_bw, T *adjsrczz_bw, T *adjsrcxz_bw, T* wsx, T* wsz, int pads, std::shared_ptr<ModelElastic2D<T>> model)
{
	if(!simage->getAllocated()) rs_error("PSmvaElastic2D<T>::calcAdjointsource: simage is not allocated.");
	int ix, iz, ihx, ihz;
	T *imagedata = simage->getImagedata();
	T msxx, mszz, msxz;
    T imagexz;
	T C44_minus;
	T C44_minus_xz;
    T* Vs = model->getVs();
    T* Rho = model->getR();
	int nhx = simage->getNhx();
	int nhz = simage->getNhz();
	int nx = simage->getNx();
	int nxs = nx+2*pads;
	int nz = simage->getNz();
    T dx = simage->getDx(); 
    T dz = simage->getDz(); 
	int hx, hz;


    float C;
    //Reset arrays
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            adjsrcxx_bw[km2D(ix,iz)] = 0.0;
            adjsrczz_bw[km2D(ix,iz)] = 0.0;
            adjsrcxz_bw[km2D(ix,iz)] = 0.0;
        }
    }
    // Calculate integral
    for (ihx=0; ihx<nhx; ihx++){
        hx= -(nhx-1)/2 + ihx;
        C = 1.0;
        if(((ihx == 0) || (ihx == nhx-1)) && (nhx > 1)) C*=0.5;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            if(((ihz == 0) || (ihz == nhz-1)) && (nhz > 1)) C*=0.5;
            for (ix=0; ix<nx; ix++){
                if( ((ix-2*hx) >= 1) && ((ix-2*hx) < nx-1) )
                {
                    for (iz=0; iz<nz; iz++){
                        if( ((iz-2*hz) >= 1) && ((iz-2*hz) < nz-1) )
                        {
							C44_minus = Rho[km2D(ix-2*hx, iz-2*hz)]*Vs[km2D(ix-2*hx, iz-2*hz)]*Vs[km2D(ix-2*hx, iz-2*hz)];
							C44_minus_xz = 0.5*Rho[km2D(ix-2*hx, iz-2*hz)]*Vs[km2D(ix-2*hx, iz-2*hz)]*Vs[km2D(ix-2*hx, iz-2*hz)];
							C44_minus_xz += 0.5*Rho[km2D(ix-2*hx+1, iz-2*hz+1)]*Vs[km2D(ix-2*hx+1, iz-2*hz+1)]*Vs[km2D(ix-2*hx+1, iz-2*hz+1)];
							msxx = (wsx[ks2D(ix-2*hx+pads, iz-2*hz+pads)] - wsx[ks2D(ix-2*hx+pads-1, iz-2*hz+pads)])/dx;
							mszz = (wsz[ks2D(ix-2*hx+pads, iz-2*hz+pads)] - wsz[ks2D(ix-2*hx+pads, iz-2*hz+pads-1)])/dz;

							msxz  = (wsx[ks2D(ix-2*hx+pads, iz-2*hz+pads+1)] - wsx[ks2D(ix-2*hx+pads, iz-2*hz+pads)])/dz;
							msxz += (wsz[ks2D(ix-2*hx+pads+1, iz-2*hz+pads)] - wsz[ks2D(ix-2*hx+pads, iz-2*hz+pads)])/dx;
                            imagexz = 0.25*(imagedata[ki2D(ix-hx,iz-hz,ihx,ihz)] + imagedata[ki2D(ix-hx+1,iz-hz+1,ihx,ihz)] + imagedata[ki2D(ix-hx+1,iz-hz,ihx,ihz)]+imagedata[ki2D(ix-hx,iz-hz+1,ihx,ihz)] );
							adjsrcxx_bw[km2D(ix,iz)] -= 2.0*C*imagedata[ki2D(ix-hx,iz-hz,ihx,ihz)]*C44_minus*mszz;
							adjsrczz_bw[km2D(ix,iz)] -= 2.0*C*imagedata[ki2D(ix-hx,iz-hz,ihx,ihz)]*C44_minus*msxx;
							adjsrcxz_bw[km2D(ix,iz)] += C*imagexz*C44_minus_xz*msxz;
                        }

                    }	
                }
            }
        }
    }
}

template<typename T>
void PSmvaElastic2D<T>::insertAdjointsource(std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, T* adjsrcxx_bw, T* adjsrczz_bw, T* adjsrcxz_bw, std::shared_ptr<ModelElastic2D<T>> model)
{
    int ix, iz;
    int nx = simage->getNx();
    int nz = simage->getNz();
    int padw = waves_bw->getLpml();
    int nxw = nx+2*padw;
	T C44;
    T* Vs = model->getVs();
    T* Rho = model->getR();
    T* wrx = waves_bw->getSxx();
    T* wrz = waves_bw->getSzz();
    T* wrxz = waves_bw->getSxz();
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            C44 = Rho[km2D(ix, iz)]*Vs[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            wrx[kw2D(ix+padw,iz+padw)] += C44*adjsrcxx_bw[km2D(ix, iz)];
            wrz[kw2D(ix+padw,iz+padw)] += C44*adjsrczz_bw[km2D(ix, iz)];
            wrxz[kw2D(ix+padw,iz+padw)] += C44*adjsrcxz_bw[km2D(ix, iz)];
        }
    }	
}

template<typename T>
int PSmvaElastic2D<T>::run(){
	int result = RTM_ERR;
	if(!vsgradset) {
		rs_warning("PSmvaElastic2D::run: No gradient set");
		return result;
	}
	int nt;
	T dt;
	T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("PSmvaElastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_adj_bw->getNx_pml(), 1, waves_adj_bw->getNz_pml(), waves_adj_bw->getDx(), 1.0, waves_adj_bw->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Fwuxsnap;
     Fwuxsnap = std::make_shared<Snapshot2D<T>>(waves_adj_bw, this->getSnapinc());
     Fwuxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
     Fwuxsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Fwuzsnap;
     Fwuzsnap = std::make_shared<Snapshot2D<T>>(waves_adj_bw, this->getSnapinc());
     Fwuzsnap->openSnap(this->getSnapfile() + "-uz", 'r');
     Fwuzsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Bwuxsnap;
     Bwuxsnap = std::make_shared<Snapshot2D<T>>(waves_adj_bw, this->getSnapinc());
     Bwuxsnap->openSnap(this->getSnapfile() + "-ux-bw", 'r');
     Bwuxsnap->allocSnap(0);

     std::shared_ptr<Snapshot2D<T>> Bwuzsnap;
     Bwuzsnap = std::make_shared<Snapshot2D<T>>(waves_adj_bw, this->getSnapinc());
     Bwuzsnap->openSnap(this->getSnapfile() + "-uz-bw", 'r');
     Bwuzsnap->allocSnap(0);

    // Create gradient array
    if(this->vsgradset) vsgrad->allocateImage();

    // Allocate memory for adjoint sources 
    T *adjsrcxx;
    T *adjsrczz;
    T *adjsrcxz;
    adjsrcxx = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));
    adjsrczz = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));
    adjsrcxz = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Elastic PS RTMVA gradient with full checkpointing.");
     this->writeLog("Doing time loop.");
    // Loop over reverse time
    int imageit=0;
    for(int it=0; it < nt; it++)
    {
    	// Time stepping stress
    	waves_adj_bw->forwardstepStress(model, der);
    
        //Read snapshots (interpolation needed here if snapinc not 1!)
        Fwuxsnap->setSnapit(imageit);
        Fwuxsnap->readSnap(it);
        Fwuzsnap->setSnapit(imageit);
        Fwuzsnap->readSnap(it);

        // Inserting adjoint sources
    	this->calcAdjointsource(adjsrcxx, adjsrczz, adjsrcxz, Fwuxsnap->getData(0), Fwuzsnap->getData(0), 0, model);
    	this->insertAdjointsource(waves_adj_bw, adjsrcxx, adjsrczz, adjsrcxz, model);

    	// Time stepping displacement
    	waves_adj_bw->forwardstepDisplacement(model, der);

        //Read snapshots
        Bwuxsnap->setSnapit(Bwuxsnap->getSnapnt() - 1 - imageit);
        Bwuxsnap->readSnap(it);
        Bwuzsnap->setSnapit(Bwuzsnap->getSnapnt() - 1 - imageit);
        Bwuzsnap->readSnap(it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Bwuxsnap->getEnddiff()) % Fwuxsnap->getSnapinc()) == 0){
            crossCorr(Bwuxsnap->getData(0), Bwuzsnap->getData(0), 0, waves_adj_bw, model, adjsrcxx, adjsrczz, adjsrcxz);
            imageit++;
        }
        
        // Roll pointers
    	waves_adj_bw->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Fwuxsnap->removeSnap();
	Fwuzsnap->removeSnap();

	Bwuxsnap->removeSnap();
	Bwuzsnap->removeSnap();

    // Free arrays
    free(adjsrcxx);
    free(adjsrczz);
    free(adjsrcxz);

    result=RTM_OK;
    return result;
}

template<typename T>
int PSmvaElastic2D<T>::run_optimal(){
     int result = MVA_ERR;
     int nt;
     T dt;
	 T ot;
     T *wsx;
     T *wrx;
     T *wsz;
     T *wrz;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("PSmvaElastic2D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw1 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw2 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_bw1 (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_adj_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
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
    if(this->vsgradset) vsgrad->allocateImage();

    // Allocate memory for adjoint sources 
    T *adjsrcxx;
    T *adjsrczz;
    T *adjsrcxz;
    adjsrcxx = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));
    adjsrczz = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));
    adjsrcxz = (T *) calloc(waves_adj_bw->getNx()*waves_adj_bw->getNz(), sizeof(T));

     this->writeLog("Running 2D Elastic PS RTMVA gradient with optimal checkpointing.");
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
                waves_fw1->forwardstepStress(model, der);
                waves_bw1->forwardstepStress(model, der);

                // Inserting pressure source 
                waves_fw1->insertPressuresource(model, source, SMAP, it);

                // Inserting data
                waves_bw1->insertPressuresource(model, dataUx, GMAP, nt-1-it);
                waves_bw1->insertPressuresource(model, dataUz, GMAP, nt-1-it);

                // Time stepping displacement
                waves_fw1->forwardstepDisplacement(model, der);
                waves_bw1->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw1->insertForcesource(model, source, SMAP, it);

                // Inserting data
                waves_bw1->insertForcesource(model, dataUx, GMAP, nt-1-it);
                waves_bw1->insertForcesource(model, dataUz, GMAP, nt-1-it);

                // Roll the pointers 
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
            waves_fw1->forwardstepStress(model, der);
            waves_bw1->forwardstepStress(model, der);

            waves_fw2->forwardstepStress(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw1->insertPressuresource(model, source, SMAP, capo);
            waves_fw2->insertPressuresource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw1->insertPressuresource(model, dataUx, GMAP, nt-1-capo);
            waves_bw1->insertPressuresource(model, dataUz, GMAP, nt-1-capo);

            // Inserting adjoint sources
            wsx = waves_fw2->getUx1();
            wsz = waves_fw2->getUz1();
            this->calcAdjointsource(adjsrcxx, adjsrczz, adjsrcxz, wsx, wsz, waves_fw2->getLpml(), model);
            this->insertAdjointsource(waves_adj_bw, adjsrcxx, adjsrczz, adjsrcxz, model);

            // Time stepping displacement
            waves_fw1->forwardstepDisplacement(model, der);
            waves_bw1->forwardstepDisplacement(model, der);
            waves_fw2->forwardstepDisplacement(model, der);
            waves_adj_bw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw1->insertForcesource(model, source, SMAP, capo);
            waves_fw2->insertForcesource(model, source, SMAP, nt-1-capo);

            // Inserting data
            waves_bw1->insertForcesource(model, dataUx, GMAP, nt-1-capo);
            waves_bw1->insertForcesource(model, dataUz, GMAP, nt-1-capo);

            /* Do Crosscorrelation */
            wrx = waves_bw1->getUx1();
            wrz = waves_bw1->getUz1();
            crossCorr(wrx, wrz, waves_bw1->getLpml(), waves_adj_bw, model, adjsrcxx, adjsrczz, adjsrcxz);

            // Roll the pointers 
            waves_fw1->roll();
            waves_bw1->roll();
            waves_fw2->roll();
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
            waves_fw2->forwardstepStress(model, der);
            waves_adj_bw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw2->insertPressuresource(model, source, SMAP, nt-1-capo);

            // Inserting adjoint sources
            wsx = waves_fw2->getUx1();
            wsz = waves_fw2->getUz1();
            this->calcAdjointsource(adjsrcxx, adjsrczz, adjsrcxz, wsx, wsz, waves_fw2->getLpml(), model);
            this->insertAdjointsource(waves_adj_bw, adjsrcxx, adjsrczz, adjsrcxz, model);

            // Time stepping displacement
            waves_fw2->forwardstepDisplacement(model, der);
            waves_adj_bw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw2->insertForcesource(model, source, SMAP, nt-1-capo);

            /* Do Crosscorrelation */
            wrx = waves_bw1->getUx1();
            wrz = waves_bw1->getUz1();
            crossCorr(wrx, wrz, waves_bw1->getLpml(), waves_adj_bw, model, adjsrcxx, adjsrczz, adjsrcxz);

            // Roll the pointers 
            waves_fw2->roll();
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

    // Free arrays
    free(adjsrcxx);
    free(adjsrczz);
    free(adjsrcxz);

    return result;
}

template<typename T>
PSmvaElastic2D<T>::~PSmvaElastic2D() {
    // Nothing here
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Mva<float>;
template class Mva<double>;

template class MvaAcoustic2D<float>;
template class MvaAcoustic2D<double>;

template class PPmvaElastic2D<float>;
template class PPmvaElastic2D<double>;

template class PSmvaElastic2D<float>;
template class PSmvaElastic2D<double>;

}
