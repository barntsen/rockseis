// Include statements
#include "rtm.h"

namespace rockseis {

// =============== ABSTRACT RTM CLASS =============== //
template<typename T>
Rtm<T>::Rtm() {
	order = 4;
    snapinc=1;
    snapmethod = FULL;
    ncheck = 0;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
}

template<typename T>
Rtm<T>::Rtm(int _order, int _snapinc) {
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
bool Rtm<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return RTM_ERR;
	}else{
		Flog.close();
		return RTM_OK;
	}
}

template<typename T>
void Rtm<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Rtm<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Rtm<T>::writeProgressbar(int x, int n, int r, int w){
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
void Rtm<T>::writeProgress(int x, int n, int r, int w){
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
Rtm<T>::~Rtm() {
    // Nothing here
}

// =============== ACOUSTIC 2D RTM CLASS =============== //

template<typename T>
RtmAcoustic2D<T>::RtmAcoustic2D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    pimageset = false;
}

template<typename T>
RtmAcoustic2D<T>::RtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<Image2D<T>> _pimage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataP, int order, int snapinc):Rtm<T>(order, snapinc){
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
bool RtmAcoustic2D<T>::checkStability(){
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
void RtmAcoustic2D<T>::crossCorr(T *ws, int pads, T* wr, int padr)
{
	if(!pimage->getAllocated()) pimage->allocateImage();
	int ix, iz, ihx, ihz;
	T *imagedata = pimage->getImagedata();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nxs = nx+2*pads;
	int nxr = nx+2*padr;
	int nz = pimage->getNz();
	int hx, hz;
	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihz=0; ihz<nhz; ihz++){
			hz= -(nhz-1)/2 + ihz;
			for (ix=0; ix<nx; ix++){
				if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
				{
					for (iz=0; iz<nz; iz++){
						if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz))
							imagedata[ki2D(ix,iz,ihx,ihz)] += ws[ks2D(ix-hx+pads, iz-hz+pads)]*wr[kr2D(ix+hx+padr, iz+hz+padr)];
					}	
				}
			}
		}
	}
}


template<typename T>
int RtmAcoustic2D<T>::run(){
     int result = RTM_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmAcoustic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap, Bwsnap;
     Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field
     

     this->writeLog("Running 2D Acoustic reverse-time migration with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

    	//Writting out results to snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field
        Psnap->writeSnap(it);

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Psnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesAcoustic2D<T>>(model, nt, dt, ot);

    // Create image
    pimage->allocateImage();

    Psnap->openSnap(this->getSnapfile(), 'r');
    Psnap->allocSnap(0);

    if(this->getRunmva()){
        Bwsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Bwsnap->openSnap(this->getSnapfile() + "-bw", 'w'); // Create a new snapshot file
        Bwsnap->setData(waves->getP1(), 0); //Set Pressure as snap field
    }

     this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataP, GMAP, (nt - 1 - it));

        
        if(this->getRunmva()){
            //Writting out backward modelled wavefield to snapshot file
            Bwsnap->setData(waves->getP1(), 0); //Set Pressure as snap field
            Bwsnap->writeSnap(it);
        }

        //Read forward snapshot
        Psnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
            T *wr = waves->getP1();
            crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
        }

        // Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }

    if(!this->getRunmva()){
        //Remove snapshot file
        Psnap->removeSnap();
    }

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmAcoustic2D<T>::run_edge(){
     int result = RTM_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmAcoustic2D::run_edge: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves_fw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), 1, waves_fw->getNz_pml(), waves_fw->getDx(), 1.0, waves_fw->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     Psnap = std::make_shared<Snapshot2D<T>>(waves_fw, this->getSnapinc());
     Psnap->setData(waves_fw->getP1(), 0); //Set Pressure as snap field
     Psnap->openEdge(this->getSnapfile(), 'w'); // Create a new snapshot file

    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves_fw->forwardstepAcceleration(model, der);
    	waves_fw->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves_fw->insertSource(model, source, SMAP, it);

    	//Writting out results to snapshot file
        Psnap->setData(waves_fw->getP1(), 0); //Set Pressure as snap field
        Psnap->writeEdge(it);

    	// Roll the pointers P1 and P2
    	waves_fw->roll();

        // Output progress to logfile
        this->writeProgress(it, 2*nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Psnap->closeSnap();

    // Reset waves_fw
    waves_fw.reset();
    waves_fw  = std::make_shared<WavesAcoustic2D<T>>(model, nt, dt, ot);
    std::shared_ptr<WavesAcoustic2D<T>> waves_bw (new WavesAcoustic2D<T>(model, nt, dt, ot));

    // Create image
    pimage->allocateImage();

    Psnap.reset();
    Psnap = std::make_shared<Snapshot2D<T>>(waves_fw, this->getSnapinc());
    Psnap->setData(waves_fw->getP2(), 0); //Set Pressure as snap field
    Psnap->openEdge(this->getSnapfile(), 'r');

    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves_fw->forwardstepAcceleration(model, der);
    	waves_fw->forwardstepStress(model, der);

    	waves_bw->forwardstepAcceleration(model, der);
    	waves_bw->forwardstepStress(model, der);

    	// Inserting source 
    	waves_bw->insertSource(model, dataP, GMAP, (nt - 1 - it));

        //Read forward edges
        Psnap->setData(waves_fw->getP2(), 0); //Set Pressure as snap field
        Psnap->readEdge(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
            T *ws = waves_fw->getP1();
            T *wr = waves_bw->getP1();
            crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());
        }

        // Roll the pointers P1 and P2
    	waves_fw->roll();
    	waves_bw->roll();

        // Output progress to logfile
        this->writeProgress(nt-1 + it, 2*nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Psnap->removeSnap();

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmAcoustic2D<T>::run_optimal(){
     int result = RTM_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmAcoustic2D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves_fw (new WavesAcoustic2D<T>(model, nt, dt, ot));
    std::shared_ptr<WavesAcoustic2D<T>> waves_bw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), 1, waves_fw->getNz_pml(), waves_fw->getDx(), 1.0, waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
     pimage->allocateImage();

     this->writeLog("Running 2D Acoustic reverse-time migration with optimal checkpointing.");
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
                waves_fw->forwardstepAcceleration(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Roll the pointers P1 and P2
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
            waves_fw->forwardstepAcceleration(model, der);
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertSource(model, dataP, GMAP, capo);

            /* Do Crosscorrelation */
            T *ws = waves_fw->getP1();
            T *wr = waves_bw->getP1();
            crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

            // Roll the pointers P1 and P2
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
            waves_bw->forwardstepAcceleration(model, der);
            waves_bw->forwardstepStress(model, der);

            // Inserting data
            waves_bw->insertSource(model, dataP, GMAP, capo);

            /* Do Crosscorrelation */
            T *ws = waves_fw->getP1();
            T *wr = waves_bw->getP1();
            crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

            // Roll the pointers P1 and P2
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

    result=RTM_OK;
    return result;
}

template<typename T>
RtmAcoustic2D<T>::~RtmAcoustic2D() {
    // Nothing here
}

// =============== ACOUSTIC 3D RTM CLASS =============== //

template<typename T>
RtmAcoustic3D<T>::RtmAcoustic3D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    pimageset = false;
}

template<typename T>
RtmAcoustic3D<T>::RtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model, std::shared_ptr<Image3D<T>> _pimage, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataP, int order, int snapinc):Rtm<T>(order, snapinc){
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
bool RtmAcoustic3D<T>::checkStability(){
    T *Vp = model->getVp();
    // Find maximum Vp
    T Vpmax;
    Vpmax=Vp[0];
    size_t n=model->getNx()*model->getNy()*model->getNz();
    for(size_t i=1; i<n; i++){
        if(Vp[i] > Vpmax){
            Vpmax = Vp[i];
        }
    }

    T dx = model->getDx();
    T dy = model->getDy();
    T dz = model->getDz();
    T dt = source->getDt();
    T dt_stab;
    dt_stab = 2.0/(3.1415*sqrt((1.0/(dx*dx))+(1/(dy*dy))+(1/(dz*dz)))*Vpmax); 
    if(dt < dt_stab){
        return true;
    }else{
        rs_warning("Modeling time interval exceeds maximum stable number of: ", std::to_string(dt_stab));
        return false;
    }
}

template<typename T>
void RtmAcoustic3D<T>::crossCorr(T *ws, int pads, T* wr, int padr)
{
	if(!pimage->getAllocated()) pimage->allocateImage();
	int ix, iy, iz, ihx, ihy, ihz;
	int nhx = pimage->getNhx();
	int nhy = pimage->getNhy();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nxs = nx + 2*pads;
	int nxr = nx + 2*padr;
	int ny = pimage->getNy();
	int nys = ny + 2*pads;
	int nyr = ny + 2*padr;
	int nz = pimage->getNz();
	int hx, hy, hz;
	T* imagedata = pimage->getImagedata();
	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihy=0; ihy<nhy; ihy++){
			hy= -(nhy-1)/2 + ihy;
			for (ihz=0; ihz<nhz; ihz++){
				hz= -(nhz-1)/2 + ihz;
				for (ix=0; ix<nx; ix++){
					if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
						for (iy=0; iy<ny; iy++){
							if( ((iy-hy) >= 0) && ((iy-hy) < ny) && ((iy+hy) >= 0) && ((iy+hy) < ny))
								for (iz=0; iz<nz; iz++){
									if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz))
										imagedata[ki3D(ix,iy,iz,ihx,ihy,ihz)] += ws[ks3D(ix-hx+pads,iy-hy+pads,iz-hz+pads)]*wr[kr3D(ix+hx+padr,iy+hy+padr,iz+hz+padr)];
								}	
						}
				}
			}
		}
	}
}


template<typename T>
int RtmAcoustic3D<T>::run(){
     int result = RTM_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmAcoustic3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

     this->writeLog("Running 3D Acoustic reverse-time migration with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

    	//Writting out results to snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field
        Psnap->writeSnap(it);

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Psnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesAcoustic3D<T>>(model, nt, dt, ot);

    // Create image
    pimage->allocateImage();

    Psnap->openSnap(this->getSnapfile(), 'r');
    Psnap->allocSnap(0);

     this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataP, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Psnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
            T *wr = waves->getP1();
            crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);

        // Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    

	//Remove snapshot file
	Psnap->removeSnap();

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmAcoustic3D<T>::run_optimal(){
     int result = RTM_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmAcoustic3D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves_fw (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic3D<T>> waves_bw (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), waves_fw->getNy_pml(), waves_fw->getNz_pml(), waves_fw->getDx(), waves_fw->getDy(), waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
     pimage->allocateImage();


     this->writeLog("Running 3D Acoustic reverse-time migration with optimal checkpointing.");
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
                waves_fw->forwardstepAcceleration(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Roll the pointers P1 and P2
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
            waves_fw->forwardstepAcceleration(model, der);
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertSource(model, dataP, GMAP, capo);

            /* Do Crosscorrelation */
            T *ws = waves_fw->getP1();
            T *wr = waves_bw->getP1();
            crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

            // Roll the pointers P1 and P2
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
            waves_bw->forwardstepAcceleration(model, der);
            waves_bw->forwardstepStress(model, der);

            // Inserting data
            waves_bw->insertSource(model, dataP, GMAP, capo);

            /* Do Crosscorrelation */
            T *ws = waves_fw->getP1();
            T *wr = waves_bw->getP1();
            crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

            // Roll the pointers P1 and P2
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

    result=RTM_OK;
    return result;
}

template<typename T>
RtmAcoustic3D<T>::~RtmAcoustic3D() {
    // Nothing here
}

// =============== ELASTIC 2D RTM CLASS =============== //

template<typename T>
RtmElastic2D<T>::RtmElastic2D(){
    sourceset = false;
    dataUxset = false;
    dataUzset = false;
    modelset = false;
    pimageset = false;
    simageset = false;
}

template<typename T>
RtmElastic2D<T>::RtmElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataUx, std::shared_ptr<Data2D<T>> _dataUz, int order, int snapinc):Rtm<T>(order, snapinc){
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
bool RtmElastic2D<T>::checkStability(){
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
void RtmElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads, T* wrx, T* wrz, int padr, T* Vp, T* Vs, T* Rho)
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

template<typename T>
int RtmElastic2D<T>::run(){
	int result = RTM_ERR;
	if(!pimageset && !simageset) {
		rs_warning("RtmElastic2D::run: No image set");
		return result;
	}
	int nt;
	T dt;
	T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmElastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Uxsnap, Bwuxsnap;
     Uxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Uxsnap->openSnap(this->getSnapfile() + "-ux", 'w'); // Create a new snapshot file
     Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field

     std::shared_ptr<Snapshot2D<T>> Uzsnap, Bwuzsnap;
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

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting source 
    	waves->insertForcesource(model, source, SMAP, it);

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

    if(this->getRunmva()){
        Bwuxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Bwuxsnap->openSnap(this->getSnapfile() + "-ux-bw", 'w'); // Create a new snapshot file
        Bwuxsnap->setData(waves->getUx1(), 0); //Set Displacement field as snap field

        Bwuzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Bwuzsnap->openSnap(this->getSnapfile() + "-uz-bw", 'w'); // Create a new snapshot file
        Bwuzsnap->setData(waves->getUz1(), 0); //Set Displacement field as snap field
    }

    // Get models for scaling
    T *Vp, *Vs, *Rho;
    Vp = model->getVp();
    Vs = model->getVs();
    Rho = model->getR();

    this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping stress 
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertPressuresource(model, dataUx, GMAP, (nt - 1 - it));
    	waves->insertPressuresource(model, dataUz, GMAP, (nt - 1 - it));

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);

    	// Inserting source 
    	waves->insertForcesource(model, dataUx, GMAP, (nt - 1 - it));
    	waves->insertForcesource(model, dataUz, GMAP, (nt - 1 - it));

        if(this->getRunmva()){
            //Writting out results to snapshot files
            Bwuxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field
            Bwuxsnap->writeSnap(it);

            Bwuzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field
            Bwuzsnap->writeSnap(it);
        }

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
    
    if(!this->getRunmva()){
        //Remove snapshot file
        Uxsnap->removeSnap();
        Uzsnap->removeSnap();
    }

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmElastic2D<T>::run_optimal(){
     int result = RTM_ERR;
     if(!pimageset && !simageset) {
         rs_warning("RtmElastic2D::run: No image set");
         return result;
     }
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmElastic2D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

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
    if(this->pimageset) pimage->allocateImage();
    if(this->simageset) simage->allocateImage();

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
                // Time stepping stress
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertPressuresource(model, source, SMAP, it);

                // Time stepping displacement
                waves_fw->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw->insertForcesource(model, source, SMAP, it);

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
            // Time stepping stress
            waves_fw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw->insertPressuresource(model, source, SMAP, capo);

            // Time stepping displacement
            waves_fw->forwardstepDisplacement(model, der);

            // Inserting force source 
            waves_fw->insertForcesource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertPressuresource(model, dataUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataUz, GMAP, capo);
            waves_bw->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataUz, GMAP, capo);

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
            // Time stepping stress
            waves_bw->forwardstepStress(model, der);

            // Inserting data
            waves_bw->insertPressuresource(model, dataUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataUz, GMAP, capo);

            // Time stepping displacement 
            waves_bw->forwardstepDisplacement(model, der);

            // Inserting data
            waves_bw->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataUz, GMAP, capo);

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

    result=RTM_OK;
    return result;
}

template<typename T>
RtmElastic2D<T>::~RtmElastic2D() {
    // Nothing here
}

// =============== ELASTIC 3D RTM CLASS =============== //

template<typename T>
RtmElastic3D<T>::RtmElastic3D(){
    sourceset = false;
    dataUxset = false;
    dataUyset = false;
    dataUzset = false;
    modelset = false;
    pimageset = false;
    simageset = false;
}

template<typename T>
RtmElastic3D<T>::RtmElastic3D(std::shared_ptr<ModelElastic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataUx, std::shared_ptr<Data3D<T>> _dataUy, std::shared_ptr<Data3D<T>> _dataUz, int order, int snapinc):Rtm<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUy = _dataUy;
    dataUz = _dataUz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUyset = true;
    dataUzset = true;
    pimageset = false;
    simageset = false;
}

template<typename T>
bool RtmElastic3D<T>::checkStability(){
    T *Vp = model->getVp();
    // Find maximum Vp
    T Vpmax;
    Vpmax=Vp[0];
    size_t n=model->getNx()*model->getNy()*model->getNz();
    for(size_t i=1; i<n; i++){
        if(Vp[i] > Vpmax){
            Vpmax = Vp[i];
        }
    }

    T dx = model->getDx();
    T dy = model->getDy();
    T dz = model->getDz();
    T dt = source->getDt();
    T dt_stab;
    dt_stab = 2.0/(3.1415*sqrt((1.0/(dx*dx))+(1/(dy*dy))+(1/(dz*dz)))*Vpmax); 
    if(dt < dt_stab){
        return true;
    }else{
        rs_warning("Modeling time interval exceeds maximum stable number of: ", std::to_string(dt_stab));
        return false;
    }
}

template<typename T>
void RtmElastic3D<T>::crossCorr(T *wsx, T*wsy, T *wsz, int pads, T* wrx, T* wry, T* wrz, int padr, T* Vp, T* Vs, T* Rho)
{
	int ix, iy, iz, ihx, ihy, ihz;
	T *pimagedata = NULL; 
	T *simagedata = NULL;
	T msxx, msyy, mszz, msyz, msxz, msxy, mrxx, mryy, mrzz, mryz, mrxz, mrxy;
	T C33_minus, C33_plus;
	T C44_minus, C44_plus;
	int nhx; 
	int nhy; 
	int nhz;
	int nx;
	int ny;
	int nz;
	T dx;
	T dy;
	T dz;
	int hx, hy,hz;

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
		nhy = pimage->getNhy();
		nhz = pimage->getNhz();
		nx = pimage->getNx();
		ny = pimage->getNy();
		nz = pimage->getNz();
		dx = pimage->getDx(); 
		dy = pimage->getDy(); 
		dz = pimage->getDz(); 
	}else{
		nhx = simage->getNhx();
		nhy = simage->getNhy();
		nhz = simage->getNhz();
		nx = simage->getNx();
		ny = simage->getNy();
		nz = simage->getNz();
		dx = simage->getDx(); 
		dy = simage->getDy(); 
		dz = simage->getDz(); 
	}

	int nxs = nx + 2*pads;
	int nxr = nx + 2*padr;
	int nys = ny + 2*pads;
	int nyr = ny + 2*padr;

	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihy=0; ihy<nhy; ihy++){
			hy= -(nhy-1)/2 + ihy;
			for (ihz=0; ihz<nhz; ihz++){
				hz= -(nhz-1)/2 + ihz;
				for (ix=0; ix<nx; ix++){
					if( ((ix-hx) >= 1) && ((ix-hx) < nx-1) && ((ix+hx) >= 1) && ((ix+hx) < nx-1))
						for (iy=0; iy<ny; iy++){
							if( ((iy-hy) >= 1) && ((iy-hy) < ny-1) && ((iy+hy) >= 1) && ((iy+hy) < ny-1))
								for (iz=0; iz<nz; iz++){
									if( ((iz-hz) >= 1) && ((iz-hz) < nz-1) && ((iz+hz) >= 1) && ((iz+hz) < nz-1))

									{
										C33_minus = Rho[km3D(ix-hx, iy-hy, iz-hz)]*Vp[km3D(ix-hx, iy-hy, iz-hz)]*Vp[km3D(ix-hx, iy-hy, iz-hz)];
										C33_plus = Rho[km3D(ix+hx, iy+hy, iz+hz)]*Vp[km3D(ix+hx, iy+hy, iz+hz)]*Vp[km3D(ix+hx, iy+hy, iz+hz)];
										C44_minus = Rho[km3D(ix-hx, iy-hy, iz-hz)]*Vs[km3D(ix-hx, iy-hy, iz-hz)]*Vs[km3D(ix-hx, iy-hy, iz-hz)];
										C44_plus = Rho[km3D(ix+hx, iy+hy, iz+hz)]*Vs[km3D(ix+hx, iy+hy, iz+hz)]*Vs[km3D(ix+hx, iy+hy, iz+hz)];

										msxx = (wsx[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)] - wsx[ks3D(ix-hx+pads-1, iy-hy+pads, iz-hz+pads)])/dx;
										msyy = (wsy[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)] - wsy[ks3D(ix-hx+pads, iy-hy+pads-1, iz-hz+pads)])/dy;
										mszz = (wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)] - wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads-1)])/dz;
										mrxx = (wrx[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)] - wrx[kr3D(ix+hx+padr-1, iy+hy+padr, iz+hz+padr)])/dx;
										mryy = (wry[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)] - wry[kr3D(ix+hx+padr, iy+hy+padr-1, iz+hz+padr)])/dy;
										mrzz = (wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)] - wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr-1)])/dz;

										if(pimageset){
											pimagedata[ki3D(ix,iy,iz,ihx,ihy,ihz)] += C33_minus*C33_plus*(msxx + msyy + mszz) * (mrxx + mryy + mrzz);
										}

										if(simageset){
											msyz = 0.5*(wsz[ks3D(ix-hx+pads, iy-hy+pads+1, iz-hz+pads)] - wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dy;
											msyz += 0.5*(wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads-1)] - wsz[ks3D(ix-hx+pads, iy-hy+pads-1, iz-hz+pads-1)])/dy;
											msyz += 0.5*(wsy[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads+1)] - wsy[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dz;
											msyz += 0.5*(wsy[ks3D(ix-hx+pads, iy-hy+pads-1, iz-hz+pads)] - wsy[ks3D(ix-hx+pads, iy-hy+pads-1, iz-hz+pads-1)])/dz;

											mryz = 0.5*(wrz[kr3D(ix+hx+padr, iy+hy+padr+1, iz+hz+padr)] - wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dy;
											mryz += 0.5*(wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr-1)] - wrz[kr3D(ix+hx+padr, iy+hy+padr-1, iz+hz+padr-1)])/dy;
											mryz += 0.5*(wry[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr+1)] - wry[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dz;
											mryz += 0.5*(wry[kr3D(ix+hx+padr, iy+hy+padr-1, iz+hz+padr)] - wry[kr3D(ix+hx+padr, iy+hy+padr-1, iz+hz+padr-1)])/dz;

											msxz = 0.5*(wsx[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads+1)] - wsx[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dz;
											msxz += 0.5*(wsx[ks3D(ix-hx+pads-1, iy-hy+pads, iz-hz+pads)] - wsx[ks3D(ix-hx+pads-1, iy-hy+pads, iz-hz+pads-1)])/dz;
											msxz += 0.5*(wsz[ks3D(ix-hx+pads+1, iy-hy+pads, iz-hz+pads)] - wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dx;
											msxz += 0.5*(wsz[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads-1)] - wsz[ks3D(ix-hx+pads-1, iy-hy+pads, iz-hz+pads-1)])/dx;

											mrxz = 0.5*(wrx[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr+1)] - wrx[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dz;
											mrxz += 0.5*(wrx[kr3D(ix+hx+padr-1, iy+hy+padr, iz+hz+padr)] - wrx[kr3D(ix+hx+padr-1, iy+hy+padr, iz+hz+padr-1)])/dz;
											mrxz += 0.5*(wrz[kr3D(ix+hx+padr+1, iy+hy+padr, iz+hz+padr)] - wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dx;
											mrxz += 0.5*(wrz[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr-1)] - wrz[kr3D(ix+hx+padr-1, iy+hy+padr, iz+hz+padr-1)])/dx;

											msxy = 0.5*(wsx[ks3D(ix-hx+pads, iy-hy+pads+1, iz-hz+pads)] - wsx[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dy;
											msxy += 0.5*(wsx[ks3D(ix-hx+pads-1, iy-hy+pads, iz-hz+pads)] - wsx[ks3D(ix-hx+pads-1, iy-hy+pads-1, iz-hz+pads)])/dy;
											msxy += 0.5*(wsy[ks3D(ix-hx+pads+1, iy-hy+pads, iz-hz+pads)] - wsy[ks3D(ix-hx+pads, iy-hy+pads, iz-hz+pads)])/dx;
											msxy += 0.5*(wsy[ks3D(ix-hx+pads, iy-hy+pads-1, iz-hz+pads)] - wsy[ks3D(ix-hx+pads-1, iy-hy+pads-1, iz-hz+pads)])/dx;

											mrxy = 0.5*(wrx[kr3D(ix+hx+padr, iy+hy+padr+1, iz+hz+padr)] - wrx[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dy;
											mrxy += 0.5*(wrx[kr3D(ix+hx+padr-1, iy+hy+padr, iz+hz+padr)] - wrx[kr3D(ix+hx+padr-1, iy+hy+padr-1, iz+hz+padr)])/dy;
											mrxy += 0.5*(wry[kr3D(ix+hx+padr+1, iy+hy+padr, iz+hz+padr)] - wry[kr3D(ix+hx+padr, iy+hy+padr, iz+hz+padr)])/dx;
											mrxy += 0.5*(wry[kr3D(ix+hx+padr, iy+hy+padr-1, iz+hz+padr)] - wry[kr3D(ix+hx+padr-1, iy+hy+padr-1, iz+hz+padr)])/dx;

											simagedata[ki3D(ix,iy,iz,ihx,ihy,ihz)] += C44_minus*C44_plus*(-2.0*(msyy*mrzz + mszz*mryy) -2.0*(msxx*mrzz + mszz*mrxx) -2.0*(msyy*mrxx + msxx*mryy) + msyz*mryz + msxz*mrxz + msxy*mrxy);
										}
									}
								}	
						}
				}
			}
		}
	}
}

template<typename T>
int RtmElastic3D<T>::run(){
     int result = RTM_ERR;
     if(!pimageset && !simageset) {
         rs_warning("RtmElastic3D::run: No image set");
         return result;
     }
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmElastic3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D_DS<T>> waves (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Uxsnap;
     Uxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uxsnap->openSnap(this->getSnapfile() + "-ux", 'w'); // Create a new snapshot file
     Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field

     std::shared_ptr<Snapshot3D<T>> Uysnap;
     Uysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uysnap->openSnap(this->getSnapfile() + "-uy", 'w'); // Create a new snapshot file
     Uysnap->setData(waves->getUy1(), 0); //Set Uy as snap field

     std::shared_ptr<Snapshot3D<T>> Uzsnap;
     Uzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uzsnap->openSnap(this->getSnapfile() + "-uz", 'w'); // Create a new snapshot file
     Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field

     this->writeLog("Running 3D Elastic reverse-time migration with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field
        Uxsnap->writeSnap(it);

        Uysnap->setData(waves->getUy1(), 0); //Set Uy as snap field
        Uysnap->writeSnap(it);

        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field
        Uzsnap->writeSnap(it);

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting source 
    	waves->insertForcesource(model, source, SMAP, it);

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Uxsnap->closeSnap();
    Uysnap->closeSnap();
    Uzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic3D_DS<T>>(model, nt, dt, ot);

    // Create image
    if(this->pimageset) pimage->allocateImage();
    if(this->simageset) simage->allocateImage();

    Uxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
    Uxsnap->allocSnap(0);

    Uysnap->openSnap(this->getSnapfile() + "-uy", 'r');
    Uysnap->allocSnap(0);

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
    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);

    	// Inserting source 
    	waves->insertForcesource(model, dataUx, GMAP, (nt - 1 - it));
    	waves->insertForcesource(model, dataUy, GMAP, (nt - 1 - it));
    	waves->insertForcesource(model, dataUz, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Uxsnap->readSnap(nt - 1 - it);
        Uysnap->readSnap(nt - 1 - it);
        Uzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Uxsnap->getEnddiff()) % Uxsnap->getSnapinc()) == 0){
            T *Uxr = waves->getUx1();
            T *Uyr = waves->getUy1();
            T *Uzr = waves->getUz1();
            crossCorr(Uxsnap->getData(0), Uysnap->getData(0), Uzsnap->getData(0), 0, Uxr, Uyr, Uzr, waves->getLpml(), Vp, Vs, Rho);
        }

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Uxsnap->removeSnap();
	Uysnap->removeSnap();
	Uzsnap->removeSnap();

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmElastic3D<T>::run_optimal(){
     int result = RTM_ERR;
     if(!pimageset && !simageset) {
         rs_warning("RtmElastic3D::run: No image set");
         return result;
     }
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("RtmElastic3D::run_optimal: Wavelet sampling interval (dt) does not match the stability criteria.");

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D_DS<T>> waves_fw (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic3D_DS<T>> waves_bw (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), waves_fw->getNy_pml(), waves_fw->getNz_pml(), waves_fw->getDx(), waves_fw->getDy(), waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');


    // Create image
    if(this->pimageset) pimage->allocateImage();
    if(this->simageset) simage->allocateImage();

     // Get models for scaling
     T *Vp, *Vs, *Rho;
     Vp = model->getVp();
     Vs = model->getVs();
     Rho = model->getR();

     this->writeLog("Running 3D Elastic reverse-time migration with optimal checkpointing.");
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
                // Time stepping stress
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertPressuresource(model, source, SMAP, it);

                // Time stepping displacement
                waves_fw->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw->insertForcesource(model, source, SMAP, it);


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
            // Time stepping stress
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertPressuresource(model, source, SMAP, capo);

            // Time stepping displacement
            waves_fw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw->insertForcesource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataUz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getUx1();
            T *wsy = waves_fw->getUy1();
            T *wsz = waves_fw->getUz1();
            T *wrx = waves_bw->getUx1();
            T *wry = waves_bw->getUy1();
            T *wrz = waves_bw->getUz1();

            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), wrx, wry, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

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
            // Time stepping stress
            waves_bw->forwardstepStress(model, der);

            // Time stepping displacement
            waves_bw->forwardstepDisplacement(model, der);

            // Inserting data
            waves_bw->insertForcesource(model, dataUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataUy, GMAP, capo);
            waves_bw->insertForcesource(model, dataUz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getUx1();
            T *wsy = waves_fw->getUy1();
            T *wsz = waves_fw->getUz1();
            T *wrx = waves_bw->getUx1();
            T *wry = waves_bw->getUy1();
            T *wrz = waves_bw->getUz1();
            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), wrx, wry, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

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

    result=RTM_OK;
    return result;
}

template<typename T>
RtmElastic3D<T>::~RtmElastic3D() {
    // Nothing here
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Rtm<float>;
template class Rtm<double>;
template class RtmAcoustic2D<float>;
template class RtmAcoustic2D<double>;
template class RtmAcoustic3D<float>;
template class RtmAcoustic3D<double>;

template class RtmElastic2D<float>;
template class RtmElastic2D<double>;
template class RtmElastic3D<float>;
template class RtmElastic3D<double>;

}
