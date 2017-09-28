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
void Rtm<T>::writeLog(char *text){
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
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

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
        this->writeProgress(it, 2*nt-1, 20, 48);
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

    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataP, GMAP, (nt - 1 - it));

        //Get forward snapshot
        Psnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
            T *wr = waves->getP1();
            crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
        }

        // Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(nt-1 + it, 2*nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Psnap->removeSnap();

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmAcoustic2D<T>::run_edge(){
     int result = RTM_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

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

        //Get forward edges
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
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

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

            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
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

        // Output progress to logfile
        //this->writeProgress(capo, nt, 20, 48);
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
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

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
        this->writeProgress(it, 2*nt-1, 20, 48);
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

    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataP, GMAP, (nt - 1 - it));

        //Get forward snapshot
        Psnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
            T *wr = waves->getP1();
            crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
        }

        // Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(nt-1 + it, 2*nt-1, 20, 48);
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
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

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

            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
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

        // Output progress to logfile
        //this->writeProgress(capo, nt, 20, 48);
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
    dataVxset = false;
    dataVzset = false;
    modelset = false;
    pimageset = false;
    simageset = false;
}

template<typename T>
RtmElastic2D<T>::RtmElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataVx, std::shared_ptr<Data2D<T>> _dataVz, int order, int snapinc):Rtm<T>(order, snapinc){
    source = _source;
    dataVx = _dataVx;
    dataVz = _dataVz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataVxset = true;
    dataVzset = true;
    pimageset = false;
    simageset = false;
}

template<typename T>
int RtmElastic2D<T>::run(){
     int result = RTM_ERR;
     if(!pimageset && !simageset) {
         rs_warning("RtmElastic2D::run: No image set");
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
     std::shared_ptr<WavesElastic2D<T>> waves (new WavesElastic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Vxsnap;
     Vxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Vxsnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field

     std::shared_ptr<Snapshot2D<T>> Vzsnap;
     Vzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Vzsnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field

    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

    	//Writting out results to snapshot files
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field
        Vxsnap->writeSnap(it);

        Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field
        Vzsnap->writeSnap(it);

        // Output progress to logfile
        this->writeProgress(it, 2*nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Vxsnap->closeSnap();
    Vzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic2D<T>>(model, nt, dt, ot);

    // Create image
    if(this->pimageset) pimage->allocateImage();
    if(this->simageset) simage->allocateImage();

    Vxsnap->openSnap(this->getSnapfile(), 'r');
    Vxsnap->allocSnap(0);

    Vzsnap->openSnap(this->getSnapfile(), 'r');
    Vzsnap->allocSnap(0);

    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);

    	// Inserting source 
    	waves->insertSource(model, dataVx, GMAP, (nt - 1 - it));
    	waves->insertSource(model, dataVz, GMAP, (nt - 1 - it));

        //Get forward snapshot
        Vxsnap->readSnap(nt - 1 - it);
        Vzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Vxsnap->getEnddiff()) % Vxsnap->getSnapinc()) == 0){
            T *Vxr = waves->getVx();
            T *Vzr = waves->getVz();
            //crossCorr_pp(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ppps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
        }

        // Output progress to logfile
        this->writeProgress(nt-1 + it, 2*nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Vxsnap->removeSnap();
	Vzsnap->removeSnap();

    result=RTM_OK;
    return result;
}

template<typename T>
int RtmElastic2D<T>::run_optimal(){
     int result = RTM_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D<T>> waves_fw (new WavesElastic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D<T>> waves_bw (new WavesElastic2D<T>(model, nt, dt, ot));
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
                waves_fw->forwardstepVelocity(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepVelocity(model, der);
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Inserting data
            waves_bw->insertSource(model, dataVx, GMAP, capo);
            waves_bw->insertSource(model, dataVz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getVx();
            T *wsz = waves_fw->getVz();
            T *wrx = waves_bw->getVx();
            T *wrz = waves_bw->getVz();

            //crossCorr_pp(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ppps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());

            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
        }
        if (whatodo == youturn)
        {
            // Time stepping
            waves_bw->forwardstepVelocity(model, der);
            waves_bw->forwardstepStress(model, der);

            // Inserting data
            waves_bw->insertSource(model, dataVx, GMAP, capo);
            waves_bw->insertSource(model, dataVz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getVx();
            T *wsz = waves_fw->getVz();
            T *wrx = waves_bw->getVx();
            T *wrz = waves_bw->getVz();
            //crossCorr_pp(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
            //crossCorr_ppps(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml());
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

        // Output progress to logfile
        //this->writeProgress(capo, nt, 20, 48);
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

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Rtm<float>;
template class Rtm<double>;
template class RtmAcoustic2D<float>;
template class RtmAcoustic2D<double>;
template class RtmAcoustic3D<float>;
template class RtmAcoustic3D<double>;

}
