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
RtmAcoustic2D<T>::RtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<ImageAcoustic2D<T>> _pimage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataP, int order, int snapinc):Rtm<T>(order, snapinc){
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
            pimage->crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
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
            pimage->crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());
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
            pimage->crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

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
            pimage->crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

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
RtmAcoustic3D<T>::RtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model, std::shared_ptr<ImageAcoustic3D<T>> _pimage, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataP, int order, int snapinc):Rtm<T>(order, snapinc){
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
            pimage->crossCorr(Psnap->getData(0), 0, wr, waves->getLpml());
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
            pimage->crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

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
            pimage->crossCorr(ws, waves_fw->getLpml(), wr, waves_bw->getLpml());

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

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Rtm<float>;
template class Rtm<double>;
template class RtmAcoustic2D<float>;
template class RtmAcoustic2D<double>;
template class RtmAcoustic3D<float>;
template class RtmAcoustic3D<double>;

}
