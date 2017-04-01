// Include statements
#include "modelling.h"

namespace rockseis {

// =============== ABSTRACT MODELLING CLASS =============== //
template<typename T>
Modelling<T>::Modelling() {
	order = 4;
    snapinc=1;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
}

template<typename T>
Modelling<T>::Modelling(int _order, int _snapinc) {
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

	prog.previous = 0;
	prog.current = 0;
    prog.persec = 1;
}

template<typename T>
bool Modelling<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return MOD_ERR;
	}else{
		Flog.close();
		return MOD_OK;
	}
}

template<typename T>
void Modelling<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Modelling<T>::writeLog(char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Modelling<T>::writeProgressbar(int x, int n, int r, int w){
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
void Modelling<T>::writeProgress(int x, int n, int r, int w){
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
Modelling<T>::~Modelling() {
    // Nothing here
}

// =============== ACOUSTIC 2D MODELLING CLASS =============== //

template<typename T>
ModellingAcoustic2D<T>::ModellingAcoustic2D(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recAxset = false;
    recAzset = false;
    snapPset = false;
    snapAxset = false;
    snapAzset = false;
}

template<typename T>
ModellingAcoustic2D<T>::ModellingAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recAxset = false;
    recAzset = false;
    snapPset = false;
    snapAxset = false;
    snapAzset = false;
}

template<typename T>
int ModellingAcoustic2D<T>::run(){
     int result = MOD_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog("log.txt");

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot, this->getSnapinc()));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

    // Create snapshot for pressure
    if(this->snapPset){ 
        waves->createSnap(this->snapP, waves->getPsnap());
    }
    if(this->snapAxset){ 
        waves->createSnap(this->snapAx, waves->getAxsnap());
    }
    if(this->snapAzset){ 
        waves->createSnap(this->snapAz, waves->getAzsnap());
    }

    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source (Pressure)
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Pressure)
        if(this->recPset){
            waves->recordData(this->recP, GMAP, it);
        }

        if(this->recAxset){
            waves->recordData(this->recAx, GMAP, it);
        }

        if(this->recAzset){
            waves->recordData(this->recAz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            waves->writeSnap(it, waves->getPsnap());
        }

        if(this->snapAxset){ 
            waves->writeSnap(it, waves->getAxsnap());
        }

        if(this->snapAzset){ 
            waves->writeSnap(it, waves->getAzsnap());
        }

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    result=MOD_OK;
    return result;
}

template<typename T>
ModellingAcoustic2D<T>::~ModellingAcoustic2D() {
    // Nothing here
}

// =============== ACOUSTIC 3D MODELLING CLASS =============== //

template<typename T>
ModellingAcoustic3D<T>::ModellingAcoustic3D(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recAxset = false;
    recAyset = false;
    recAzset = false;
    snapPset = false;
    snapAxset = false;
    snapAyset = false;
    snapAzset = false;
}

template<typename T>
ModellingAcoustic3D<T>::ModellingAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recAxset = false;
    recAyset = false;
    recAzset = false;
    snapPset = false;
    snapAxset = false;
    snapAyset = false;
    snapAzset = false;
}

template<typename T>
int ModellingAcoustic3D<T>::run(){
     int result = MOD_ERR;
     int nt;
     float dt;
     float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves (new WavesAcoustic3D<T>(model, nt, dt, ot, this->getSnapinc()));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

	// Create log file
	this->createLog("log.txt");

    // Create snapshot for pressure
    if(this->snapPset){ 
        waves->createSnap(this->snapP, waves->getPsnap());
    }
    if(this->snapAxset){ 
        waves->createSnap(this->snapAx, waves->getAxsnap());
    }
    if(this->snapAyset){ 
        waves->createSnap(this->snapAy, waves->getAysnap());
    }
    if(this->snapAzset){ 
        waves->createSnap(this->snapAz, waves->getAzsnap());
    }

    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source (Pressure)
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Pressure)
        if(this->recPset){
            waves->recordData(this->recP, GMAP, it);
        }

        if(this->recAxset){
            waves->recordData(this->recAx, GMAP, it);
        }

        if(this->recAyset){
            waves->recordData(this->recAy, GMAP, it);
        }

        if(this->recAzset){
            waves->recordData(this->recAz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            waves->writeSnap(it, waves->getPsnap());
        }

        if(this->snapAxset){ 
            waves->writeSnap(it, waves->getAxsnap());
        }

        if(this->snapAyset){ 
            waves->writeSnap(it, waves->getAysnap());
        }

        if(this->snapAzset){ 
            waves->writeSnap(it, waves->getAzsnap());
        }

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    result=MOD_OK;
    return result;
}

template<typename T>
ModellingAcoustic3D<T>::~ModellingAcoustic3D() {
    // Nothing here
}


// =============== ELASTIC 2D MODELLING CLASS =============== //
template<typename T>
ModellingElastic2D<T>::ModellingElastic2D(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recVxset = false;
    recVzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapVxset = false;
    snapVzset = false;
}

template<typename T>
ModellingElastic2D<T>::ModellingElastic2D(std::shared_ptr<ModelElastic2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recVxset = false;
    recVzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapVxset = false;
    snapVzset = false;
}

template<typename T>
int ModellingElastic2D<T>::run(){
     int result = MOD_ERR;
     int nt;
     float dt;
     float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     // Create the classes 
     std::shared_ptr<WavesElastic2D<T>> waves (new WavesElastic2D<T>(model, nt, dt, ot, this->getSnapinc()));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

    // Create snapshot for pressure
    if(this->snapPset){ 
        waves->createSnap(this->snapP, waves->getPsnap());
    }
    if(this->snapVxset){ 
        waves->createSnap(this->snapVx, waves->getVxsnap());
    }
    if(this->snapVzset){ 
        waves->createSnap(this->snapVz, waves->getVzsnap());
    }

    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source (Pressure)
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Pressure)
        if(this->recPset){
            waves->recordData(this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(this->recVx, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(this->recVz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            waves->writeSnap(it, waves->getPsnap());
        }

        if(this->snapVxset){ 
            waves->writeSnap(it, waves->getVxsnap());
        }

        if(this->snapVzset){ 
            waves->writeSnap(it, waves->getVzsnap());
        }
    }	
    
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingElastic2D<T>::~ModellingElastic2D() {
    // Nothing here
}

// =============== ELASTIC 3D MODELLING CLASS =============== //
template<typename T>
ModellingElastic3D<T>::ModellingElastic3D(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recVxset = false;
    recVyset = false;
    recVzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSyyset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapSyzset = false;
    snapSxyset = false;
    snapVxset = false;
    snapVyset = false;
    snapVzset = false;
}

template<typename T>
ModellingElastic3D<T>::ModellingElastic3D(std::shared_ptr<ModelElastic3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recVxset = false;
    recVyset = false;
    recVzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSyyset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapSyzset = false;
    snapSxyset = false;
    snapVxset = false;
    snapVyset = false;
    snapVzset = false;
}

template<typename T>
int ModellingElastic3D<T>::run(){
     int result = MOD_ERR;
     int nt;
     float dt;
     float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     // Create the classes 
     std::shared_ptr<WavesElastic3D<T>> waves (new WavesElastic3D<T>(model, nt, dt, ot, this->getSnapinc()));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

	// Create log file
	this->createLog("log.txt");

    // Create snapshot for pressure
    if(this->snapPset){ 
        waves->createSnap(this->snapP, waves->getPsnap());
    }
    if(this->snapVxset){ 
        waves->createSnap(this->snapVx, waves->getVxsnap());
    }
    if(this->snapVyset){ 
        waves->createSnap(this->snapVy, waves->getVysnap());
    }
    if(this->snapVzset){ 
        waves->createSnap(this->snapVz, waves->getVzsnap());
    }

    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source (Pressure)
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Pressure)
        if(this->recPset){
            waves->recordData(this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(this->recVx, GMAP, it);
        }

        // Recording data (Vy)
        if(this->recVyset){
            waves->recordData(this->recVy, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(this->recVz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            waves->writeSnap(it, waves->getPsnap());
        }

        if(this->snapVxset){ 
            waves->writeSnap(it, waves->getVxsnap());
        }

        if(this->snapVyset){ 
            waves->writeSnap(it, waves->getVysnap());
        }

        if(this->snapVzset){ 
            waves->writeSnap(it, waves->getVzsnap());
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingElastic3D<T>::~ModellingElastic3D() {
    // Nothing here
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Modelling<float>;
template class ModellingAcoustic2D<float>;
template class ModellingAcoustic3D<float>;
template class ModellingElastic2D<float>;
template class ModellingElastic3D<float>;

}
