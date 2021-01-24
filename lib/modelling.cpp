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
void Modelling<T>::writeLog(const char *text){
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
T ModellingAcoustic2D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingAcoustic2D<T>::checkStability(){
    T Vpmax = this->getVpmax(); 
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
int ModellingAcoustic2D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
	 T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingAcoustic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     this->createLog(this->getLogfile());

     // Create the finite difference modelling classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

    // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     std::shared_ptr<Snapshot2D<T>> Axsnap;
     std::shared_ptr<Snapshot2D<T>> Azsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as the field to snap
    }
    if(this->snapAxset){ 
        Axsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Axsnap->openSnap(this->snapAx, 'w'); // Create a new snapshot file
        Axsnap->setData(waves->getAx(), 0); //Set Ax as the field to snap
    }
    if(this->snapAzset){ 
        Azsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Azsnap->openSnap(this->snapAz, 'w'); // Create a new snapshot file
        Azsnap->setData(waves->getAz(), 0); //Set Az as the field to snap
    }

     this->writeLog("Running 2D Acoustic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
        if((model->getDomain()->getStatus())){
            if((model->getDomain())->getLow() > 0){
                (model->getDomain())->copyFromboundary(0, waves->getP2());
                (model->getDomain())->copyToboundary(0, waves->getP2());
            }

            if((model->getDomain())->getHigh() > 0){
                (model->getDomain())->copyFromboundary(1, waves->getP2());
                (model->getDomain())->copyToboundary(1, waves->getP2());
            }
        }
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
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
            Psnap->writeSnap(it);
        }

        if(this->snapAxset){ 
            Axsnap->writeSnap(it);
        }

        if(this->snapAzset){ 
            Azsnap->writeSnap(it);
        }

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	

    this->writeLog("Modelling is complete.");
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
T ModellingAcoustic3D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingAcoustic3D<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingAcoustic3D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingAcoustic3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

	// Create log file
    this->createLog(this->getLogfile());

    // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     std::shared_ptr<Snapshot3D<T>> Axsnap;
     std::shared_ptr<Snapshot3D<T>> Aysnap;
     std::shared_ptr<Snapshot3D<T>> Azsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as the field to snap
    }
    if(this->snapAxset){ 
        Axsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Axsnap->openSnap(this->snapAx, 'w'); // Create a new snapshot file
        Axsnap->setData(waves->getAx(), 0); //Set Ax as the field to snap
    }
    if(this->snapAyset){ 
        Aysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Aysnap->openSnap(this->snapAy, 'w'); // Create a new snapshot file
        Aysnap->setData(waves->getAy(), 0); //Set Ay as the field to snap
    }
    if(this->snapAzset){ 
        Azsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Azsnap->openSnap(this->snapAz, 'w'); // Create a new snapshot file
        Azsnap->setData(waves->getAz(), 0); //Set Az as the field to snap
    }

     this->writeLog("Running 3D Acoustic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
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
            Psnap->writeSnap(it);
        }

        if(this->snapAxset){ 
            Axsnap->writeSnap(it);
        }

        if(this->snapAyset){ 
            Aysnap->writeSnap(it);
        }

        if(this->snapAzset){ 
            Azsnap->writeSnap(it);
        }

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
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
T ModellingElastic2D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingElastic2D<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingElastic2D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingElastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D<T>> waves (new WavesElastic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

    // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     std::shared_ptr<Snapshot2D<T>> Vxsnap;
     std::shared_ptr<Snapshot2D<T>> Vzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSzz(), 1); //Set Stress as second field
    }
    if(this->snapVxset){ 
        Vxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as field to snap
    }
    if(this->snapVzset){ 
        Vzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
        Vzsnap->setData(waves->getVz(), 0); //Set Vz as field to snap
    }

     this->writeLog("Running 2D Elastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->recPset){
            waves->recordData(model,this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(model,this->recVx, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(model,this->recVz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapVxset){ 
            Vxsnap->writeSnap(it);
        }

        if(this->snapVzset){ 
            Vzsnap->writeSnap(it);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
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
T ModellingElastic3D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingElastic3D<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingElastic3D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingElastic3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     // Create the classes 
     std::shared_ptr<WavesElastic3D<T>> waves (new WavesElastic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

	// Create log file
     this->createLog(this->getLogfile());

    // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     std::shared_ptr<Snapshot3D<T>> Vxsnap;
     std::shared_ptr<Snapshot3D<T>> Vysnap;
     std::shared_ptr<Snapshot3D<T>> Vzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSyy(), 1); //Set Stress as second field 
        Psnap->setData(waves->getSzz(), 2); //Set Stress as third field
    }
    if(this->snapVxset){ 
        Vxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as field to snap
    }
    if(this->snapVyset){ 
        Vysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vysnap->openSnap(this->snapVy, 'w'); // Create a new snapshot file
        Vysnap->setData(waves->getVy(), 0); //Set Vy as field to snap
    }
    if(this->snapVzset){ 
        Vzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
        Vzsnap->setData(waves->getVz(), 0); //Set Vz as field to snap
    }

     this->writeLog("Running 3D Elastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->recPset){
            waves->recordData(model, this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(model, this->recVx, GMAP, it);
        }

        // Recording data (Vy)
        if(this->recVyset){
            waves->recordData(model, this->recVy, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(model, this->recVz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapVxset){ 
            Vxsnap->writeSnap(it);
        }

        if(this->snapVyset){ 
            Vysnap->writeSnap(it);
        }

        if(this->snapVzset){ 
            Vzsnap->writeSnap(it);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingElastic3D<T>::~ModellingElastic3D() {
    // Nothing here
}

// =============== ELASTIC 2D MODELLING CLASS =============== //
template<typename T>
ModellingElastic2D_DS<T>::ModellingElastic2D_DS(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recUxset = false;
    recUzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapUxset = false;
    snapUzset = false;
}

template<typename T>
ModellingElastic2D_DS<T>::ModellingElastic2D_DS(std::shared_ptr<ModelElastic2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recUxset = false;
    recUzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapUxset = false;
    snapUzset = false;
}

template<typename T>
T ModellingElastic2D_DS<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingElastic2D_DS<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingElastic2D_DS<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingElastic2D_DS::run: Wavelet sampling interval (dt) does not match the stability criteria.");
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

    // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     std::shared_ptr<Snapshot2D<T>> Uxsnap;
     std::shared_ptr<Snapshot2D<T>> Uzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSzz(), 1); //Set Stress as second field
    }
    if(this->snapUxset){ 
        Uxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Uxsnap->openSnap(this->snapUx, 'w'); // Create a new snapshot file
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as field to snap
    }
    if(this->snapUzset){ 
        Uzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Uzsnap->openSnap(this->snapUz, 'w'); // Create a new snapshot file
        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as field to snap
    }

     this->writeLog("Running 2D Elastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

        // Inserting pressure source 
        waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
        // Inserting force source 
        waves->insertForcesource(model, source, SMAP, it);

        // Recording data 
        if(this->recPset){
            waves->recordData(model, this->recP, GMAP, it);
        }

        // Recording data (Ux)
        if(this->recUxset){
            waves->recordData(model, this->recUx, GMAP, it);
        }

        // Recording data (Uz)
        if(this->recUzset){
            waves->recordData(model, this->recUz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapUxset){ 
            Uxsnap->writeSnap(it);
        }

        if(this->snapUzset){ 
            Uzsnap->writeSnap(it);
        }

    	// Roll the pointers 
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingElastic2D_DS<T>::~ModellingElastic2D_DS() {
    // Nothing here
}

// =============== ELASTIC 3D DISPLACEMENT STRESS MODELLING CLASS =============== //
template<typename T>
ModellingElastic3D_DS<T>::ModellingElastic3D_DS(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recUxset = false;
    recUyset = false;
    recUzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSyyset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapSyzset = false;
    snapSxyset = false;
    snapUxset = false;
    snapUyset = false;
    snapUzset = false;
}

template<typename T>
ModellingElastic3D_DS<T>::ModellingElastic3D_DS(std::shared_ptr<ModelElastic3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recUxset = false;
    recUyset = false;
    recUzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSyyset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapSyzset = false;
    snapSxyset = false;
    snapUxset = false;
    snapUyset = false;
    snapUzset = false;
}

template<typename T>
T ModellingElastic3D_DS<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingElastic3D_DS<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingElastic3D_DS<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingElastic3D_DS::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     // Create the classes 
     std::shared_ptr<WavesElastic3D_DS<T>> waves (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

	// Create log file
     this->createLog(this->getLogfile());

    // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     std::shared_ptr<Snapshot3D<T>> Uxsnap;
     std::shared_ptr<Snapshot3D<T>> Uysnap;
     std::shared_ptr<Snapshot3D<T>> Uzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSyy(), 1); //Set Stress as second field 
        Psnap->setData(waves->getSzz(), 2); //Set Stress as third field
    }
    if(this->snapUxset){ 
        Uxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Uxsnap->openSnap(this->snapUx, 'w'); // Create a new snapshot file
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as field to snap
    }
    if(this->snapUyset){ 
        Uysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Uysnap->openSnap(this->snapUy, 'w'); // Create a new snapshot file
        Uysnap->setData(waves->getUy1(), 0); //Set Uy as field to snap
    }
    if(this->snapUzset){ 
        Uzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Uzsnap->openSnap(this->snapUz, 'w'); // Create a new snapshot file
        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as field to snap
    }

     this->writeLog("Running 3D Elastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting pressure source 
    	waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting force source 
    	waves->insertForcesource(model, source, SMAP, it);

        // Recording data 
        if(this->recPset){
            waves->recordData(model, this->recP, GMAP, it);
        }

        // Recording data (Ux)
        if(this->recUxset){
            waves->recordData(model, this->recUx, GMAP, it);
        }

        // Recording data (Uy)
        if(this->recUyset){
            waves->recordData(model, this->recUy, GMAP, it);
        }

        // Recording data (Uz)
        if(this->recUzset){
            waves->recordData(model, this->recUz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapUxset){ 
            Uxsnap->writeSnap(it);
        }

        if(this->snapUyset){ 
            Uysnap->writeSnap(it);
        }

        if(this->snapUzset){ 
            Uzsnap->writeSnap(it);
        }
        
        //Roll the pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingElastic3D_DS<T>::~ModellingElastic3D_DS() {
    // Nothing here
}


// =============== VISCOELASTIC 2D MODELLING CLASS =============== //
template<typename T>
ModellingViscoelastic2D<T>::ModellingViscoelastic2D(){
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
ModellingViscoelastic2D<T>::ModellingViscoelastic2D(std::shared_ptr<ModelViscoelastic2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
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
T ModellingViscoelastic2D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingViscoelastic2D<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingViscoelastic2D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingViscoelastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesViscoelastic2D<T>> waves (new WavesViscoelastic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

    // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     std::shared_ptr<Snapshot2D<T>> Vxsnap;
     std::shared_ptr<Snapshot2D<T>> Vzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSzz(), 1); //Set Stress as second field
    }
    if(this->snapVxset){ 
        Vxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as field to snap
    }
    if(this->snapVzset){ 
        Vzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
        Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
        Vzsnap->setData(waves->getVz(), 0); //Set Vz as field to snap
    }

     this->writeLog("Running 2D Viscoelastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepStress(model, der);
    	waves->forwardstepVelocity(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->recPset){
            waves->recordData(model,this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(model,this->recVx, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(model,this->recVz, GMAP, it);
        }
    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapVxset){ 
            Vxsnap->writeSnap(it);
        }

        if(this->snapVzset){ 
            Vzsnap->writeSnap(it);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingViscoelastic2D<T>::~ModellingViscoelastic2D() {
    // Nothing here
}

// =============== VISCOELASTIC 3D MODELLING CLASS =============== //
template<typename T>
ModellingViscoelastic3D<T>::ModellingViscoelastic3D(){
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
ModellingViscoelastic3D<T>::ModellingViscoelastic3D(std::shared_ptr<ModelViscoelastic3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
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
T ModellingViscoelastic3D<T>::getVpmax(){
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
    return Vpmax;
}

template<typename T>
bool ModellingViscoelastic3D<T>::checkStability(){
    T Vpmax = this->getVpmax();
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
int ModellingViscoelastic3D<T>::run(){
     int result = MOD_ERR;
     int nt;
     T dt;
     T ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     if(!this->checkStability()) rs_error("ModellingViscoelastic3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

     // Create the classes 
     std::shared_ptr<WavesViscoelastic3D<T>> waves (new WavesViscoelastic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     (waves->getPml())->setSmax(-this->getVpmax()*4*log(1e-6)/(2*waves->getLpml()*waves->getDx()));
     (waves->getPml())->setSmax(SMAX);
     (waves->getPml())->computeABC();

	// Create log file
     this->createLog(this->getLogfile());

    // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     std::shared_ptr<Snapshot3D<T>> Vxsnap;
     std::shared_ptr<Snapshot3D<T>> Vysnap;
     std::shared_ptr<Snapshot3D<T>> Vzsnap;
    if(this->snapPset){ 
        Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Psnap->openSnap(this->snapP, 'w'); // Create a new snapshot file
        Psnap->setData(waves->getSxx(), 0); //Set Stress as first field 
        Psnap->setData(waves->getSyy(), 1); //Set Stress as second field 
        Psnap->setData(waves->getSzz(), 2); //Set Stress as third field
    }
    if(this->snapVxset){ 
        Vxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as field to snap
    }
    if(this->snapVyset){ 
        Vysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vysnap->openSnap(this->snapVy, 'w'); // Create a new snapshot file
        Vysnap->setData(waves->getVy(), 0); //Set Vy as field to snap
    }
    if(this->snapVzset){ 
        Vzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
        Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
        Vzsnap->setData(waves->getVz(), 0); //Set Vz as field to snap
    }

     this->writeLog("Running 3D Viscoelastic modelling.");
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping displacement
    	waves->forwardstepVelocity(model, der);

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting pressure source 
    	waves->insertSource(model, source, SMAP, it);
    
        // Recording data 
        if(this->recPset){
            waves->recordData(model, this->recP, GMAP, it);
        }

        // Recording data (Vx)
        if(this->recVxset){
            waves->recordData(model, this->recVx, GMAP, it);
        }

        // Recording data (Vy)
        if(this->recVyset){
            waves->recordData(model, this->recVy, GMAP, it);
        }

        // Recording data (Vz)
        if(this->recVzset){
            waves->recordData(model, this->recVz, GMAP, it);
        }

    
    	//Writting out results to snapshot file
        if(this->snapPset){ 
            Psnap->writeSnap(it);
        }

        if(this->snapVxset){ 
            Vxsnap->writeSnap(it);
        }

        if(this->snapVyset){ 
            Vysnap->writeSnap(it);
        }

        if(this->snapVzset){ 
            Vzsnap->writeSnap(it);
        }
        
        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }	
    
    this->writeLog("Modelling is complete.");
    result=MOD_OK;
    return result;
}


template<typename T>
ModellingViscoelastic3D<T>::~ModellingViscoelastic3D() {
    // Nothing here
}




// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Modelling<float>;
template class Modelling<double>;
template class ModellingAcoustic2D<float>;
template class ModellingAcoustic2D<double>;
template class ModellingAcoustic3D<float>;
template class ModellingAcoustic3D<double>;
template class ModellingElastic2D<float>;
template class ModellingElastic2D<double>;
template class ModellingElastic3D<float>;
template class ModellingElastic3D<double>;
template class ModellingElastic2D_DS<float>;
template class ModellingElastic2D_DS<double>;

template class ModellingElastic3D_DS<float>;
template class ModellingElastic3D_DS<double>;

template class ModellingViscoelastic2D<float>;
template class ModellingViscoelastic2D<double>;

template class ModellingViscoelastic3D<float>;
template class ModellingViscoelastic3D<double>;


}
