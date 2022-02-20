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
    recVxset = false;
    recVzset = false;
    snapPset = false;
    snapVxset = false;
    snapVzset = false;
}

template<typename T>
ModellingAcoustic2D<T>::ModellingAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recVxset = false;
    recVzset = false;
    snapPset = false;
    snapVxset = false;
    snapVzset = false;
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

   // Set PML constants 
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
      Psnap->setData(waves->getP(), 0); //Set Pressure as the field to snap
   }
   if(this->snapVxset){ 
      Vxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
      Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
      Vxsnap->setData(waves->getVx(), 0); //Set Vx as the field to snap
   }
   if(this->snapVzset){ 
      Vzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
      Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
      Vzsnap->setData(waves->getVz(), 0); //Set Vz as the field to snap
   }

   this->writeLog("Running 2D Acoustic modelling.");
   // Loop over time
   for(int it=0; it < nt; it++)
   {
      // Time stepping
      waves->forwardstepVelocity(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getP());
      }

      // Inserting source 
      waves->insertSource(model, source, SMAP, it);

      // Recording data 
      if(this->recPset){
         waves->recordData(this->recP, GMAP, it);
      }

      if(this->recVxset){
         waves->recordData(this->recVx, GMAP, it);
      }

      if(this->recVzset){
         waves->recordData(this->recVz, GMAP, it);
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
ModellingAcoustic2D<T>::~ModellingAcoustic2D() {
    // Nothing here
}

// =============== ACOUSTIC 3D MODELLING CLASS =============== //

template<typename T>
ModellingAcoustic3D<T>::ModellingAcoustic3D(){
    sourceset = false;
    modelset = false;
    recPset = false;
    recVxset = false;
    recVyset = false;
    recVzset = false;
    snapPset = false;
    snapVxset = false;
    snapVyset = false;
    snapVzset = false;
}

template<typename T>
ModellingAcoustic3D<T>::ModellingAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    model = _model;
    sourceset = true;
    modelset = true;
    recPset = false;
    recVxset = false;
    recVyset = false;
    recVzset = false;
    snapPset = false;
    snapVxset = false;
    snapVyset = false;
    snapVzset = false;
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

   // Set PML constants
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
      Psnap->setData(waves->getP(), 0); //Set Pressure as the field to snap
   }
   if(this->snapVxset){ 
      Vxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
      Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
      Vxsnap->setData(waves->getVx(), 0); //Set Vx as the field to snap
   }
   if(this->snapVyset){ 
      Vysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
      Vysnap->openSnap(this->snapVy, 'w'); // Create a new snapshot file
      Vysnap->setData(waves->getVy(), 0); //Set Vy as the field to snap
   }
   if(this->snapVzset){ 
      Vzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
      Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
      Vzsnap->setData(waves->getVz(), 0); //Set Vz as the field to snap
   }

   this->writeLog("Running 3D Acoustic modelling.");
   // Loop over time
   for(int it=0; it < nt; it++)
   {
      // Time stepping
      waves->forwardstepVelocity(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVy());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getP());
      }

      // Inserting source
      waves->insertSource(model, source, SMAP, it);

      // Recording data 
      if(this->recPset){
         waves->recordData(this->recP, GMAP, it);
      }

      if(this->recVxset){
         waves->recordData(this->recVx, GMAP, it);
      }

      if(this->recVyset){
         waves->recordData(this->recVy, GMAP, it);
      }

      if(this->recVzset){
         waves->recordData(this->recVz, GMAP, it);
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

   // Set PML constants
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
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getSxx());
         (model->getDomain())->shareEdges3D(waves->getSzz());
         (model->getDomain())->shareEdges3D(waves->getSxz());
      }

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

   // Set PML constants
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
         // Time stepping velocity
      waves->forwardstepVelocity(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVy());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
         // Time stepping stress
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getSxx());
         (model->getDomain())->shareEdges3D(waves->getSyy());
         (model->getDomain())->shareEdges3D(waves->getSzz());
         (model->getDomain())->shareEdges3D(waves->getSxz());
         (model->getDomain())->shareEdges3D(waves->getSyz());
         (model->getDomain())->shareEdges3D(waves->getSxy());
      }

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
       waves->forwardstepVelocity(model, der);
       if((model->getDomain()->getStatus())){
          (model->getDomain())->shareEdges3D(waves->getVx());
          (model->getDomain())->shareEdges3D(waves->getVz());
       }
       waves->forwardstepStress(model, der);
       if((model->getDomain()->getStatus())){
          (model->getDomain())->shareEdges3D(waves->getSxx());
          (model->getDomain())->shareEdges3D(waves->getSzz());
          (model->getDomain())->shareEdges3D(waves->getSxz());
          (model->getDomain())->shareEdges3D(waves->getMxx());
          (model->getDomain())->shareEdges3D(waves->getMzz());
          (model->getDomain())->shareEdges3D(waves->getMxz());
       }

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
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVy());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }


      // Time stepping stress
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getSxx());
         (model->getDomain())->shareEdges3D(waves->getSyy());
         (model->getDomain())->shareEdges3D(waves->getSzz());
         (model->getDomain())->shareEdges3D(waves->getSxz());
         (model->getDomain())->shareEdges3D(waves->getSyz());
         (model->getDomain())->shareEdges3D(waves->getSxy());

         (model->getDomain())->shareEdges3D(waves->getMxx());
         (model->getDomain())->shareEdges3D(waves->getMyy());
         (model->getDomain())->shareEdges3D(waves->getMzz());
         (model->getDomain())->shareEdges3D(waves->getMxz());
         (model->getDomain())->shareEdges3D(waves->getMyz());
         (model->getDomain())->shareEdges3D(waves->getMxy());
      }


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

// =============== ELASTIC 2D MODELLING CLASS =============== //
template<typename T>
ModellingVti2D<T>::ModellingVti2D(){
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
ModellingVti2D<T>::ModellingVti2D(std::shared_ptr<ModelVti2D<T>> _model,std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
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
T ModellingVti2D<T>::getVpmax(){
    T *c11 = model->getC11();
    T *c33 = model->getC33();
    T *rho = model->getRx();
    // Find maximum Vp
    T Vpmax;
    T Vp1;
    T Vp2;
    T Vp;
    Vp1=sqrt(c11[0]*rho[0]);
    Vp2=sqrt(c33[0]*rho[0]);
    Vpmax = MAX(Vp1,Vp2);
    size_t n=model->getNx()*model->getNz();
    for(size_t i=1; i<n; i++){
        Vp1=sqrt(c11[i]*rho[i]);
        Vp2=sqrt(c33[i]*rho[i]);
        Vp = MAX(Vp1,Vp2);
        if(Vp > Vpmax){
            Vpmax = Vp;
        }
    }
    return Vpmax;
}

template<typename T>
bool ModellingVti2D<T>::checkStability(){
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
int ModellingVti2D<T>::run(){
   int result = MOD_ERR;
   int nt;
   T dt;
   T ot;

   nt = source->getNt();
   dt = source->getDt();
   ot = source->getOt();

   if(!this->checkStability()) rs_error("ModellingVti2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");
   this->createLog(this->getLogfile());

   // Create the classes 
   std::shared_ptr<WavesVti2D<T>> waves (new WavesVti2D<T>(model, nt, dt, ot));
   std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

   // Set PML constants
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

   this->writeLog("Running 2D Vti modelling.");
   // Loop over time
   for(int it=0; it < nt; it++)
   {
      // Time stepping
      waves->forwardstepVelocity(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getSxx());
         (model->getDomain())->shareEdges3D(waves->getSzz());
         (model->getDomain())->shareEdges3D(waves->getSxz());
      }

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
ModellingVti2D<T>::~ModellingVti2D() {
    // Nothing here
}

// =============== ORTHOROMBIC 3D MODELLING CLASS =============== //
template<typename T>
ModellingOrtho3D<T>::ModellingOrtho3D(){
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
ModellingOrtho3D<T>::ModellingOrtho3D(std::shared_ptr<ModelOrtho3D<T>> _model,std::shared_ptr<Data3D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
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
T ModellingOrtho3D<T>::getVpmax(){
    T *c11 = model->getC11();
    T *c22 = model->getC22();
    T *c33 = model->getC33();
    T *rho = model->getRx();
    // Find maximum Vp
    T Vpmax;
    T Vp1;
    T Vp2;
    T Vp3;
    T Vp;
    Vp1=sqrt(c11[0]*rho[0]);
    Vp2=sqrt(c22[0]*rho[0]);
    Vp3=sqrt(c33[0]*rho[0]);
    Vpmax = MAX(Vp1,MAX(Vp2,Vp3));
    size_t n=model->getNx()*model->getNz();
    for(size_t i=1; i<n; i++){
        Vp1=sqrt(c11[i]*rho[i]);
        Vp2=sqrt(c22[i]*rho[i]);
        Vp3=sqrt(c33[i]*rho[i]);
        Vp = MAX(Vp1,MAX(Vp2,Vp3));
        if(Vp > Vpmax){
            Vpmax = Vp;
        }
    }
    return Vpmax;
}

template<typename T>
bool ModellingOrtho3D<T>::checkStability(){
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
int ModellingOrtho3D<T>::run(){
   int result = MOD_ERR;
   int nt;
   T dt;
   T ot;

   nt = source->getNt();
   dt = source->getDt();
   ot = source->getOt();

   if(!this->checkStability()) rs_error("ModellingOrtho3D::run: Wavelet sampling interval (dt) does not match the stability criteria.");

   // Create the classes 
   std::shared_ptr<WavesOrtho3D<T>> waves (new WavesOrtho3D<T>(model, nt, dt, ot));
   std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

   // Set PML constants
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

   this->writeLog("Running 3D Ortho modelling.");
   // Loop over time
   for(int it=0; it < nt; it++)
   {
         // Time stepping velocity
      waves->forwardstepVelocity(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getVx());
         (model->getDomain())->shareEdges3D(waves->getVy());
         (model->getDomain())->shareEdges3D(waves->getVz());
      }
         // Time stepping stress
      waves->forwardstepStress(model, der);
      if((model->getDomain()->getStatus())){
         (model->getDomain())->shareEdges3D(waves->getSxx());
         (model->getDomain())->shareEdges3D(waves->getSyy());
         (model->getDomain())->shareEdges3D(waves->getSzz());
         (model->getDomain())->shareEdges3D(waves->getSxz());
         (model->getDomain())->shareEdges3D(waves->getSyz());
         (model->getDomain())->shareEdges3D(waves->getSxy());
      }

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
ModellingOrtho3D<T>::~ModellingOrtho3D() {
    // Nothing here
}

// =============== POROELASTIC 2D MODELLING CLASS =============== //
template<typename T>
ModellingPoroelastic2D<T>::ModellingPoroelastic2D(){
    sourceset = false;
    acu_modelset = false;
    poro_modelset = false;
    recSzzset = false;
    recVxset = false;
    recVzset = false;
    recPset = false;
    recQxset = false;
    recQzset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapVxset = false;
    snapVzset = false;
    snapQxset = false;
    snapQzset = false;
}

template<typename T>
ModellingPoroelastic2D<T>::ModellingPoroelastic2D(std::shared_ptr<ModelAcoustic2D<T>> _acu_model, std::shared_ptr<ModelPoroelastic2D<T>> _poro_model, std::shared_ptr<Data2D<T>> _source, int order, int snapinc):Modelling<T>(order, snapinc){
    source = _source;
    acu_model = _acu_model;
    poro_model = _poro_model;
    sourceset = true;
    acu_modelset = true;
    poro_modelset = true;
    recSzzset = false;
    recVxset = false;
    recVzset = false;
    recQxset = false;
    recQzset = false;
    recPset = false;
    snapPset = false;
    snapSxxset = false;
    snapSzzset = false;
    snapSxzset = false;
    snapVxset = false;
    snapVzset = false;
    snapQxset = false;
    snapQzset = false;
}

template<typename T>
T ModellingPoroelastic2D<T>::getAcu_vpmax(){
    T *Vp = acu_model->getVp();
    // Find maximum Vp
    T Vpmax;
    Vpmax=Vp[0];
    size_t n=acu_model->getNx()*acu_model->getNz();
    for(size_t i=1; i<n; i++){
        if(Vp[i] > Vpmax){
            Vpmax = Vp[i];
        }
    }
    return Vpmax;
}

template<typename T>
T ModellingPoroelastic2D<T>::getPoro_vpmax(){
    T *LuM = poro_model->getLuM();
    T *Rho = poro_model->getRho_x();
    // Find maximum Vp
    T Vpmax;
    T Vp;
    Vpmax=sqrt(LuM[0]/Rho[0]);
    size_t n=poro_model->getNx()*poro_model->getNz();
    for(size_t i=1; i<n; i++){
        Vp = sqrt(LuM[i]/Rho[i]);
        if(Vp > Vpmax){
            Vpmax = Vp;
        }
    }
    return Vpmax;
}

template<typename T>
bool ModellingPoroelastic2D<T>::checkStability(){
    T Vpmax = this->getPoro_vpmax();
    T dx = poro_model->getDx();
    T dz = poro_model->getDz();
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
void ModellingPoroelastic2D<T>::Velocity_BC(T *wz, T *p1, T *vz, T *qz, T *p2, T *szz, int nx, int nz_acu, int lpml){
        int i;
        off_t  poro_index;
        off_t  acu_index;
        T wtmp, qtmp;
        for (i=0; i < nx; i++){
            acu_index = (nz_acu-lpml-1)*nx + i;
            poro_index = (lpml)*nx + i;
            wtmp = wz[acu_index];
            qtmp = qz[poro_index];
            wz[acu_index] = qz[poro_index] + vz[poro_index];
            //qz[poro_index] = wtmp - vz[poro_index];
            //vz[poro_index] = wtmp - qtmp;
        }
}

template<typename T>
void ModellingPoroelastic2D<T>::Stress_BC(T *wz, T *p1, T *vz, T *qz, T *p2, T *szz, int nx, int nz_acu, int lpml){
        int i;
        off_t  poro_index;
        off_t  acu_index;
        T ptmp;
        for (i=0; i < nx; i++){
            acu_index = (nz_acu-lpml-1)*nx + i;
            poro_index = (lpml)*nx + i;
            ptmp = p1[acu_index];
            p1[acu_index] = -szz[poro_index];
            p2[poro_index] = ptmp;
            szz[poro_index] = -ptmp;
        }
}

template<typename T>
int ModellingPoroelastic2D<T>::run(){
   int result = MOD_ERR;
   int nt;
   T dt;
   T ot;

   nt = source->getNt();
   dt = source->getDt();
   ot = source->getOt();

   if(!this->checkStability()) rs_error("ModellingPoroelastic2D::run: Wavelet sampling interval (dt) does not match the stability criteria.");
   this->createLog(this->getLogfile());

   // Create the classes 
   std::shared_ptr<WavesAcoustic2D<T>> acu_waves (new WavesAcoustic2D<T>(acu_model, nt, dt, ot));
   std::shared_ptr<Der<T>> acu_der (new Der<T>(acu_waves->getNx_pml(), 1, acu_waves->getNz_pml(), acu_waves->getDx(), 1.0, acu_waves->getDz(), this->getOrder()));
   std::shared_ptr<WavesPoroelastic2D<T>> poro_waves (new WavesPoroelastic2D<T>(poro_model, nt, dt, ot));
   std::shared_ptr<Der<T>> poro_der (new Der<T>(poro_waves->getNx_pml(), 1, poro_waves->getNz_pml(), poro_waves->getDx(), 1.0, poro_waves->getDz(), this->getOrder()));

   // Set PML constants
   (acu_waves->getPml())->setSmax((poro_model->getF0()*5*this->getAcu_vpmax())/(2*acu_waves->getLpml()));
   (acu_waves->getPml())->setKmax(2);
   (acu_waves->getPml())->setAmax(3.1415*poro_model->getF0());

   (poro_waves->getPml())->setSmax((poro_model->getF0()*5*this->getPoro_vpmax())/(2*poro_waves->getLpml()));
   (poro_waves->getPml())->setKmax(2);
   (poro_waves->getPml())->setAmax(3.1415*poro_model->getF0());
   (poro_waves->getPml())->computeABC();

   // Create snapshots
   std::shared_ptr<Snapshot2D<T>> Szzsnap;
   std::shared_ptr<Snapshot2D<T>> Vxsnap;
   std::shared_ptr<Snapshot2D<T>> Vzsnap;
   std::shared_ptr<Snapshot2D<T>> acu_Psnap;
   std::shared_ptr<Snapshot2D<T>> poro_Psnap;
   std::shared_ptr<Snapshot2D<T>> Qxsnap;
   std::shared_ptr<Snapshot2D<T>> Qzsnap;

   if(this->snapSzzset){ 
      Szzsnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
      Szzsnap->openSnap(this->snapSzz, 'w'); // Create a new snapshot file
      Szzsnap->setData(poro_waves->getSzz(), 0); //Set Stress as second field
   }
   if(this->snapVxset){ 
      Vxsnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
      Vxsnap->openSnap(this->snapVx, 'w'); // Create a new snapshot file
      Vxsnap->setData(poro_waves->getVx(), 0); //Set Vx as field to snap
   }
   if(this->snapVzset){ 
      Vzsnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
      Vzsnap->openSnap(this->snapVz, 'w'); // Create a new snapshot file
      Vzsnap->setData(poro_waves->getVz(), 0); //Set Vz as field to snap
   }
   if(this->snapPset){ 
      acu_Psnap = std::make_shared<Snapshot2D<T>>(acu_waves, this->getSnapinc());
      acu_Psnap->openSnap("Acu"+this->snapP, 'w'); // Create a new snapshot file
      acu_Psnap->setData(acu_waves->getP(), 0); //Set Pore pressure as field
      poro_Psnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
      poro_Psnap->openSnap("Poro"+this->snapP, 'w'); // Create a new snapshot file
      poro_Psnap->setData(poro_waves->getP(), 0); //Set Pore pressure as field
   }

   if(this->snapQxset){ 
       Qxsnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
       Qxsnap->openSnap(this->snapQx, 'w'); // Create a new snapshot file
       Qxsnap->setData(poro_waves->getQx(), 0); //Set Qx as field to snap
   }
   if(this->snapQzset){ 
       Qzsnap = std::make_shared<Snapshot2D<T>>(poro_waves, this->getSnapinc());
       Qzsnap->openSnap(this->snapQz, 'w'); // Create a new snapshot file
       Qzsnap->setData(poro_waves->getQz(), 0); //Set Qz as field to snap
   }

   this->writeLog("Running 2D Poroelastic modelling.");
   // Loop over time
   for(int it=0; it < nt; it++)
   {
      // Time stepping
      acu_waves->forwardstepVelocity(acu_model, acu_der);
      if((acu_model->getDomain()->getStatus())){
         (acu_model->getDomain())->shareEdges3D(acu_waves->getVx());
         (acu_model->getDomain())->shareEdges3D(acu_waves->getVz());
      }
      poro_waves->forwardstepVelocity(poro_model, poro_der);
      if((poro_model->getDomain()->getStatus())){
         (poro_model->getDomain())->shareEdges3D(poro_waves->getVx());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getQx());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getVz());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getQz());
      }

      acu_waves->forwardstepStress(acu_model, acu_der);
      if((acu_model->getDomain()->getStatus())){
         (acu_model->getDomain())->shareEdges3D(acu_waves->getP());
      }

      poro_waves->forwardstepStress(poro_model, poro_der);
      if((poro_model->getDomain()->getStatus())){
         (poro_model->getDomain())->shareEdges3D(poro_waves->getP());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getSxx());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getSzz());
         (poro_model->getDomain())->shareEdges3D(poro_waves->getSxz());
      }
      this->Stress_BC(acu_waves->getVz(), acu_waves->getP(), poro_waves->getVz(), poro_waves->getQz(), poro_waves->getP(), poro_waves->getSzz(), acu_waves->getNx_pml(), acu_waves->getNz_pml(), acu_waves->getLpml());

      // Inserting source 
      poro_waves->insertSource(poro_model, source, SMAP, it);

      // Recording data 
      if(this->recPset){
         poro_waves->recordData(poro_model,this->recP, GMAP, it);
      }

      // Recording data 
      if(this->recSzzset){
         poro_waves->recordData(poro_model,this->recSzz, GMAP, it);
      }

      // Recording data (Vx)
      if(this->recVxset){
         poro_waves->recordData(poro_model,this->recVx, GMAP, it);
      }

      // Recording data (Vz)
      if(this->recVzset){
         poro_waves->recordData(poro_model,this->recVz, GMAP, it);
      }

      // Recording data (Qx)
      if(this->recQxset){
         poro_waves->recordData(poro_model,this->recQx, GMAP, it);
      }

      // Recording data (Qz)
      if(this->recQzset){
         poro_waves->recordData(poro_model,this->recQz, GMAP, it);
      }

      //Writting out results to snapshot file
      if(this->snapPset){ 
         acu_Psnap->writeSnap(it);
         poro_Psnap->writeSnap(it);
      }

      if(this->snapSzzset){ 
         Szzsnap->writeSnap(it);
      }

      if(this->snapVxset){ 
         Vxsnap->writeSnap(it);
      }

      if(this->snapVzset){ 
         Vzsnap->writeSnap(it);
      }

      if(this->snapQxset){ 
         Qxsnap->writeSnap(it);
      }

      if(this->snapQzset){ 
         Qzsnap->writeSnap(it);
      }

      // Output progress to logfile
      this->writeProgress(it, nt-1, 20, 48);
   }	

   this->writeLog("Modelling is complete.");
   result=MOD_OK;
   return result;
}


template<typename T>
ModellingPoroelastic2D<T>::~ModellingPoroelastic2D() {
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

template class ModellingPoroelastic2D<float>;
template class ModellingPoroelastic2D<double>;

template class ModellingVti2D<float>;
template class ModellingVti2D<double>;

template class ModellingOrtho3D<float>;
template class ModellingOrtho3D<double>;

}
