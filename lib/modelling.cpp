// Include statements
#include "modelling.h"

namespace rockseis {

// =============== ABSTRACT MODELLING CLASS =============== //
template<typename T>
Modelling<T>::Modelling() {
	order = 4;
    snapinc=1;
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
