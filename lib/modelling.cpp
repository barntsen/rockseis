// Include statements
#include "modelling.h"

namespace rockseis {

// =============== ABSTRACT MODELLING CLASS =============== //
template<typename T>
Modelling<T>::Modelling() {
	order = 4;
}

template<typename T>
Modelling<T>::Modelling(int _order) {
	if(_order > 1 && _order < 9)
	{
		order = _order;
	}else{
		order = 4;
	}
}

template<typename T>
int Modelling<T>::Acoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model,std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> recP){
     int result = MOD_ERR;
     int nt;
     float dt;
     float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     int nx = model->getNx();
     int nz = model->getNz();
     float dx = (float) model->getDx();
     float dz = (float) model->getDz();
     float ox = (float) model->getOx();
     float oz = (float) model->getOz();

     // Other parameters
     int lpml = model->getLpml();

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(nx+2*lpml, 1, nz+2*lpml, dx, 1.0, dz, this->order));


    // Output snapshots to a binary file
    std::fstream myFile;
    myFile.open ("data.bin", std::ios::out | std::ios::binary);
    
    float *Szz;
    Szz = waves->getP2();
    // Loop over time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source (Pressure)
    	waves->insertSource(model, source, 2, 0, it);

    	// Recording data (Pressure)
    	waves->recordData(recP, 0, 1, it);
    
    	//Writting out results to binary file
    	myFile.write (reinterpret_cast<char *> (Szz), (nx+2*lpml) * (nz+2*lpml) * sizeof(float));
    	
    	// Roll the pointers P1 and P2
    	waves->roll();
    }	
    myFile.close();
    
    result=MOD_OK;
    return result;
}

template<typename T>
Modelling<T>::~Modelling() {
    // Nothing here
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Modelling<float>;

}
