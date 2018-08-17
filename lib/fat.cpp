// Include statements
#include "fat.h"

namespace rockseis {

    // Elastic: Derivate V Source gradients inside scaleGrad
// =============== ABSTRACT FAT CLASS =============== //
template<typename T>
Fat<T>::Fat() {
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
    misfit = 0.0;
}

template<typename T>
bool Fat<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return FAT_ERR;
	}else{
		Flog.close();
		return FAT_OK;
	}
}

template<typename T>
void Fat<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Fat<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Fat<T>::writeProgressbar(int x, int n, int r, int w){
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
void Fat<T>::writeProgress(int x, int n, int r, int w){
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
Fat<T>::~Fat() {
    // Nothing here
}

// =============== ACOUSTIC 2D FAT CLASS =============== //

template<typename T>
FatAcoustic2D<T>::FatAcoustic2D(){
    sourceset = false;
    Tdataset = false;
    modelset = false;
    vpgradset = false;
    Tmodset = false;
    Tresset = false;
    Tweightset = false;
}

template<typename T>
FatAcoustic2D<T>::FatAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _Tdata):Fat<T>(){
    source = _source;
    Tdata = _Tdata;
    model = _model;
    sourceset = true;
    Tdataset = true;
    modelset = true;
    vpgradset = false;
    Tmodset = false;
    Tresset = false;
    Tweightset = false;
}

template<typename T>
int FatAcoustic2D<T>::run()
{
     int result = FAT_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();
     if(!this->Tmodset || !this->Tresset) rs_error("FatAcoustic2D::run: Tmod and Tres must be set before running the simulation.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<RaysAcoustic2D<T>> rays (new RaysAcoustic2D<T>(model));

     this->writeLog("Running 2D Acoustic first arrival tomography gradient with fast sweeping method.");

     this->writeLog("Doing forward Loop.");
     // Inserting source 
     rays->insertSource(source, SMAP);

     // Running sweeping method
     rays->solve();

     // Recording data 
     if(this->Tmodset){
         rays->recordData(Tmod, GMAP);
     }

     // Compute misfit
     computeMisfit();

     /////// Run adjoint modelling
     this->writeLog("Doing adjoint Loop.");
     rays->createRecmask(Tres, GMAP);
     rays->insertResiduals(Tres, GMAP);
     rays->solve_adj();

     // Calculate gradient
     vpgrad->allocateImage();
     this->scaleGrad(model, rays->getLam(), vpgrad->getImagedata());
   
    result=FAT_OK;
    return result;

}

template<typename T>
void FatAcoustic2D<T>::scaleGrad(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, T *lam, T *grad) {
    int nx, nz;
    nx = model->getNx();
    nz = model->getNz();
    T * vp = model->getVp();
    for (int i=0; i<nx*nz; i++){
        grad[i] = -1.0*lam[i]/CUB(vp[i]);
    }
}

template<typename T>
void FatAcoustic2D<T>::computeMisfit()
{
    size_t ntr = Tmod->getNtrace();
    if(Tdata->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    size_t nt = Tmod->getNt();
    if(Tdata->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");

    
    T* mod = Tmod->getData();
    T* res = Tres->getData();
    T* obs = Tdata->getData();
    T* wei = NULL;
    if(Tweightset) 
    {
        wei = Tweight->getData();
    }
    size_t itr;
    T deltaT = 0.0;
    T misfit = 0.0;

    Index I(1, ntr);
    for(itr=0; itr<ntr; itr++){
        deltaT = mod[I(1, itr)] - obs[I(1, itr)];
        if(Tweightset)
        {
            deltaT *= wei[I(1, itr)];
        }
        misfit += 0.5*deltaT*deltaT;
        res[I(1, itr)] = deltaT;

        if(Tweightset)
        {
            res[I(1, itr)] *= wei[I(1, itr)];
        }
    }

    // Set the final misfit value
    this->setMisfit(misfit);
}


template<typename T>
FatAcoustic2D<T>::~FatAcoustic2D() {
    // Nothing here
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Fat<float>;
template class Fat<double>;
template class FatAcoustic2D<float>;
template class FatAcoustic2D<double>;
}
