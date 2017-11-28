// Include statements
#include "fwi.h"

namespace rockseis {

// =============== ABSTRACT FWI CLASS =============== //
template<typename T>
Fwi<T>::Fwi() {
	order = 4;
    snapinc=1;
    snapmethod = FULL;
    ncheck = 0;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
    misfit_type = DIFFERENCE;
    misfit = 0.0;
    noreverse = false;
}

template<typename T>
Fwi<T>::Fwi(int _order, int _snapinc) {
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
    misfit_type = DIFFERENCE;
    misfit = 0.0;
    noreverse = false;
}

template<typename T>
bool Fwi<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return FWI_ERR;
	}else{
		Flog.close();
		return FWI_OK;
	}
}

template<typename T>
void Fwi<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Fwi<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Fwi<T>::writeProgressbar(int x, int n, int r, int w){
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
void Fwi<T>::writeProgress(int x, int n, int r, int w){
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
Fwi<T>::~Fwi() {
    // Nothing here
}

// =============== ACOUSTIC 2D FWI CLASS =============== //

template<typename T>
FwiAcoustic2D<T>::FwiAcoustic2D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
}

template<typename T>
FwiAcoustic2D<T>::FwiAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataP, int order, int snapinc):Fwi<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    sourceset = true;
    dataPset = true;
    modelset = true;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
}

template<typename T>
void FwiAcoustic2D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T* Vp, T* Rho)
{
    if(!vpgradset && !rhogradset) rs_error("FwiAcoustic2D:crossCorr: No gradient set in fwi class");
    if(!vpgrad->getAllocated()) vpgrad->allocateImage();
    if(!rhograd->getAllocated()) rhograd->allocateImage();
    int ix, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T *rhograddata = rhograd->getImagedata();
    T mspx, mspz;
    T mrpx, mrpz;
    T mrxx, mrzz;
    T L;
    T vpscale, rhoscale1, rhoscale2;
    int nx = vpgrad->getNx();
    int nz = vpgrad->getNz();
    int dx = vpgrad->getDx();
    int dz = vpgrad->getDz();
    int nxs = nx+2*pads;
    int nxr = nx+2*padr;
    for (ix=1; ix<nx; ix++){
        {
            for (iz=1; iz<nz; iz++){
                L = Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)];
                vpscale = -2.0/(Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)]);
                rhoscale1 = -1.0/(Rho[km2D(ix, iz)]*Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)]);
                rhoscale2 = -1.0/(Rho[km2D(ix, iz)]*Rho[km2D(ix, iz)]);
                mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
                mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
                vpgraddata[ki2D(ix,iz)] -= vpscale*wsp[ks2D(ix+pads, iz+pads)]*L*(mrxx + mrzz);
                rhograddata[ki2D(ix,iz)] -= rhoscale1*wsp[ks2D(ix+pads, iz+pads)]*L*(mrxx + mrzz);
                mspx = (wsp[ks2D(ix+pads, iz+pads)] - wsp[ks2D(ix+pads-1, iz+pads)])/dx;
                mspz = (wsp[ks2D(ix+pads, iz+pads)] - wsp[ks2D(ix+pads, iz+pads-1)])/dz;
                mrpx = (wrp[kr2D(ix+padr, iz+padr)] - wrp[kr2D(ix+padr-1, iz+padr)])/dx;
                mrpz = (wrp[kr2D(ix+padr, iz+padr)] - wrp[kr2D(ix+padr, iz+padr-1)])/dz;
                rhograddata[ki2D(ix,iz)] -= rhoscale2*(mspx*mrpx + mspz*mrpz);
            }	
        }
    }
}

template<typename T>
void FwiAcoustic2D<T>::computeMisfit(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    T res = 0.0;
    T misfit = 0.0;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;

                for(it=0; it<nt; it++){
                    res = (-1.0)*(mod[I(it, itr)]*rec[I(it, itr)]/(norm1*norm2));
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += res;
                }
            }

            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
    }
    // Set the final misfit value
    this->setMisfit(misfit);
}

template<typename T>
void FwiAcoustic2D<T>::computeResiduals(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* res = dataresP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2, norm3;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                norm3 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                   norm3 += mod[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;
                norm3 /= (norm1*norm2);

                for(it=0; it<nt; it++){
                    res[I(it, itr)]=(-1.0)*((rec[I(it, itr)]/(norm1*norm2)) - (mod[I(it, itr)]/(norm1*norm1))*norm3);
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }

            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int FwiAcoustic2D<T>::run(){
     int result = FWI_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();
     if(!this->datamodPset || !this->dataresPset) rs_error("FwiAcoustic2D::run: datamodP and dataresP must be set before running the simulation.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

     this->writeLog("Running 2D Acoustic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->datamodPset){
            waves->recordData(this->datamodP, GMAP, it);
        }

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

    // Create images 
    vpgrad->allocateImage();
    rhograd->allocateImage();

    // Compute msifit
    computeMisfit();

    if(!this->getNoreverse()){
        // Compute residuals
        computeResiduals();

        Psnap->openSnap(this->getSnapfile(), 'r');
        Psnap->allocSnap(0);

        this->writeLog("\nDoing reverse-time Loop.");
        // Loop over reverse time
        for(int it=0; it < nt; it++)
        {
            // Time stepping
            waves->forwardstepAcceleration(model, der);
            waves->forwardstepStress(model, der);

            // Inserting residuals
            waves->insertSource(model, dataresP, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Psnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
                T *wrp = waves->getP1();
                T* wrx = waves->getAx(); 
                T* wrz = waves->getAz(); 
                crossCorr(Psnap->getData(0), 0, wrp, wrx, wrz, waves->getLpml(), model->getVp(), model->getR());
            }
            // Record wavelet gradient
            waves->recordData(this->wavgrad, SMAP, nt-1-it);

            // Roll the pointers P1 and P2
            waves->roll();

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");

    } // End of reverse loop

    //Remove snapshot file
    Psnap->removeSnap();

    result=FWI_OK;
    return result;
}

template<typename T>
int FwiAcoustic2D<T>::run_optimal(){
     int result = FWI_ERR;
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
     vpgrad->allocateImage();

     this->writeLog("Running 2D Acoustic full-waveform inversion gradient with optimal checkpointing.");
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

                // Recording data 
                if(this->datamodPset && !reverse){
                    waves_fw->recordData(this->datamodP, GMAP, it);
                }

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

            // Recording data 
            if(this->datamodPset){
                waves_fw->recordData(this->datamodP, GMAP, capo);
            }

            // Compute msifit
            computeMisfit();

            // Compute residuals
            computeResiduals();

            // Inserting data
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

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
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

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
     this->writeLog("\nGradient computation completed.");


	//Remove snapshot file
	optimal->removeCheck();

    result=FWI_OK;
    return result;
}

template<typename T>
FwiAcoustic2D<T>::~FwiAcoustic2D() {
    // Nothing here
}

// =============== ACOUSTIC 3D FWI CLASS =============== //

template<typename T>
FwiAcoustic3D<T>::FwiAcoustic3D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
}

template<typename T>
FwiAcoustic3D<T>::FwiAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataP, int order, int snapinc):Fwi<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    modelset = true;
    sourceset = true;
    dataPset = true;
    datamodPset = false;
    dataresPset = false;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    dataweightset = false;
}

template<typename T>
void FwiAcoustic3D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wry, T*wrz, int padr, T* Vp, T* Rho)
{
    if(!vpgradset && !rhogradset) rs_error("FwiAcoustic3D:crossCorr: No gradient set in fwi class");
	if(!vpgrad->getAllocated()) vpgrad->allocateImage();
    if(!rhograd->getAllocated()) rhograd->allocateImage();
	int ix, iy, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T *rhograddata = rhograd->getImagedata();
    T mspx, mspy, mspz;
    T mrpx, mrpy, mrpz;
    T mrxx, mryy, mrzz;
    T L;
    T vpscale, rhoscale1, rhoscale2;
	int nx = vpgrad->getNx();
	int ny = vpgrad->getNy();
	int nz = vpgrad->getNz();
    int dx = vpgrad->getDx();
    int dy = vpgrad->getDy();
    int dz = vpgrad->getDz();
	int nxs = nx + 2*pads;
	int nxr = nx + 2*padr;
	int nys = ny + 2*pads;
	int nyr = ny + 2*padr;
    for (ix=1; ix<nx; ix++){
        for (iy=1; iy<ny; iy++){
            for (iz=1; iz<nz; iz++){
                L = Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                vpscale = -2.0/(Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]);
                rhoscale1 = -1.0/(Rho[km3D(ix, iy, iz)]*Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]);
                rhoscale2 = -1.0/(Rho[km3D(ix, iy, iz)]*Rho[km3D(ix, iy, iz)]);
                mrxx = (wrx[kr3D(ix+padr, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mryy = (wry[kr3D(ix+padr, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrzz = (wrz[kr3D(ix+padr, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                mspx = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads-1, iy+pads, iz+pads)])/dx;
                mspy = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads, iy+pads-1, iz+pads)])/dy;
                mspz = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads, iy+pads, iz+pads-1)])/dz;
                mrpx = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mrpy = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrpz = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                vpgraddata[ki3D(ix,iy,iz)] -= vpscale*wsp[ks3D(ix+pads, iy+pads, iz+pads)]*L*(mrxx + mryy + mrzz);
                rhograddata[ki3D(ix,iy,iz)] -= rhoscale1*wsp[ks3D(ix+pads, iy+pads, iz+pads)]*L*(mrxx + mryy + mrzz);
                rhograddata[ki3D(ix,iy,iz)] -= rhoscale2*(mspx*mrpx + mspy*mrpy + mspz*mrpz);
            }	
        }
    }
}

template<typename T>
void FwiAcoustic3D<T>::computeResiduals(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* res = dataresP->getData();
    T *wei = NULL;
    if(dataweightset)
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2, norm3;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                norm3 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                   norm3 += mod[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;
                norm3 /= (norm1*norm2);

                for(it=0; it<nt; it++){
                    res[I(it, itr)]=(-1.0)*((rec[I(it, itr)]/(norm1*norm2)) - (mod[I(it, itr)]/(norm1*norm1))*norm3);
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }

            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int FwiAcoustic3D<T>::run(){
     int result = FWI_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

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

     this->writeLog("Running 3D Acoustic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->datamodPset){
            waves->recordData(this->datamodP, GMAP, it);
        }

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
    vpgrad->allocateImage();
    rhograd->allocateImage();

    if(!this->getNoreverse()){
        // Compute residuals
        computeResiduals();

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
            waves->insertSource(model, dataresP, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Psnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
                T *wrp = waves->getP1();
                T* wrx = waves->getAx(); 
                T* wry = waves->getAy(); 
                T* wrz = waves->getAz(); 
                crossCorr(Psnap->getData(0), 0, wrp, wrx, wry, wrz, waves->getLpml(), model->getVp(), model->getR());
            }

            // Record wavelet gradient
            waves->recordData(this->wavgrad, SMAP, nt-1-it);

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);

            // Roll the pointers P1 and P2
            waves->roll();

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");
    } // End of reverse loop

	//Remove snapshot file
	Psnap->removeSnap();

    result=FWI_OK;
    return result;
}

template<typename T>
int FwiAcoustic3D<T>::run_optimal(){
     int result = FWI_ERR;
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
     vpgrad->allocateImage();


     this->writeLog("Running 3D Acoustic full-waveform inversion gradient with optimal checkpointing.");
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

                // Recording data 
                if(this->datamodPset && !reverse){
                    waves_fw->recordData(this->datamodP, GMAP, it);
                }

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

            // Recording data 
            if(this->datamodPset){
                waves_fw->recordData(this->datamodP, GMAP, capo);
            }

            // Compute residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wry = waves_bw->getAy(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wry, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

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

            // Inserting residuals
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wry = waves_bw->getAy(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wry, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

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
     this->writeLog("\nGradient computation completed.");


	//Remove snapshot file
	optimal->removeCheck();

    result=FWI_OK;
    return result;
}

template<typename T>
FwiAcoustic3D<T>::~FwiAcoustic3D() {
    // Nothing here
}

// =============== ELASTIC 2D FWI CLASS =============== //

template<typename T>
FwiElastic2D<T>::FwiElastic2D(){
    sourceset = false;
    dataVxset = false;
    dataVzset = false;
    modelset = false;
    vpgradset = false;
    vsgradset = false;
    datamodVxset = false;
    datamodVzset = false;
    dataresVxset = false;
    dataresVzset = false;
    dataweightset = false;
}

template<typename T>
FwiElastic2D<T>::FwiElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataVx, std::shared_ptr<Data2D<T>> _dataVz, int order, int snapinc):Fwi<T>(order, snapinc){
    source = _source;
    dataVx = _dataVx;
    dataVz = _dataVz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataVxset = true;
    dataVzset = true;
    vpgradset = false;
    vsgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodVxset = false;
    datamodVzset = false;
    dataresVxset = false;
    dataresVzset = false;
    dataweightset = false;
}

template<typename T>
void FwiElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it)
{
    if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) rs_error("FwiElastic2D<T>::crossCorr: No gradient set for computation.");
	int ix, iz;
    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getVx();
    T* wrz = waves_bw->getVz();
    T* rsxx = waves_bw->getSxx();
    T* rszz = waves_bw->getSzz();
    T* rsxz = waves_bw->getSxz();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *wavgraddata = NULL;
	T *rhograddata = NULL;
	T msxx=0, mszz=0, msxz=0, mrxx=0, mrzz=0, mrxz=0, uderx=0, uderz=0;
	int nx;
	int nz;
	T dx;
	T dz;

    T* Rx = model->getRx();
    T* Rz = model->getRz();

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
	if(rhogradset){
		if(!rhograd->getAllocated()){
			rhograd->allocateImage();
		}
		rhograddata = rhograd->getImagedata();
	}

    int nt=0;
    int ntrace=0;
    int i=0;
    Point2D<int> *map=NULL;
    
	if(wavgradset){
		wavgraddata = wavgrad->getData();
        nt = wavgrad->getNt();
        ntrace = wavgrad->getNtrace();
        map = (wavgrad->getGeom())->getSmap();
	}

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

            if(vpgradset || rhogradset){
                vpgraddata[ki2D(ix,iz)] += (msxx + mszz) * (mrxx + mrzz);
            }

            if(vsgradset || rhogradset){
                msxz = 0.5*(wsx[ks2D(ix+pads, iz+pads+1)] - wsx[ks2D(ix+pads, iz+pads)])/dz;
                msxz += 0.5*(wsx[ks2D(ix+pads-1, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads-1)])/dz;
                msxz += 0.5*(wsz[ks2D(ix+pads+1, iz+pads)] - wsz[ks2D(ix+pads, iz+pads)])/dx;
                msxz += 0.5*(wsz[ks2D(ix+pads, iz+pads-1)] - wsz[ks2D(ix+pads-1, iz+pads-1)])/dx;

                mrxz = 0.5*(wrx[kr2D(ix+padr, iz+padr+1)] - wrx[kr2D(ix+padr, iz+padr)])/dz;
                mrxz += 0.5*(wrx[kr2D(ix+padr-1, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr-1)])/dz;
                mrxz += 0.5*(wrz[kr2D(ix+padr+1, iz+padr)] - wrz[kr2D(ix+padr, iz+padr)])/dx;
                mrxz += 0.5*(wrz[kr2D(ix+padr, iz+padr-1)] - wrz[kr2D(ix+padr-1, iz+padr-1)])/dx;

                vsgraddata[ki2D(ix,iz)] += (2.0*msxx*mrxx + 2.0*mszz*mrzz + msxz*mrxz);
            }
            if(wavgradset){
                switch(wavgrad->getField()){
                    case PRESSURE:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = (mrxx + mrzz);
                            }
                        }
                        break;
                    case VX:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = 0.5*(wrx[kr2D(ix+padr, iz+padr)] + wrx[kr2D(ix+padr-1, iz+padr)]);
                            }
                        }
                        break;
                    case VZ:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = 0.5*(wrz[kr2D(ix+padr, iz+padr)] + wrz[kr2D(ix+padr, iz+padr-1)]);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }

            if(rhogradset){
                uderx = 0.25*wsx[ks2D(ix+pads, iz+pads)]*Rx[kr2D(ix+padr, iz+padr)]*(rsxx[kr2D(ix+padr+1, iz+padr)] - rsxx[kr2D(ix+padr, iz+padr)])/dx;
                uderx += 0.25*wsx[ks2D(ix+pads-1, iz+pads)]*Rx[kr2D(ix+padr-1, iz+padr)]*(rsxx[kr2D(ix+padr, iz+padr)] - rsxx[kr2D(ix+padr-1, iz+padr)])/dx;
                uderx += 0.25*wsx[ks2D(ix+pads, iz+pads)]*Rx[kr2D(ix+padr, iz+padr)]*(rsxz[kr2D(ix+padr, iz+padr)] - rsxz[kr2D(ix+padr, iz+padr-1)])/dz;
                uderx += 0.25*wsx[ks2D(ix+pads-1, iz+pads)]*Rx[kr2D(ix+padr-1, iz+padr)]*(rsxz[kr2D(ix+padr-1, iz+padr)] - rsxz[kr2D(ix+padr-1, iz+padr-1)])/dz;

                uderz = 0.25*wsz[ks2D(ix+pads, iz+pads)]*Rz[kr2D(ix+padr, iz+padr)]*(rsxz[kr2D(ix+padr, iz+padr)] - rsxz[kr2D(ix+padr-1, iz+padr)])/dx;
                uderz += 0.25*wsz[ks2D(ix+pads, iz+pads-1)]*Rz[kr2D(ix+padr, iz+padr-1)]*(rsxz[kr2D(ix+padr, iz+padr-1)] - rsxz[kr2D(ix+padr-1, iz+padr-1)])/dx;
                uderz += 0.25*wsz[ks2D(ix+pads, iz+pads)]*Rz[kr2D(ix+padr, iz+padr)]*(rszz[kr2D(ix+padr, iz+padr+1)] - rszz[kr2D(ix+padr, iz+padr)])/dz;
                uderz += 0.25*wsz[ks2D(ix+pads, iz+pads-1)]*Rz[kr2D(ix+padr, iz+padr-1)]*(rszz[kr2D(ix+padr, iz+padr)] - rszz[kr2D(ix+padr, iz+padr-1)])/dz;
                rhograddata[ki2D(ix,iz)] -= (uderx + uderz);
            }
        }
    }	
}

template<typename T>
void FwiElastic2D<T>::scaleGrad(std::shared_ptr<ModelElastic2D<T>> model)
{
	int ix, iz;
    T* Vp = model->getVp();
    T* Vs = model->getVs();
    T* Rho = model->getR();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *rhograddata = NULL;
	T vpscale;
	T vsscale;
    T rhoscale1;
    T rhoscale2;
	int nx;
	int nz;

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
	if(rhogradset){
		if(!rhograd->getAllocated()){
			rhograd->allocateImage();
		}
		rhograddata = rhograd->getImagedata();
	}

        
    // Getting sizes
    nx = model->getNx();
    nz = model->getNz();

    T lambda = 0.0;
    T mu = 0.0;
    T rho = 0.0;

    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            vpscale = 2.0*Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)];
            vsscale = 2.0*Rho[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            rhoscale1 = Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)] -2.0*Vs[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            rhoscale2 = Vs[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            lambda = vpgraddata[ki2D(ix,iz)];
            if(vsgradset || rhogradset){
                mu = vsgraddata[ki2D(ix,iz)];
            }

            if(vpgradset){
                vpgraddata[ki2D(ix,iz)] = vpscale*lambda;
            }

            if(vsgradset){
                vsgraddata[ki2D(ix,iz)] = vsscale*mu -2.0*vsscale*lambda; 
            }

            if(rhogradset){
                rho = rhograddata[ki2D(ix,iz)];
                rhograddata[ki2D(ix,iz)] = rhoscale1*lambda + rhoscale2*mu + rho;
            }

        }
    }	
}

template<typename T>
void FwiElastic2D<T>::computeMisfit(){
    size_t ntr = datamodVx->getNtrace();
    if(dataVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodVx->getNt();
    if(dataVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodVx->getData();
    T* recx = dataVx->getData();

    T* modz = datamodVz->getData();
    T* recz = dataVz->getData();
    T resx = 0.0;
    T resz = 0.0;
    T *wei = NULL;
    if(dataweightset)
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    T misfit = 0.0;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightset)
                   {
                       resx *= wei[I(it, itr)];
                       resz *= wei[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resz*resz);
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2;
            T znorm1, znorm2;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                
                znorm1 = 0.0;
                znorm2 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;

                for(it=0; it<nt; it++){
                    resx=(-1.0)*(modx[I(it, itr)]*recx[I(it, itr)]/(xnorm1*xnorm2));
                    resz=((-1.0)*(modz[I(it, itr)]*recz[I(it, itr)]/(znorm1*znorm2)));
                   if(dataweightset)
                   {
                       resx *= wei[I(it, itr)];
                       resz *= wei[I(it, itr)];
                   }

                   misfit += (resx + resz);
                }
            }
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightset)
                   {
                       resx *= wei[I(it, itr)];
                       resz *= wei[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resz*resz);
                }
            }
            break;
    }

    // Set the final misfit value
    this->setMisfit(misfit);
}

template<typename T>
void FwiElastic2D<T>::computeResiduals(){
    size_t ntr = datamodVx->getNtrace();
    if(dataVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodVx->getNt();
    if(dataVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodVx->getData();
    T* recx = dataVx->getData();
    T* resx = dataresVx->getData();

    T* modz = datamodVz->getData();
    T* recz = dataVz->getData();
    T* resz = dataresVz->getData();
    T *wei = NULL;
    if(dataweightset)
    {
        wei = dataweight->getData();
    }
    T dt = dataVx->getDt();
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=1; it<nt; it++){
                   resx[I(it, itr)] = dt*(modx[I(it, itr)] - recx[I(it, itr)]) + resx[I(it-1, itr)];
                   resz[I(it, itr)] = dt*(modz[I(it, itr)] - recz[I(it, itr)]) + resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2, xnorm3;
            T znorm1, znorm2, znorm3;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                xnorm3 = 0.0;
                
                znorm1 = 0.0;
                znorm2 = 0.0;
                znorm3 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];
                    xnorm3 += modx[I(it, itr)]*recx[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                    znorm3 += modz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;
                xnorm3 /= (xnorm1*xnorm2);

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;
                znorm3 /= (znorm1*znorm2);

                for(it=1; it<nt; it++){
                    resx[I(it, itr)]=dt*((-1.0)*((recx[I(it, itr)]/(xnorm1*xnorm2)) - (modx[I(it, itr)]/(xnorm1*xnorm1))*xnorm3)) + resx[I(it-1, itr)];
                    resz[I(it, itr)]=dt*((-1.0)*((recz[I(it, itr)]/(znorm1*znorm2)) - (modz[I(it, itr)]/(znorm1*znorm1))*znorm3)) + resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=1; it<nt; it++){
                   resx[I(it, itr)] = dt*(modx[I(it, itr)] - recx[I(it, itr)]) + resx[I(it-1, itr)];
                   resz[I(it, itr)] = dt*(modz[I(it, itr)] - recz[I(it, itr)]) + resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int FwiElastic2D<T>::run(){
	int result = FWI_ERR;
	if(!vpgradset && !vsgradset && !wavgradset) {
		rs_warning("FwiElastic2D::run: No image set");
		return result;
	}
	int nt;
	float dt;
	float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D<T>> waves (new WavesElastic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Vxsnap;
     Vxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Vxsnap->openSnap(this->getSnapfile() + "-vx", 'w'); // Create a new snapshot file
     Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field

     std::shared_ptr<Snapshot2D<T>> Vzsnap;
     Vzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Vzsnap->openSnap(this->getSnapfile() + "-vz", 'w'); // Create a new snapshot file
     Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field

     this->writeLog("Running 2D Elastic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field
        Vxsnap->writeSnap(it);

        Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field
        Vzsnap->writeSnap(it);

    	// Time stepping
    	waves->forwardstepStress(model, der);
    	waves->forwardstepVelocity(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Vx)
        if(this->datamodVxset){
            waves->recordData(this->datamodVx, GMAP, it);
        }

        // Recording data (Vz)
        if(this->datamodVzset){
            waves->recordData(this->datamodVz, GMAP, it);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Vxsnap->closeSnap();
    Vzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic2D<T>>(model, nt, dt, ot);

    // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

    Vxsnap->openSnap(this->getSnapfile() + "-vx", 'r');
    Vxsnap->allocSnap(0);

    Vzsnap->openSnap(this->getSnapfile() + "-vz", 'r');
    Vzsnap->allocSnap(0);

    // Compute msifit
    computeMisfit();

    // Compute Residuals
    computeResiduals();

     this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping 
    	waves->forwardstepStress(model, der);
    	waves->forwardstepVelocity(model, der);

    	// Inserting residuals
    	waves->insertSource(model, dataresVx, GMAP, (nt - 1 - it));
    	waves->insertSource(model, dataresVz, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Vxsnap->readSnap(nt - 1 - it);
        Vzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        crossCorr(Vxsnap->getData(0), Vzsnap->getData(0), 0, waves, model, (nt - 1 - it));

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
     this->writeLog("\nGradient computation completed.");

    // Scale gradients
    this->scaleGrad(model);
    
	//Remove snapshot file
	Vxsnap->removeSnap();
	Vzsnap->removeSnap();

    result=FWI_OK;
    return result;
}

template<typename T>
int FwiElastic2D<T>::run_optimal(){
     int result = FWI_ERR;
     if(!vpgradset && !vsgradset && !wavgradset) {
         rs_warning("FwiElastic2D::run: No image set");
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
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

     this->writeLog("Running 2D Elastic full-waveform inversion gradient with optimal checkpointing.");
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
                waves_fw->forwardstepStress(model, der);
                waves_fw->forwardstepVelocity(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Recording data (Vx)
                if(this->datamodVxset && !reverse){
                    waves_fw->recordData(this->datamodVx, GMAP, it);
                }

                // Recording data (Vz)
                if(this->datamodVzset && !reverse){
                    waves_fw->recordData(this->datamodVz, GMAP, it);
                }

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }

        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepStress(model, der);
            waves_fw->forwardstepVelocity(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Recording data (Vx)
            if(this->datamodVxset){
                waves_fw->recordData(this->datamodVx, GMAP, capo);
            }

            // Recording data (Vz)
            if(this->datamodVzset){
                waves_fw->recordData(this->datamodVz, GMAP, capo);
            }

            // Compute msifit
            computeMisfit();

            // Compute Residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertSource(model, dataresVx, GMAP, capo);
            waves_bw->insertSource(model, dataresVz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getVx();
            T *wsz = waves_fw->getVz();
            crossCorr(wsx, wsz, waves_fw->getLpml(), waves_bw, model, capo);

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
            waves_bw->forwardstepStress(model, der);
            waves_bw->forwardstepVelocity(model, der);

            // Inserting residuals
            waves_bw->insertSource(model, dataresVx, GMAP, capo);
            waves_bw->insertSource(model, dataresVz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getVx();
            T *wsz = waves_fw->getVz();
            crossCorr(wsx, wsz, waves_fw->getLpml(), waves_bw, model, capo);

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
     this->writeLog("\nGradient computation completed.");
    
    // Scale gradients
    this->scaleGrad(model);

	//Remove snapshot file
	optimal->removeCheck();

    result=FWI_OK;
    return result;
}

template<typename T>
FwiElastic2D<T>::~FwiElastic2D() {
    // Nothing here
}

// =============== ELASTIC 3D FWI CLASS =============== //

template<typename T>
FwiElastic3D<T>::FwiElastic3D(){
    sourceset = false;
    dataVxset = false;
    dataVyset = false;
    dataVzset = false;
    modelset = false;
    rhogradset = false;
    vpgradset = false;
    vsgradset = false;
    wavgradset = false;
    datamodVxset = false;
    datamodVyset = false;
    datamodVzset = false;
    dataresVxset = false;
    dataresVyset = false;
    dataresVzset = false;
    dataweightset = false;
}

template<typename T>
FwiElastic3D<T>::FwiElastic3D(std::shared_ptr<ModelElastic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataVx, std::shared_ptr<Data3D<T>> _dataVy, std::shared_ptr<Data3D<T>> _dataVz, int order, int snapinc):Fwi<T>(order, snapinc){
    source = _source;
    dataVx = _dataVx;
    dataVy = _dataVy;
    dataVz = _dataVz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataVxset = true;
    dataVyset = true;
    dataVzset = true;
    rhogradset = false;
    vpgradset = false;
    vsgradset = false;
    wavgradset = false;
    datamodVxset = false;
    datamodVyset = false;
    datamodVzset = false;
    dataresVxset = false;
    dataresVyset = false;
    dataresVzset = false;
    dataweightset = false;
}

template<typename T>
void FwiElastic3D<T>::crossCorr(T *wsx, T*wsy, T *wsz, int pads, std::shared_ptr<WavesElastic3D<T>> waves_bw, std::shared_ptr<ModelElastic3D<T>> model, int it)
{
    if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) rs_error("FwiElastic3D<T>::crossCorr: No gradient set for computation.");
	int ix, iy, iz;

    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getVx();
    T* wry = waves_bw->getVy();
    T* wrz = waves_bw->getVz();
    T* rsxx = waves_bw->getSxx();
    T* rsyy = waves_bw->getSyy();
    T* rszz = waves_bw->getSzz();
    T* rsyz = waves_bw->getSyz();
    T* rsxz = waves_bw->getSxz();
    T* rsxy = waves_bw->getSxy();

    T* Vp = model->getVp();
    T* Vs = model->getVs();
    T* Rho = model->getR();
    T* Rx = model->getRx();
    T* Ry = model->getRy();
    T* Rz = model->getRz();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *wavgraddata = NULL;
	T *rhograddata = NULL;

	T msxx=0, msyy=0, mszz=0, msyz=0, msxz=0, msxy=0, mrxx=0, mryy=0, mrzz=0, mryz=0, mrxz=0, mrxy=0, uderx=0, udery=0, uderz=0;
	T vpscale;
	T vsscale;
    T rhoscale1;
    T rhoscale2;
	int nx;
	int ny;
	int nz;
	T dx;
	T dy;
	T dz;

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
    if(rhogradset){
        if(!rhograd->getAllocated()){
            rhograd->allocateImage();
        }
        rhograddata = rhograd->getImagedata();
    }


    int nt=0;
    int ntrace=0;
    int i=0;
    Point3D<int> *map=NULL;
    
	if(wavgradset){
		wavgraddata = wavgrad->getData();
        nt = wavgrad->getNt();
        ntrace = wavgrad->getNtrace();
        map = (wavgrad->getGeom())->getSmap();
	}

    // Getting sizes
    nx = waves_bw->getNx();
    ny = waves_bw->getNy();
    nz = waves_bw->getNz();
    dx = waves_bw->getDx(); 
    dy = waves_bw->getDy(); 
    dz = waves_bw->getDz(); 

    int nxs = nx + 2*pads;
    int nxr = nx + 2*padr;
    int nys = ny + 2*pads;
    int nyr = ny + 2*padr;

    for (ix=1; ix<nx-1; ix++){
        for (iy=1; iy<ny-1; iy++){
            for (iz=1; iz<nz-1; iz++){
                vpscale = 2.0*Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                vsscale = 2.0*Rho[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];
                rhoscale1 = Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                rhoscale2 = Vs[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];

                msxx = (wsx[ks3D(ix+pads, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads)])/dx;
                msyy = (wsy[ks3D(ix+pads, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads)])/dy;
                mszz = (wsz[ks3D(ix+pads, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads-1)])/dz;
                mrxx = (wrx[kr3D(ix+padr, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mryy = (wry[kr3D(ix+padr, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrzz = (wrz[kr3D(ix+padr, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;

                if(vpgradset){
                    vpgraddata[ki3D(ix,iy,iz)] += vpscale*(msxx + msyy + mszz) * (mrxx + mryy + mrzz);
                }

                if(vsgradset || rhogradset){
                    msyz = 0.25*(wsz[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msyz += 0.25*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dy;
                    msyz += 0.25*(wsy[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msyz += 0.25*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dz;

                    mryz = 0.25*(wrz[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mryz += 0.25*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dy;
                    mryz += 0.25*(wry[kr3D(ix+padr, iy+padr, iz+padr+1)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mryz += 0.25*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dz;

                    msxz = 0.25*(wsx[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msxz += 0.25*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dz;
                    msxz += 0.25*(wsz[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxz += 0.25*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dx;

                    mrxz = 0.25*(wrx[kr3D(ix+padr, iy+padr, iz+padr+1)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mrxz += 0.25*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dz;
                    mrxz += 0.25*(wrz[kr3D(ix+padr+1, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxz += 0.25*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dx;

                    msxy = 0.25*(wsx[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msxy += 0.25*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dy;
                    msxy += 0.25*(wsy[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxy += 0.25*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dx;

                    mrxy = 0.25*(wrx[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mrxy += 0.25*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dy;
                    mrxy += 0.25*(wry[kr3D(ix+padr+1, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxy += 0.25*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dx;
                }
                if(vsgradset){
                    vsgraddata[ki3D(ix,iy,iz)] += vsscale*(-2.0*(msyy*mrzz + mszz*mryy) -2.0*(msxx*mrzz + mszz*mrxx) -2.0*(msyy*mrxx + msxx*mryy) + msyz*mryz + msxz*mrxz + msxy*mrxy);
                }
                if(wavgradset){
                    switch(wavgrad->getField()){
                        case PRESSURE:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = (mrxx + mryy + mrzz);
                                }
                            }
                            break;
                        case VX:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wrx[kr3D(ix+padr, iy+padr, iz+padr)] + wrx[kr3D(ix+padr-1, iy+padr, iz+padr)]);
                                }
                            }
                            break;
                        case VY:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wry[kr3D(ix+padr, iy+padr, iz+padr)] + wry[kr3D(ix+padr, iy+padr-1, iz+padr)]);
                                }
                            }
                            break;
                        case VZ:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr)] + wrz[kr3D(ix+padr, iy+padr, iz+padr-1)]);
                                }
                            }
                            break;
                        default:
                            break;
                    }
                }

            if(rhogradset){
                uderx = (1.0/6.0)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxx[kr3D(ix+padr+1, iy+padr, iz+padr)] - rsxx[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                uderx += (1.0/6.0)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxx[kr3D(ix+padr, iy+padr, iz+padr)] - rsxx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                uderx += (1.0/6.0)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                uderx += (1.0/6.0)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr-1, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dy;
                uderx += (1.0/6.0)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                uderx += (1.0/6.0)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr-1, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dz;

                udery = (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                udery += (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr-1, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dx;
                udery += (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyy[kr3D(ix+padr, iy+padr+1, iz+padr)] - rsyy[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                udery += (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsyy[kr3D(ix+padr, iy+padr, iz+padr)] - rsyy[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                udery += (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                udery += (1.0/6.0)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr-1, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dz;

                uderz = (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                uderz += (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr-1)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dx;
                uderz += (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                uderz += (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr-1)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dy;
                uderz += (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rszz[kr3D(ix+padr, iy+padr, iz+padr+1)] - rszz[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                uderz += (1.0/6.0)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rszz[kr3D(ix+padr, iy+padr, iz+padr)] - rszz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;

                rhograddata[ki3D(ix,iy,iz)] -= (uderx + udery + uderz);
                rhograddata[ki3D(ix,iy,iz)] += rhoscale1*(msxx + msyy + mszz) * (mrxx + mryy+ mrzz);
                rhograddata[ki3D(ix,iy,iz)] += rhoscale2*(-2.0*(msyy*mrzz + mszz*mryy) -2.0*(msxx*mrzz + mszz*mrxx) -2.0*(msyy*mrxx + msxx*mryy) + msyz*mryz + msxz*mrxz + msxy*mrxy);
            }

            }
        }	
    }
}

template<typename T>
void FwiElastic3D<T>::computeResiduals(){
    size_t ntr = datamodVx->getNtrace();
    if(dataVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodVx->getNt();
    if(dataVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataVy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataVy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresVz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresVz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodVx->getData();
    T* recx = dataVx->getData();
    T* resx = dataresVx->getData();

    T* mody = datamodVy->getData();
    T* recy = dataVy->getData();
    T* resy = dataresVy->getData();

    T* modz = datamodVz->getData();
    T* recz = dataVz->getData();
    T* resz = dataresVz->getData();
    T *wei = NULL;
    if(dataweightset)
    {
        wei = dataweight->getData();
    }
    T dt = dataVx->getDt();
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = dt*(modx[I(it, itr)] - recx[I(it, itr)]) + resx[I(it-1, itr)];
                   resy[I(it, itr)] = dt*(mody[I(it, itr)] - recy[I(it, itr)]) + resy[I(it-1, itr)];
                   resz[I(it, itr)] = dt*(modz[I(it, itr)] - recz[I(it, itr)]) + resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resy[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2, xnorm3;
            T ynorm1, ynorm2, ynorm3;
            T znorm1, znorm2, znorm3;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                xnorm3 = 0.0;

                ynorm1 = 0.0;
                ynorm2 = 0.0;
                ynorm3 = 0.0;

                znorm1 = 0.0;
                znorm2 = 0.0;
                znorm3 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];
                    xnorm3 += modx[I(it, itr)]*recx[I(it, itr)];

                    ynorm1 += mody[I(it, itr)]*mody[I(it, itr)];
                    ynorm2 += recy[I(it, itr)]*recy[I(it, itr)];
                    ynorm3 += mody[I(it, itr)]*recy[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                    znorm3 += modz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;
                xnorm3 /= (xnorm1*xnorm2);

                ynorm1 = sqrt(ynorm1);
                ynorm2 = sqrt(ynorm2);
                if(ynorm1 ==0 ) ynorm1= 1.0;
                if(ynorm2 ==0 ) ynorm2= 1.0;
                ynorm3 /= (ynorm1*ynorm2);

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;
                znorm3 /= (znorm1*znorm2);

                for(it=1; it<nt; it++){
                    resx[I(it, itr)]=dt*((-1.0)*((recx[I(it, itr)]/(xnorm1*xnorm2)) - (modx[I(it, itr)]/(xnorm1*xnorm1))*xnorm3)) + resx[I(it-1, itr)];
                    resy[I(it, itr)]=dt*((-1.0)*((recy[I(it, itr)]/(ynorm1*ynorm2)) - (mody[I(it, itr)]/(ynorm1*ynorm1))*ynorm3)) + resy[I(it-1, itr)];
                    resz[I(it, itr)]=dt*((-1.0)*((recz[I(it, itr)]/(znorm1*znorm2)) - (modz[I(it, itr)]/(znorm1*znorm1))*znorm3)) +resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resy[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=1; it<nt; it++){
                   resx[I(it, itr)] = dt*(modx[I(it, itr)] - recx[I(it, itr)]) + resx[I(it-1, itr)];
                   resy[I(it, itr)] = dt*(mody[I(it, itr)] - recy[I(it, itr)]) + resy[I(it-1, itr)];
                   resz[I(it, itr)] = dt*(modz[I(it, itr)] - recz[I(it, itr)]) + resz[I(it-1, itr)];
                   if(dataweightset)
                   {
                       resx[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resy[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                       resz[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}


template<typename T>
int FwiElastic3D<T>::run(){
     int result = FWI_ERR;
     if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) {
         rs_warning("FwiElastic3D::run: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D<T>> waves (new WavesElastic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Vxsnap;
     Vxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Vxsnap->openSnap(this->getSnapfile() + "-vx", 'w'); // Create a new snapshot file
     Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field

     std::shared_ptr<Snapshot3D<T>> Vysnap;
     Vysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Vysnap->openSnap(this->getSnapfile() + "-vy", 'w'); // Create a new snapshot file
     Vysnap->setData(waves->getVy(), 0); //Set Vy as snap field

     std::shared_ptr<Snapshot3D<T>> Vzsnap;
     Vzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Vzsnap->openSnap(this->getSnapfile() + "-vz", 'w'); // Create a new snapshot file
     Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field

     this->writeLog("Running 3D Elastic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Vxsnap->setData(waves->getVx(), 0); //Set Vx as snap field
        Vxsnap->writeSnap(it);

        Vysnap->setData(waves->getVy(), 0); //Set Vy as snap field
        Vysnap->writeSnap(it);

        Vzsnap->setData(waves->getVz(), 0); //Set Vz as snap field
        Vzsnap->writeSnap(it);

    	// Time stepping
    	waves->forwardstepStress(model, der);
    	waves->forwardstepVelocity(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Vx)
        if(this->datamodVxset){
            waves->recordData(this->datamodVx, GMAP, it);
        }

        // Recording data (Vy)
        if(this->datamodVyset){
            waves->recordData(this->datamodVy, GMAP, it);
        }

        // Recording data (Vz)
        if(this->datamodVzset){
            waves->recordData(this->datamodVz, GMAP, it);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Vxsnap->closeSnap();
    Vysnap->closeSnap();
    Vzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic3D<T>>(model, nt, dt, ot);

    // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

    Vxsnap->openSnap(this->getSnapfile() + "-vx", 'r');
    Vxsnap->allocSnap(0);

    Vysnap->openSnap(this->getSnapfile() + "-vy", 'r');
    Vysnap->allocSnap(0);

    Vzsnap->openSnap(this->getSnapfile() + "-vz", 'r');
    Vzsnap->allocSnap(0);

    if(!this->getNoreverse()){
        // Compute Residuals
        computeResiduals();

        this->writeLog("\nDoing reverse-time Loop.");
        // Loop over reverse time
        for(int it=0; it < nt; it++)
        {
            // Time stepping 
            waves->forwardstepStress(model, der);
            waves->forwardstepVelocity(model, der);

            // Inserting residuals
            waves->insertSource(model, dataresVx, GMAP, (nt - 1 - it));
            waves->insertSource(model, dataresVy, GMAP, (nt - 1 - it));
            waves->insertSource(model, dataresVz, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Vxsnap->readSnap(nt - 1 - it);
            Vysnap->readSnap(nt - 1 - it);
            Vzsnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Vxsnap->getEnddiff()) % Vxsnap->getSnapinc()) == 0){
                crossCorr(Vxsnap->getData(0), Vysnap->getData(0), Vzsnap->getData(0), 0, waves, model, (nt - 1 - it));
            }

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");
    } // End of reverse loop
    
	//Remove snapshot file
	Vxsnap->removeSnap();
	Vysnap->removeSnap();
	Vzsnap->removeSnap();

    result=FWI_OK;
    return result;
}

template<typename T>
int FwiElastic3D<T>::run_optimal(){
     int result = FWI_ERR;
     if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) {
         rs_warning("FwiElastic3D::run: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D<T>> waves_fw (new WavesElastic3D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic3D<T>> waves_bw (new WavesElastic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), waves_fw->getNy_pml(), waves_fw->getNz_pml(), waves_fw->getDx(), waves_fw->getDy(), waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

     this->writeLog("Running 3D Elastic full-waveform inversion gradient with optimal checkpointing.");
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
                waves_fw->forwardstepStress(model, der);
                waves_fw->forwardstepVelocity(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Recording data (Vx)
                if(this->datamodVxset && !reverse){
                    waves_fw->recordData(this->datamodVx, GMAP, it);
                }

                // Recording data (Vy)
                if(this->datamodVyset && !reverse){
                    waves_fw->recordData(this->datamodVy, GMAP, it);
                }

                // Recording data (Vz)
                if(this->datamodVzset && !reverse){
                    waves_fw->recordData(this->datamodVz, GMAP, it);
                }

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepStress(model, der);
            waves_fw->forwardstepVelocity(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Recording data (Vx)
            if(this->datamodVxset){
                waves_fw->recordData(this->datamodVx, GMAP, capo);
            }

            // Recording data (Vy)
            if(this->datamodVyset){
                waves_fw->recordData(this->datamodVy, GMAP, capo);
            }

            // Recording data (Vz)
            if(this->datamodVzset){
                waves_fw->recordData(this->datamodVz, GMAP, capo);
            }

            // Compute Residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertSource(model, dataresVx, GMAP, capo);
            waves_bw->insertSource(model, dataresVy, GMAP, capo);
            waves_bw->insertSource(model, dataresVz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getVx();
            T *wsy = waves_fw->getVy();
            T *wsz = waves_fw->getVz();

            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), waves_bw, model, capo);

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
            waves_bw->forwardstepStress(model, der);
            waves_bw->forwardstepVelocity(model, der);

            // Inserting residuals
            waves_bw->insertSource(model, dataresVx, GMAP, capo);
            waves_bw->insertSource(model, dataresVy, GMAP, capo);
            waves_bw->insertSource(model, dataresVz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getVx();
            T *wsy = waves_fw->getVy();
            T *wsz = waves_fw->getVz();
            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), waves_bw, model, capo);

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
     this->writeLog("\nGradient computation completed.");

	//Remove snapshot file
	optimal->removeCheck();

    result=FWI_OK;
    return result;
}

template<typename T>
FwiElastic3D<T>::~FwiElastic3D() {
    // Nothing here
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Fwi<float>;
template class Fwi<double>;
template class FwiAcoustic2D<float>;
template class FwiAcoustic2D<double>;
template class FwiAcoustic3D<float>;
template class FwiAcoustic3D<double>;

template class FwiElastic2D<float>;
template class FwiElastic2D<double>;
template class FwiElastic3D<float>;
template class FwiElastic3D<double>;

}
