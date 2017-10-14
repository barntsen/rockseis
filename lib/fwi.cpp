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
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
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
                }
            }

            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
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
}

template<typename T>
FwiAcoustic3D<T>::FwiAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataP, int order, int snapinc):Fwi<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    sourceset = true;
    dataPset = true;
    modelset = true;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
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
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
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
                }
            }

            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
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
}

template<typename T>
void FwiElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads, T* wrx, T* wrz, int padr, T* Vp, T* Vs, T* Rho, int it)
{
	int ix, iz;
	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *wavgraddata = NULL;
	//T *rhograddata = NULL;
	T msxx, mszz, msxz, mrxx, mrzz, mrxz;
	T vpscale;
	T vsscale;
	int nx;
	T dx;
	T dz;
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

    int nt;
    int ntrace;
    int i;
    Point2D<int> *map;
    
	if(wavgradset){
		wavgraddata = wavgrad->getData();
        nt = wavgrad->getNt();
        ntrace = wavgrad->getNtrace();
        map = (wavgrad->getGeom())->getSmap();
	}

	// Getting sizes
	if(vpgradset) {
		nx = vpgrad->getNx();
		nz = vpgrad->getNz();
		dx = vpgrad->getDx(); 
		dz = vpgrad->getDz(); 
	}else if(vsgradset){
		nx = vsgrad->getNx();
		nz = vsgrad->getNz();
		dx = vsgrad->getDx(); 
		dz = vsgrad->getDz(); 
    }else{
		nx = rhograd->getNx();
		nz = rhograd->getNz();
		dx = rhograd->getDx(); 
		dz = rhograd->getDz(); 
	}

	int nxs = nx+2*pads;
	int nxr = nx+2*padr;


    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            vpscale = -2.0*Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)];
            vsscale = 2.0*Rho[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            msxx = (wsx[ks2D(ix+pads, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads)])/dx;
            mszz = (wsz[ks2D(ix+pads, iz+pads)] - wsz[ks2D(ix+pads, iz+pads-1)])/dz;
            mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
            mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;

            if(vpgradset){
                vpgraddata[ki2D(ix,iz)] += vpscale*(msxx + mszz) * (mrxx + mrzz);
            }

            if(vsgradset){
                msxz = 0.5*(wsx[ks2D(ix+pads, iz+pads+1)] - wsx[ks2D(ix+pads, iz+pads)])/dz;
                msxz += 0.5*(wsx[ks2D(ix+pads-1, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads-1)])/dz;
                msxz += 0.5*(wsz[ks2D(ix+pads+1, iz+pads)] - wsz[ks2D(ix+pads, iz+pads)])/dx;
                msxz += 0.5*(wsz[ks2D(ix+pads, iz+pads-1)] - wsz[ks2D(ix+pads-1, iz+pads-1)])/dx;

                mrxz = 0.5*(wrx[kr2D(ix+padr, iz+padr+1)] - wrx[kr2D(ix+padr, iz+padr)])/dz;
                mrxz += 0.5*(wrx[kr2D(ix+padr-1, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr-1)])/dz;
                mrxz += 0.5*(wrz[kr2D(ix+padr+1, iz+padr)] - wrz[kr2D(ix+padr, iz+padr)])/dx;
                mrxz += 0.5*(wrz[kr2D(ix+padr, iz+padr-1)] - wrz[kr2D(ix+padr-1, iz+padr-1)])/dx;
                vsgraddata[ki2D(ix,iz)] += vsscale*(-2.0*msxx*mrzz + -2.0*mszz*mrxx + msxz*mrxz);
            }
            if(wavgradset){
                for (i=0; i < ntrace; i++) 
                {
                    if((map[i].x == ix) && (map[i].y == iz))
                    {
                        //wavgraddata[kwav(it,i)] = (mrxx + mrzz);
                    }
                }
            }
        }
    }	
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
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = modx[I(it, itr)] - recx[I(it, itr)];
                   resz[I(it, itr)] = modz[I(it, itr)] - recz[I(it, itr)];
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


                for(it=0; it<nt; it++){
                    resx[I(it, itr)]=(-1.0)*((recx[I(it, itr)]/(xnorm1*xnorm2)) - (modx[I(it, itr)]/(xnorm1*xnorm1))*xnorm3);
                    resz[I(it, itr)]=(-1.0)*((recz[I(it, itr)]/(znorm1*znorm2)) - (modz[I(it, itr)]/(znorm1*znorm1))*znorm3);
                }
            }
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = modx[I(it, itr)] - recx[I(it, itr)];
                   resz[I(it, itr)] = modz[I(it, itr)] - recz[I(it, itr)];
                }
            }
            break;
    }
}

template<typename T>
int FwiElastic2D<T>::run(){
	int result = FWI_ERR;
	if(!vpgradset && !vsgradset) {
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
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data (Vx)
        if(this->datamodVxset){
            waves->recordData(this->datamodVx, GMAP, it);
        }

        // Recording data (Vy)
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

    // Compute Residuals
    computeResiduals();

    // Get models for scaling
    T *Vp, *Vs, *Rho;
    Vp = model->getVp();
    Vs = model->getVs();
    Rho = model->getR();

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
        if((((nt - 1 - it)-Vxsnap->getEnddiff()) % Vxsnap->getSnapinc()) == 0){
            T *Vxr = waves->getVx();
            T *Vzr = waves->getVz();
            crossCorr(Vxsnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vzr, waves->getLpml(), Vp, Vs, Rho, (nt - 1 - it));
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
	//Remove snapshot file
	Vxsnap->removeSnap();
	Vzsnap->removeSnap();

    result=FWI_OK;
    return result;
}

template<typename T>
int FwiElastic2D<T>::run_optimal(){
     int result = FWI_ERR;
     if(!vpgradset && !vsgradset) {
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
     vpgrad->allocateImage();
     vsgrad->allocateImage();

     // Get models for scaling
     T *Vp, *Vs, *Rho;
     Vp = model->getVp();
     Vs = model->getVs();
     Rho = model->getR();

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
                waves_fw->forwardstepVelocity(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Recording data (Vx)
                if(this->datamodVxset && !reverse){
                    waves_fw->recordData(this->datamodVx, GMAP, it);
                }

                // Recording data (Vy)
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
            if(this->datamodVzset){
                waves_fw->recordData(this->datamodVz, GMAP, capo);
            }

            // Compute Residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertSource(model, dataresVx, GMAP, capo);
            waves_bw->insertSource(model, dataresVz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getVx();
            T *wsz = waves_fw->getVz();
            T *wrx = waves_bw->getVx();
            T *wrz = waves_bw->getVz();

            crossCorr(wsx, wsz, waves_fw->getLpml(), wrx, wrz, waves_bw->getLpml(), Vp, Vs, Rho, capo);

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
            T *wrx = waves_bw->getVx();
            T *wrz = waves_bw->getVz();
            crossCorr(wsx, wsz, waves_fw->getLpml(), wrx, wrz, waves_bw->getLpml(), Vp, Vs, Rho, capo);

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
}

template<typename T>
void FwiElastic3D<T>::crossCorr(T *wsx, T*wsy, T *wsz, int pads, T* wrx, T* wry, T* wrz, int padr, T* Vp, T* Vs, T* Rho)
{
	int ix, iy, iz;
	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T msxx, msyy, mszz, msyz, msxz, msxy, mrxx, mryy, mrzz, mryz, mrxz, mrxy;
	T C33_minus, C33_plus;
	T C44_minus, C44_plus;
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
	// Getting sizes
	if(vpgradset) {
		nx = vpgrad->getNx();
		ny = vpgrad->getNy();
		nz = vpgrad->getNz();
		dx = vpgrad->getDx(); 
		dy = vpgrad->getDy(); 
		dz = vpgrad->getDz(); 
	}else{
		nx = vsgrad->getNx();
        ny = vsgrad->getNy();
        nz = vsgrad->getNz();
        dx = vsgrad->getDx(); 
        dy = vsgrad->getDy(); 
        dz = vsgrad->getDz(); 
    }

    int nxs = nx + 2*pads;
    int nxr = nx + 2*padr;
    int nys = ny + 2*pads;
    int nyr = ny + 2*padr;

    for (ix=0; ix<nx; ix++){
        for (iy=0; iy<ny; iy++){
            for (iz=0; iz<nz; iz++){
                C33_minus = Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                C33_plus = Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                C44_minus = Rho[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];
                C44_plus = Rho[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];

                msxx = (wsx[ks3D(ix+pads, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads)])/dx;
                msyy = (wsy[ks3D(ix+pads, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads)])/dy;
                mszz = (wsz[ks3D(ix+pads, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads-1)])/dz;
                mrxx = (wrx[kr3D(ix+padr, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mryy = (wry[kr3D(ix+padr, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrzz = (wrz[kr3D(ix+padr, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;

                if(vpgradset){
                    vpgraddata[ki3D(ix,iy,iz)] += C33_minus*C33_plus*(msxx + msyy + mszz) * (mrxx + mryy + mrzz);
                }

                if(vsgradset){
                    msyz = 0.5*(wsz[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msyz += 0.5*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dy;
                    msyz += 0.5*(wsy[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msyz += 0.5*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dz;

                    mryz = 0.5*(wrz[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mryz += 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dy;
                    mryz += 0.5*(wry[kr3D(ix+padr, iy+padr, iz+padr+1)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mryz += 0.5*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dz;

                    msxz = 0.5*(wsx[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msxz += 0.5*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dz;
                    msxz += 0.5*(wsz[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxz += 0.5*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dx;

                    mrxz = 0.5*(wrx[kr3D(ix+padr, iy+padr, iz+padr+1)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mrxz += 0.5*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dz;
                    mrxz += 0.5*(wrz[kr3D(ix+padr+1, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxz += 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dx;

                    msxy = 0.5*(wsx[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msxy += 0.5*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dy;
                    msxy += 0.5*(wsy[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxy += 0.5*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dx;

                    mrxy = 0.5*(wrx[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mrxy += 0.5*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dy;
                    mrxy += 0.5*(wry[kr3D(ix+padr+1, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxy += 0.5*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dx;

                    vsgraddata[ki3D(ix,iy,iz)] += C44_minus*C44_plus*(-2.0*(msyy*mrzz + mszz*mryy) -2.0*(msxx*mrzz + mszz*mrxx) -2.0*(msyy*mrxx + msxx*mryy) + msyz*mryz + msxz*mrxz + msxy*mrxy);
                }
            }
        }	
    }
}

template<typename T>
int FwiElastic3D<T>::run(){
     int result = FWI_ERR;
     if(!vpgradset && !vsgradset) {
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
    	waves->forwardstepVelocity(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

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

    // Get models for scaling
    T *Vp, *Vs, *Rho;
    Vp = model->getVp();
    Vs = model->getVs();
    Rho = model->getR();

     this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping 
    	waves->forwardstepStress(model, der);
    	waves->forwardstepVelocity(model, der);

    	// Inserting residuals
    	waves->insertSource(model, dataVx, GMAP, (nt - 1 - it));
    	waves->insertSource(model, dataVy, GMAP, (nt - 1 - it));
    	waves->insertSource(model, dataVz, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Vxsnap->readSnap(nt - 1 - it);
        Vysnap->readSnap(nt - 1 - it);
        Vzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        if((((nt - 1 - it)-Vxsnap->getEnddiff()) % Vxsnap->getSnapinc()) == 0){
            T *Vxr = waves->getVx();
            T *Vyr = waves->getVy();
            T *Vzr = waves->getVz();
            crossCorr(Vxsnap->getData(0), Vysnap->getData(0), Vzsnap->getData(0), 0, Vxr, Vyr, Vzr, waves->getLpml(), Vp, Vs, Rho);
        }

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
    
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
     if(!vpgradset && !vsgradset) {
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
     vpgrad->allocateImage();
     vsgrad->allocateImage();

     // Get models for scaling
     T *Vp, *Vs, *Rho;
     Vp = model->getVp();
     Vs = model->getVs();
     Rho = model->getR();

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
                waves_fw->forwardstepVelocity(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

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

            // Inserting data
            waves_bw->insertSource(model, dataVx, GMAP, capo);
            waves_bw->insertSource(model, dataVz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getVx();
            T *wsy = waves_fw->getVy();
            T *wsz = waves_fw->getVz();
            T *wrx = waves_bw->getVx();
            T *wry = waves_bw->getVy();
            T *wrz = waves_bw->getVz();

            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), wrx, wry, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

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

            // Inserting data
            waves_bw->insertSource(model, dataVx, GMAP, capo);
            waves_bw->insertSource(model, dataVy, GMAP, capo);
            waves_bw->insertSource(model, dataVz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getVx();
            T *wsy = waves_fw->getVy();
            T *wsz = waves_fw->getVz();
            T *wrx = waves_bw->getVx();
            T *wry = waves_bw->getVy();
            T *wrz = waves_bw->getVz();
            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), wrx, wry, wrz, waves_bw->getLpml(), Vp, Vs, Rho);

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
