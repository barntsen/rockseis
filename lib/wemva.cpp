#include "wemva.h"


namespace rockseis {

// =============== ABSTRACT WEMVA CLASS =============== //
template<typename T>
Wemva<T>::Wemva() {
    //Set default parameters
    fs = false;
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    snapmethod = OPTIMAL; 
    paramtype = PAR_GRID;
    misfit_type = SI;
    dtx = -1;
    dty = -1;
    dtz = -1;
    nhx = 1;
    nhy = 1;
    nhz = 1;
    fnorm = 0.0;
    if(createLog() == WVA_ERR)
    {
        rs_error("Wemva<T>::Wemva(): Error creating logfile for writting.");
    }
    if(createProglog() == WVA_ERR)
    {
        rs_error("Wemva<T>::Wemva(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
Wemva<T>::Wemva(MPImodeling *_mpi) {
    mpi = _mpi;

    //Set default parameters
    fs = false;
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    snapmethod = OPTIMAL; 
    paramtype = PAR_GRID;
    misfit_type = SI;
    dtx = -1;
    dty = -1;
    dtz = -1;
    nhx = 1;
    nhy = 1;
    nhz = 1;
    fnorm = 0.0;
    if(createLog() == WVA_ERR)
    {
        rs_error("Wemva<T>::Wemva(): Error creating logfile for writting.");
    }
    if(createProglog() == WVA_ERR)
    {
        rs_error("Wemva<T>::Wemva(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
void Wemva<T>::normalize(double *v, double *f, int n){
	int i;
    if(fnorm == 0.0) {
        fnorm = 1.0;
    }
    *f /= fabs(fnorm);
	for(i=0; i<n; i++) {
		v[i] /= fabs(fnorm);
	}
}

template<typename T>
double Wemva<T>::vector_norm(double *v, const int type, const int n){
	// Variables
	int i;
	double norm;

	norm = 0.0;
	if(type == 2) {
		for(i=0; i<n; i++) {
			norm += v[i]*v[i];
		}
		norm = sqrt(norm);
	}
	else if(type == 99999) {
		for(i=0; i<n; i++) {
			if(fabs(v[i]) >= norm) {
				norm = fabs(v[i]);
			}
		}
	}

	return norm;
}

template<typename T>
bool Wemva<T>::createLog(){
	logfile = LOGFILE;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return WVA_ERR;
	}else{
		Flog.close();
		return WVA_OK;
	}
}

template<typename T>
bool Wemva<T>::createProglog(){
	progresslogfile = PROGLOGFILE;
	Flog.open(progresslogfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return WVA_ERR;
	}else{
		Flog.close();
		return WVA_OK;
	}
}

template<typename T>
void Wemva<T>::writeLog(std::string text){
    char tempo[256];
	time_t now;
    now = time(NULL);
    strcpy(tempo, ctime(&now));
    tempo[strlen(tempo)-1] = '\0';
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << tempo <<": " << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Wemva<T>::writeProgress(std::string text){
    if(!progresslogfile.empty()){
        Flog.open(progresslogfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Wemva<T>::createResult(){
    struct stat s;
    // Checking if result folder is present, and creates it if not
    if(stat(RESULTDIR,&s) != 0) {
        int mkdir_return = mkdir(RESULTDIR,0777);
        if(mkdir_return != 0) rs_error("Not able to create result directory: ", RESULTDIR);
    }
}

template<typename T>
Wemva<T>::~Wemva() {
    //Do nothing
}

// =============== 2D ACOUSTIC WEMVA CLASS =============== //
//
template<typename T>
WemvaAcoustic2D<T>::WemvaAcoustic2D() {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;

    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
WemvaAcoustic2D<T>::WemvaAcoustic2D(MPImodeling *mpi): Wemva<T>(mpi) {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;
    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
WemvaAcoustic2D<T>::~WemvaAcoustic2D() {
    //Do nothing
}

template<typename T>
void WemvaAcoustic2D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> shot2D;
    std::shared_ptr<rockseis::Data2D<T>> shot2Di;
    std::shared_ptr<rockseis::Image2D<T>> pimage;
    std::shared_ptr<rockseis::Image2D<T>> lpimage;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> gmodel (new rockseis::ModelAcoustic2D<T>(Vpfile, Rhofile, this->getLpml() ,this->getFs()));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data2D<T>> source (new rockseis::Data2D<T>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<T>> interp (new rockseis::Interp<T>(SINC));

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());

	if(mpi->getRank() == 0) {
		// Master

        // Get shot map
        Sort->readKeymap();
        Sort->readSortmap();
        size_t ngathers =  Sort->getNensemb();

        // Create a data class for the recorded data
        std::shared_ptr<rockseis::Data2D<T>> shot2D (new rockseis::Data2D<T>(Precordfile));
       
        //Migration
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

        // Image
        pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile, gmodel, this->getNhx(), this->getNhz());
        pimage->createEmpty();
		for(long int i=0; i<ngathers; i++) {
            pimage->stackImage(Pimagefile + "-" + std::to_string(i));
            remove_file(Pimagefile + "-" + std::to_string(i));
        }

        //Calculate and output misfit
        this->computeMisfit(pimage);

		//Clear work vector 
		mpi->clearWork();

        // Gradient computation
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

        vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile, gmodel, 1, 1);
        vpgrad->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            remove_file(Vpgradfile + "-" + std::to_string(i));
        }

		//Clear work vector 
		mpi->clearWork();

    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::RtmAcoustic2D<T>> rtm;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do migration
                Sort->readKeymap();
                Sort->readSortmap();

                // Get the shot
                Sort->setDatafile(Precordfile);
                shot2D = Sort->get2DGather(work.id);
                size_t ntr = shot2D->getNtrace();

                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot2D, shot2Di);
                shot2Di->makeMap(lmodel->getGeom(), GMAP);

                // Creating image object
                pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile + "-" + std::to_string(work.id), lmodel, this->getNhx(), this->getNhz());

                // Create rtm object
                rtm = std::make_shared<rockseis::RtmAcoustic2D<T>>(lmodel, pimage, source, shot2Di, this->getOrder(), this->getSnapinc());

                // Setting Snapshot file 
                rtm->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting MVA flag
                rtm ->setRunmva(true);

                // Setting Snapshot parameters
                rtm->setNcheck(this->getNsnaps());
                rtm->setIncore(this->getIncore());

                // Set logfile
                rtm->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        rtm->run();
                        break;
                    case rockseis::OPTIMAL:
                        rtm->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output image
                pimage->write();

                
                // Reset all classes
                shot2D.reset();
                shot2Di.reset();
                lmodel.reset();
                pimage.reset();
                rtm.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }

        std::shared_ptr<rockseis::MvaAcoustic2D<T>> mva;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Compute gradient
                // Get the shot
                Sort->setDatafile(Precordfile);
                shot2D = Sort->get2DGather(work.id);
                size_t ntr = shot2D->getNtrace();

                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot2D, shot2Di);
                shot2Di->makeMap(lmodel->getGeom(), GMAP);

                // Create image object
                pimage = std::make_shared<rockseis::Image2D<T>>(PIMAGERESFILE);
                lpimage = pimage->getLocal(shot2D, apertx, SMAP);

                // Create mva object
                mva = std::make_shared<rockseis::MvaAcoustic2D<T>>(lmodel, lpimage, source, shot2Di, this->getOrder(), this->getSnapinc());

                // Setting misfit type
                mva->setMisfit_type(this->getMisfit_type());

                // Creating gradient object
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                // Setting up gradient objects in wemvafwi class
                mva->setVpgrad(vpgrad);

                // Setting Snapshot files
                mva->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                mva->setNcheck(this->getNsnaps());
                mva->setIncore(this->getIncore());

                // Set logfile
                mva->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        mva->run();
                        break;
                    case rockseis::OPTIMAL:
                        mva->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output gradient
                vpgrad->write();

                // Reset all classes
                shot2D.reset();
                shot2Di.reset();
                lmodel.reset();
                pimage.reset();
                lpimage.reset();
                vpgrad.reset();
                mva.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);
            }

        }
    }
}

template<typename T>
void WemvaAcoustic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;

    if(Mutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad (new rockseis::ModelAcoustic2D<T>(vpgradfile, vpgradfile, this->getLpml() ,this->getFs()));

	// Read model
	grad->readModel();
	
	T *vpgrad;
	vpgrad = grad->getVp();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl2D<T>> spline (new rockseis::Bspl2D<T>(grad->getNx(), grad->getNz(), grad->getDx(), grad->getDz(), this->getDtx(), this->getDtz(), 3, 3));
    int nc = spline->getNc();

	/* Allocating projection arrays */
	float *vpproj= (float *) calloc(nc, sizeof(float));
	if(vpproj==NULL){
		rs_error("WemvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
	}

    if(mpi->getRank() == 0) {
		// Master

        mpi->setVerbose(false); // Turn off queue printing
		// Create work queue
		for(long int i=0; i<nc; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

		//Clear work vector 
		mpi->clearWork();

		global_stack= (float *) calloc(nc, sizeof(float));
		if(global_stack==NULL){
			rs_error("WemvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
		MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   

		/* Output spline */
        std::shared_ptr<File> Fout (new File());
        Fout->output(VPPROJGRADFILE);
        Fout->setN(1,nc);
        Fout->setD(1,1.0);
        Fout->setData_format(sizeof(float));
        Fout->write(global_stack, nc, 0);
        Fout->close();

        mpi->setVerbose(true); // Turn on queue printing

       }else {
        /* Slave */
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do work
                c = spline->getSpline();
				c[work.id]=1.0; // Projection point
				spline->bisp(); // Evaluate spline for this coefficient
                wrk = spline->getMod();
				vpsum = 0.0;
				for(long int i=0; i<grad->getNx()*grad->getNz(); i++){
						vpsum += wrk[i]*vpgrad[i];
				}
				vpproj[work.id]=vpsum;
				c[work.id]=0.0; // Reset coefficient to 0
			}

			// Send result back
			work.status = WORK_FINISHED;
			mpi->sendResult(work);		
		}

		global_stack= (float *) calloc(nc, sizeof(float));
		if(global_stack==NULL){
			rs_error("WemvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
        MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

	   }
    // Free variables
    free(global_stack);
    free(vpproj);
}

template<typename T>
int WemvaAcoustic2D<T>::setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile)
{
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model_in (new rockseis::ModelAcoustic2D<T>(vpfile, rhofile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> source_in (new rockseis::Data2D<T>(sourcefile));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    model_in->readModel();  
    // Write initial model files
    model_in->setVpfile(VP0FILE);
    model_in->setRfile(RHO0FILE);
    model_in->writeModel();
    // Write linesearch model files
    model_in->setVpfile(VPLSFILE);
    model_in->setRfile(RHOLSFILE);
    model_in->writeModel();
    source_in->read();
    // Write source initial model
    source_in->setFile(SOURCE0FILE);
    source_in->write();
    // Write source linesearch file
    source_in->setFile(SOURCELSFILE);
    source_in->write();
    int N=0;
    switch(this->getParamtype()){
        case PAR_GRID:
            N = (model_in->getGeom())->getNtot();
            break;
        case PAR_BSPLINE:
            spline = std::make_shared<rockseis::Bspl2D<T>>(model_in->getNx(), model_in->getNz(), model_in->getDx(), model_in->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N=spline->getNc();
            break;
        default:
            rs_error("WemvaAcoustic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return N;
}

template<typename T>
void WemvaAcoustic2D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcels (new rockseis::Data2D<T>(SOURCELSFILE));
    lsmodel->readModel();
    name = dir + "/" + VP_UP + "-" + std::to_string(iter);
    lsmodel->setVpfile(name);
    lsmodel->writeVp();

    // Write out image
    std::shared_ptr<rockseis::Image2D<T>> pimage;
    pimage = std::make_shared<rockseis::Image2D<T>>(PIMAGEFILE);
    if(!pimage->getAllocated()) pimage->allocateImage();
    pimage->read();
    name = dir + "/" + PIMAGE_UP + "-" + std::to_string(iter);
    pimage->setImagefile(name);
    pimage->write();
}

template<typename T>
void WemvaAcoustic2D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model0 (new rockseis::ModelAcoustic2D<T>(VP0FILE, RHO0FILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> source0 (new rockseis::Data2D<T>(SOURCE0FILE));
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> lssource (new rockseis::Data2D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute;

    // Write linesearch model
    model0->readModel();
    lsmodel->readModel();
    source0->read();
    lssource->read();
    T *vp0, *rho0, *wav0, *vpls, *rhols, *wavls;
    T *c, *mod;
    T *vpmutedata;
    vp0 = model0->getVp(); 
    rho0 = model0->getR(); 
    wav0 = source0->getData();
    vpls = lsmodel->getVp(); 
    rhols = lsmodel->getR(); 
    wavls = lssource->getData();
    int i;
    int N, Ns, Nmod;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Mutefile, Mutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("WemvaAcoustic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readModel();
        vpmutedata = mute->getVp();
        N = (lsmodel->getGeom())->getNtot();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
        }
    }

    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
                vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
                rhols[i] = rho0[i];
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = wav0[i];
            }
            lssource->write();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            spline = std::make_shared<rockseis::Bspl2D<T>>(model0->getNx(), model0->getNz(), model0->getDx(), model0->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            c = spline->getSpline();
            for(i=0; i< N; i++)
            {
                c[i] = x[i];
            }
            spline->bisp();
            mod = spline->getMod();

            for(i=0; i< Nmod; i++)
            {
                vpls[i] = vp0[i] + mod[i]*vpmutedata[i]*kvp;
                rhols[i] = rho0[i];
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = wav0[i];
            }
            lssource->write();
            break;
        default:
            rs_error("WemvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Mutefile.empty()){
        free(vpmutedata);
    }
}

template<typename T>
void WemvaAcoustic2D<T>::computeMisfit(std::shared_ptr<rockseis::Image2D<T>> pimage)
{
    T f=0;
    if(!pimage->getAllocated()) pimage->allocateImage();
    int ix, iz, ihx, ihz;
    // Read image data
    pimage->read();
    T *imagedata = pimage->getImagedata();
    int nhx = pimage->getNhx();
    int nhz = pimage->getNhz();
    int nx = pimage->getNx();
    int nz = pimage->getNz();
    T *wrk = (T *) calloc(nz, sizeof(T));
    int hx, hz;
    T G1 = 0.;
    T G2 = 0.;
    T f1 = 0.;
    T f2 = 0.;

    switch(this->getMisfit_type()){
        case SI:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        // Misfit
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f -= 0.5*G2*wrk[iz]*wrk[iz];
                        }
                        // Residual
                        for (iz=1; iz<nz-1; iz++){
						    wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = G2*wrk[iz];
                        }
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        case DS:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = (hx*hx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1 + (hz*hz);
                    for (ix=0; ix<nx; ix++){
                        // Misfit
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f += 0.5*G2*wrk[iz]*wrk[iz];
                        }
                        // Residual
                        for (iz=1; iz<nz-1; iz++){
						    wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = -1.0*G2*wrk[iz];
                        }
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        case DS_PLUS_SI:
            // Misfit
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f1 += 0.5*((hx*hx) +  (hz*hz))*wrk[iz]*wrk[iz];
                            f2 += G2*wrk[iz]*wrk[iz];
                        }
                    }
                }
            }
            if(f2 != 0) {
                f = f1/f2;
            }else {
                f = f1;
            }
            //Residual
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        for (iz=1; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }	
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = (f1/(f2*f2))*G2*wrk[iz] - (((hx*hx) +  (hz*hz))*1.0/f2)*wrk[iz]; 
                        }	
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        default:
            f = 0;
            break;
    }

    // Free work array
    free(wrk);

    pimage->setImagefile(PIMAGERESFILE);
    pimage->write();

    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->output(MISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&f, 1, 0);
}

template<typename T>
void WemvaAcoustic2D<T>::readMisfit(double *f)
{
    // Data misfit
    *f = 0.0;
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->input(MISFITFILE);
    T val;
    for(int i=0; i<Fmisfit->getN(1); i++){
        Fmisfit->read(&val, 1); 
        *f += val;
    }
    Fmisfit->close();

    // Regularisation misfit
    Fmisfit->input(VPREGMISFITFILE);
    Fmisfit->read(&val, 1, 0); 
    *f += reg_alpha[0]*val;
    Fmisfit->close();
}

template<typename T>
void WemvaAcoustic2D<T>::readGrad(double *g)
{
    int i;
    int N;
    float *g_in;
    T *gvp;
    std::string vpgradfile;
    if(Mutefile.empty()){
        vpgradfile = VPGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelAcoustic2D<T>> modelgrad (new rockseis::ModelAcoustic2D<T>(vpgradfile, vpgradfile, 1 ,0));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readModel();
            N = (modelgrad->getGeom())->getNtot();
            gvp = modelgrad->getVp(); 
            for(i=0; i< N; i++)
            {
                g[i] = gvp[i]*kvp;
            }
            break;
        case PAR_BSPLINE:
           spline = std::make_shared<rockseis::Bspl2D<T>>(modelgrad->getNx(), modelgrad->getNz(), modelgrad->getDx(), modelgrad->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            g_in = (float *) calloc(2*N, sizeof(float));
            Fgrad = std::make_shared<rockseis::File>();
            Fgrad->input(VPPROJGRADFILE);
            Fgrad->read(&g_in[0], N, 0);
            Fgrad->close();
            for(i=0; i< N; i++)
            {
                g[i] = g_in[i]*kvp;
            }
            // Free temporary array
            free(g_in);
            break;
        default:
            rs_error("WemvaAcoustic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void WemvaAcoustic2D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADFILE, VPGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPREGGRADFILE, VPREGGRADFILE, 1 ,0);

    // Read gradients
    grad->readModel();
    reggrad->readModel();
    T *vp, *regvp;
    vp = grad->getVp(); 
    regvp = reggrad->getVp(); 
    int i;
    int N;

    N = (grad->getGeom())->getNtot();
    // Compute 
    for(i=0; i< N; i++)
    {
        vp[i] = vp[i] + reg_alpha[0]*regvp[i];
    }
    grad->setVpfile(VPGRADCOMBFILE);
    grad->setRfile(VPGRADCOMBFILE);
    grad->writeVp();
}

template<typename T>
void WemvaAcoustic2D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADCOMBFILE, VPGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute (new rockseis::ModelAcoustic2D<T>(Mutefile, Mutefile, 1 ,0));

        // Mute model and write
        model->readModel();
        mute->readModel();
        T *vp, *vpmute;
        vp = model->getVp(); 
        vpmute = mute->getVp(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
        }
        model->setVpfile(VPGRADMUTEFILE);
        model->setRfile(VPGRADMUTEFILE);
        model->writeVp();
    }
}

template<typename T>
void WemvaAcoustic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), 1, model->getNz(), model->getDx(), 1.0, model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx, *dvpdz;
    T *vpgrad;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    gwrk = (double *) calloc(Nmod, sizeof(double));
    model->readModel();
    model->setVpfile(VPREGGRADFILE);
    model->setRfile(VPREGGRADFILE);
    vpgrad = model->getVp();
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (model->getGeom())->getNtot();
            der->ddx_fw(x);
            for(i=0; i< N; i++)
            {
                dvpdx[i] = df[i]*kvp;
            }
            der->ddz_fw(x);
            for(i=0; i< N; i++)
            {
                dvpdz[i] = df[i]*kvp;
            }
            
            break;
        case PAR_BSPLINE:
            Nmod = (model->getGeom())->getNtot();
            spline = std::make_shared<rockseis::Bspl2D<double>>(model->getNx(), model->getNz(), model->getDx(), model->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            c = spline->getSpline();
            for(i=0; i< N; i++)
            {
                c[i] = x[i];
            }
            spline->bisp();
            mod = spline->getMod();
            der->ddx_fw(mod);
            for(i=0; i< Nmod; i++)
            {
                dvpdx[i] = df[i]*kvp;
            }
            der->ddz_fw(mod);
            for(i=0; i< Nmod; i++)
            {
                dvpdz[i] = df[i]*kvp;
            }

            break;
        default:
            rs_error("WemvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    // Computing misfit
    double M; 
    T fvp = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i];
        M = sqrt(M);
        fvp += M;
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdx[i]/M;
    }
    der->ddx_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] -= df[i];
    }


    /* Write out misfit */
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->output(VPREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&fvp,1,0);
    Fmisfit->close();

    /* Write out gradient */
    model->writeVp();

    // Free variables
    free(dvpdx);
    free(dvpdz);
    free(gwrk);
}

// =============== 2D ELASTIC WEMVA CLASS =============== //
template<typename T>
WemvaElastic2D<T>::WemvaElastic2D() {
    // Set default parameters
    apertx = -1;
    sourcetype = 0;

    kvp = 1.0;
    kvs = 1.0;

    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;

    update_vp = false;
    update_vs = false;
    wavemode = 0;
}

template<typename T>
WemvaElastic2D<T>::WemvaElastic2D(MPImodeling *mpi): Wemva<T>(mpi) {
    // Set default parameters
    apertx = -1;
    sourcetype = 0;

    kvp = 1.0;
    kvs = 1.0;
    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;

    update_vp = false;
    update_vs = false;
    wavemode = 0;
}

template<typename T>
WemvaElastic2D<T>::~WemvaElastic2D() {
    //Do nothing
}

template<typename T>
void WemvaElastic2D<T>::setWavemode(int mode) {
    if(mode != 0 && mode != 1) rs_error("WemvaElastic2D<T>::setWavemod: Invalid wavemode.");
    this->wavemode = mode;
    if(this->wavemode == 1){
        update_vp = false;
        update_vs = true;
    }else{
        update_vp = true;
        update_vs = false;
    }
}

template<typename T>
void WemvaElastic2D<T>::runPPgrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2Di;
    std::shared_ptr<rockseis::Image2D<T>> pimage;
    std::shared_ptr<rockseis::Image2D<T>> lpimage;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::Image2D<T>> vsgrad;
    std::shared_ptr<rockseis::ModelElastic2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Uxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelElastic2D<T>> gmodel (new rockseis::ModelElastic2D<T>(Vpfile, Vsfile, Rhofile, this->getLpml() ,this->getFs()));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data2D<T>> source (new rockseis::Data2D<T>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<T>> interp (new rockseis::Interp<T>(SINC));

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());

	if(mpi->getRank() == 0) {
		// Master

        // Get shot map
        Sort->readKeymap();
        Sort->readSortmap();
        size_t ngathers =  Sort->getNensemb();
        if(Sort->getReciprocity()){
            ngathers *= 2;
        }

        //Migration
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

        // Image
        pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile, gmodel, this->getNhx(), this->getNhz());
        pimage->createEmpty();
		for(long int i=0; i<ngathers; i++) {
            pimage->stackImage(Pimagefile + "-" + std::to_string(i));
            remove_file(Pimagefile + "-" + std::to_string(i));
        }

        //Calculate and output misfit
        this->computeMisfit(pimage, PIMAGERESFILE);

		//Clear work vector 
		mpi->clearWork();

        // Gradient computation
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

       // Images
       vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile, gmodel, 1, 1);
       vpgrad->createEmpty();

       vsgrad = std::make_shared<rockseis::Image2D<T>>(Vsgradfile, gmodel, 1, 1);
       vsgrad->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            if(update_vp){
                vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
                remove_file(Vpgradfile + "-" + std::to_string(i));
            }
        }


		//Clear work vector 
		mpi->clearWork();

    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::RtmElastic2D<T>> rtm;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do migration
                Sort->readKeymap();
                Sort->readSortmap();

                int gatherid;
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/2;
                }

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata2D = Sort->get2DGather(gatherid);
                size_t ntr = Uxdata2D->getNtrace();

                Sort->setDatafile(Uzrecordfile);
                Uzdata2D = Sort->get2DGather(gatherid);

                // Creating local model
                lmodel = gmodel->getLocal(Uxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                if(!Sort->getReciprocity()){
                    switch(this->getSourcetype()){
                        case 0:
                            source->setField(PRESSURE);
                            break;
                        case 1:
                            source->setField(VX);
                            break;
                        case 3:
                            source->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }else{
                    if(work.id % 2 == 0){
                            source->setField(VX);
                    }else{
                            source->setField(VZ);
                    }
                }

                // Interpolate shot
                Uxdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata2D, Uxdata2Di);
                Uxdata2Di->makeMap(lmodel->getGeom(), GMAP);

                Uzdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata2D, Uzdata2Di);
                Uzdata2Di->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdata2Di->setField(VX);
                    Uzdata2Di->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdata2Di->setField(PRESSURE);
                            Uzdata2Di->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdata2Di->setField(VX);
                            Uzdata2Di->setField(VX);
                            break;
                        case 3:
                            Uxdata2Di->setField(VZ);
                            Uzdata2Di->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }


                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Creating image object
                pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile + "-" + std::to_string(work.id), lmodel, this->getNhx(), this->getNhz());

                // Create rtm object
                rtm = std::make_shared<rockseis::RtmElastic2D<T>>(lmodel, source, Uxdata2Di, Uzdata2Di, this->getOrder(), this->getSnapinc());

                // Setting Snapshot file 
                rtm->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting MVA flag
                rtm ->setRunmva(true);

                // Setting Image objects
                rtm->setPimage(pimage);

                // Setting Snapshot parameters
                rtm->setNcheck(this->getNsnaps());
                rtm->setIncore(this->getIncore());

                // Set logfile
                rtm->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        rtm->run();
                        break;
                    case rockseis::OPTIMAL:
                        rtm->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output image
                pimage->write();
                
                // Reset all classes
                Uxdata2D.reset();
                Uxdata2Di.reset();
                Uzdata2D.reset();
                Uzdata2Di.reset();
                lmodel.reset();
                pimage.reset();
                rtm.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }
        // Compute gradient 
        std::shared_ptr<rockseis::PPmvaElastic2D<T>> mva;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                int gatherid;
                // Compute gradient
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/2;
                }

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata2D = Sort->get2DGather(gatherid);
                size_t ntr = Uxdata2D->getNtrace();

                Sort->setDatafile(Uzrecordfile);
                Uzdata2D = Sort->get2DGather(gatherid);

                // Creating local model
                lmodel = gmodel->getLocal(Uxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                if(!Sort->getReciprocity()){
                    switch(this->getSourcetype()){
                        case 0:
                            source->setField(PRESSURE);
                            break;
                        case 1:
                            source->setField(VX);
                            break;
                        case 3:
                            source->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }else{
                    if(work.id % 2 == 0){
                        source->setField(VX);
                    }else{
                        source->setField(VZ);
                    }
                }

                // Interpolate shot
                Uxdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata2D, Uxdata2Di);
                Uxdata2Di->makeMap(lmodel->getGeom(), GMAP);

                Uzdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata2D, Uzdata2Di);
                Uzdata2Di->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdata2Di->setField(VX);
                    Uzdata2Di->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdata2Di->setField(PRESSURE);
                            Uzdata2Di->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdata2Di->setField(VX);
                            Uzdata2Di->setField(VX);
                            break;
                        case 3:
                            Uxdata2Di->setField(VZ);
                            Uzdata2Di->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }

                // Create image object
                pimage = std::make_shared<rockseis::Image2D<T>>(PIMAGERESFILE);
                lpimage = pimage->getLocal(Uxdata2D, apertx, SMAP);

                // Create mva object
                mva = std::make_shared<rockseis::PPmvaElastic2D<T>>(lmodel, lpimage, source, Uxdata2Di, Uzdata2Di, this->getOrder(), this->getSnapinc());

                // Setting misfit type
                mva->setMisfit_type(this->getMisfit_type());

                // Creating gradient object
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                // Setting up gradient objects in wemvafwi class
                mva->setVpgrad(vpgrad);

                // Setting Snapshot files
                mva->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                mva->setNcheck(this->getNsnaps());
                mva->setIncore(this->getIncore());

                // Set logfile
                mva->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        mva->run();
                        break;
                    case rockseis::OPTIMAL:
                        mva->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output gradient
                vpgrad->write();

                // Reset all classes
                Uxdata2D.reset();
                Uxdata2Di.reset();
                Uzdata2D.reset();
                Uzdata2Di.reset();
                lmodel.reset();
                pimage.reset();
                lpimage.reset();
                vpgrad.reset();
                mva.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);
            }
        }
    }
}

template<typename T>
void WemvaElastic2D<T>::runPSgrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2Di;
    std::shared_ptr<rockseis::Image2D<T>> simage;
    std::shared_ptr<rockseis::Image2D<T>> lsimage;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::Image2D<T>> vsgrad;
    std::shared_ptr<rockseis::ModelElastic2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Uxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelElastic2D<T>> gmodel (new rockseis::ModelElastic2D<T>(Vpfile, Vsfile, Rhofile, this->getLpml() ,this->getFs()));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data2D<T>> source (new rockseis::Data2D<T>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<T>> interp (new rockseis::Interp<T>(SINC));

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());

	if(mpi->getRank() == 0) {
		// Master

        // Get shot map
        Sort->readKeymap();
        Sort->readSortmap();
        size_t ngathers =  Sort->getNensemb();
        if(Sort->getReciprocity()){
            ngathers *= 2;
        }

        //Migration
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

        // Image
        simage = std::make_shared<rockseis::Image2D<T>>(Simagefile, gmodel, this->getNhx(), this->getNhz());
        simage->createEmpty();
		for(long int i=0; i<ngathers; i++) {
            simage->stackImage(Simagefile + "-" + std::to_string(i));
            remove_file(Simagefile + "-" + std::to_string(i));
        }

        //Calculate and output misfit
        this->computeMisfit(simage, SIMAGERESFILE);

		//Clear work vector 
		mpi->clearWork();

        // Gradient computation
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

       // Images
       vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile, gmodel, 1, 1);
       vpgrad->createEmpty();

       vsgrad = std::make_shared<rockseis::Image2D<T>>(Vsgradfile, gmodel, 1, 1);
       vsgrad->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            if(update_vs){
                vsgrad->stackImage(Vsgradfile + "-" + std::to_string(i));
                remove_file(Vsgradfile + "-" + std::to_string(i));
            }
        }


		//Clear work vector 
		mpi->clearWork();

    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::RtmElastic2D<T>> rtm;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do migration
                Sort->readKeymap();
                Sort->readSortmap();

                int gatherid;
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/2;
                }

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata2D = Sort->get2DGather(gatherid);
                size_t ntr = Uxdata2D->getNtrace();

                Sort->setDatafile(Uzrecordfile);
                Uzdata2D = Sort->get2DGather(gatherid);

                // Creating local model
                lmodel = gmodel->getLocal(Uxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                if(!Sort->getReciprocity()){
                    switch(this->getSourcetype()){
                        case 0:
                            source->setField(PRESSURE);
                            break;
                        case 1:
                            source->setField(VX);
                            break;
                        case 3:
                            source->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }else{
                    if(work.id % 2 == 0){
                            source->setField(VX);
                    }else{
                            source->setField(VZ);
                    }
                }

                // Interpolate shot
                Uxdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata2D, Uxdata2Di);
                Uxdata2Di->makeMap(lmodel->getGeom(), GMAP);

                Uzdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata2D, Uzdata2Di);
                Uzdata2Di->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdata2Di->setField(VX);
                    Uzdata2Di->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdata2Di->setField(PRESSURE);
                            Uzdata2Di->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdata2Di->setField(VX);
                            Uzdata2Di->setField(VX);
                            break;
                        case 3:
                            Uxdata2Di->setField(VZ);
                            Uzdata2Di->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }


                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Creating image object
                simage = std::make_shared<rockseis::Image2D<T>>(Simagefile + "-" + std::to_string(work.id), lmodel, this->getNhx(), this->getNhz());

                // Create rtm object
                rtm = std::make_shared<rockseis::RtmElastic2D<T>>(lmodel, source, Uxdata2Di, Uzdata2Di, this->getOrder(), this->getSnapinc());

                // Setting Snapshot file 
                rtm->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting MVA flag
                rtm ->setRunmva(true);

                // Setting Image objects
                rtm->setSimage(simage);

                // Setting Snapshot parameters
                rtm->setNcheck(this->getNsnaps());
                rtm->setIncore(this->getIncore());

                // Set logfile
                rtm->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        rtm->run();
                        break;
                    case rockseis::OPTIMAL:
                        rtm->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output image
                simage->write();
                
                // Reset all classes
                Uxdata2D.reset();
                Uxdata2Di.reset();
                Uzdata2D.reset();
                Uzdata2Di.reset();
                lmodel.reset();
                simage.reset();
                rtm.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }
        // Compute gradient 
        std::shared_ptr<rockseis::PSmvaElastic2D<T>> mva;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                int gatherid;
                // Compute gradient
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/2;
                }

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata2D = Sort->get2DGather(gatherid);
                size_t ntr = Uxdata2D->getNtrace();

                Sort->setDatafile(Uzrecordfile);
                Uzdata2D = Sort->get2DGather(gatherid);

                // Creating local model
                lmodel = gmodel->getLocal(Uxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                if(!Sort->getReciprocity()){
                    switch(this->getSourcetype()){
                        case 0:
                            source->setField(PRESSURE);
                            break;
                        case 1:
                            source->setField(VX);
                            break;
                        case 3:
                            source->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }else{
                    if(work.id % 2 == 0){
                        source->setField(VX);
                    }else{
                        source->setField(VZ);
                    }
                }

                // Interpolate shot
                Uxdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata2D, Uxdata2Di);
                Uxdata2Di->makeMap(lmodel->getGeom(), GMAP);

                Uzdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata2D, Uzdata2Di);
                Uzdata2Di->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdata2Di->setField(VX);
                    Uzdata2Di->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdata2Di->setField(PRESSURE);
                            Uzdata2Di->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdata2Di->setField(VX);
                            Uzdata2Di->setField(VX);
                            break;
                        case 3:
                            Uxdata2Di->setField(VZ);
                            Uzdata2Di->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }

                // Create image object
                simage = std::make_shared<rockseis::Image2D<T>>(SIMAGERESFILE);
                lsimage = simage->getLocal(Uxdata2D, apertx, SMAP);

                // Create mva object
                mva = std::make_shared<rockseis::PSmvaElastic2D<T>>(lmodel, lsimage, source, Uxdata2Di, Uzdata2Di, this->getOrder(), this->getSnapinc());

                // Setting misfit type
                mva->setMisfit_type(this->getMisfit_type());

                // Creating gradient object
                vsgrad = std::make_shared<rockseis::Image2D<T>>(Vsgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                // Setting up gradient objects in wemvafwi class
                mva->setVsgrad(vsgrad);

                // Setting Snapshot files
                mva->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                mva->setNcheck(this->getNsnaps());
                mva->setIncore(this->getIncore());

                // Set logfile
                mva->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        mva->run();
                        break;
                    case rockseis::OPTIMAL:
                        mva->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output gradient
                vsgrad->write();

                // Reset all classes
                Uxdata2D.reset();
                Uxdata2Di.reset();
                Uzdata2D.reset();
                Uzdata2Di.reset();
                lmodel.reset();
                simage.reset();
                lsimage.reset();
                vsgrad.reset();
                mva.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);
            }
        }
    }
}


template<typename T>
void WemvaElastic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	T vssum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;
    std::string vsgradfile;

    if(Mutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
        vsgradfile = VSGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelElastic2D<T>> grad (new rockseis::ModelElastic2D<T>(vpgradfile, vsgradfile, vpgradfile, this->getLpml() ,this->getFs()));
    
    //Read gradients
    grad->readModel();

	T *vpgrad, *vsgrad;
	vpgrad = grad->getVp();
	vsgrad = grad->getVs();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl2D<T>> spline (new rockseis::Bspl2D<T>(grad->getNx(), grad->getNz(), grad->getDx(), grad->getDz(), this->getDtx(), this->getDtz(), 3, 3));
    int nc = spline->getNc();

    float *vpproj, *vsproj;
    /* Allocating projection arrays */
    vpproj= (float *) calloc(nc, sizeof(float));
    if(vpproj==NULL){
        rs_error("WemvaElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
    }
    vsproj= (float *) calloc(nc, sizeof(float));
    if(vsproj==NULL){
        rs_error("WemvaElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vsproj)");
    }

    if(mpi->getRank() == 0) {
		// Master

        mpi->setVerbose(false); // Turn off queue printing
		// Create work queue
		for(long int i=0; i<nc; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

		//Clear work vector 
		mpi->clearWork();

		global_stack= (float *) calloc(nc, sizeof(float));
		if(global_stack==NULL){
			rs_error("WemvaElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

        std::shared_ptr<File> Fout (new File());
        if(update_vp){
            /* Starting reduce operation */
            MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   
            /* Output spline */
            Fout->output(VPPROJGRADFILE);
            Fout->setN(1,nc);
            Fout->setD(1,1.0);
            Fout->setData_format(sizeof(float));
            Fout->write(global_stack, nc, 0);
            Fout->close();

            for(long int i=0; i< nc; i++){
                global_stack[i] = 0.0;
            }
        }

        if(update_vs){
            /* Starting reduce operation */
            MPI_Reduce(vsproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   
            /* Output spline */
            Fout->output(VSPROJGRADFILE);
            Fout->setN(1,nc);
            Fout->setD(1,1.0);
            Fout->setData_format(sizeof(float));
            Fout->write(global_stack, nc, 0);
            Fout->close();

            for(long int i=0; i< nc; i++){
                global_stack[i] = 0.0;
            }
        }

        mpi->setVerbose(true); // Turn on queue printing

       }else {
        /* Slave */
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do work
                c = spline->getSpline();
                c[work.id]=1.0; // Projection point
                spline->bisp(); // Evaluate spline for this coefficient
                wrk = spline->getMod();
                vpsum = 0.0;
                vssum = 0.0;
                for(long int i=0; i<grad->getNx()*grad->getNz(); i++){
                    if(update_vp){
                        vpsum += wrk[i]*vpgrad[i];
                    }
                    if(update_vs){
                        vssum += wrk[i]*vsgrad[i];
                    }
                }
                vpproj[work.id]=vpsum;
                vsproj[work.id]=vssum;
                c[work.id]=0.0; // Reset coefficient to 0
            }

            // Send result back
            work.status = WORK_FINISHED;
            mpi->sendResult(work);		
        }

        global_stack= (float *) calloc(nc, sizeof(float));
        if(global_stack==NULL){
            rs_error("WemvaElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
        }

        /* Starting reduce operation */
        if(update_vp){
            MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
        if(update_vs){
            MPI_Reduce(vsproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
       }
    // Free allocated variables
    free(vpproj);
    free(vsproj);
    free(global_stack);
}
		
template<typename T>
int WemvaElastic2D<T>::setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile)
{
    std::shared_ptr<rockseis::ModelElastic2D<T>> model_in (new rockseis::ModelElastic2D<T>(vpfile, vsfile, rhofile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> source_in (new rockseis::Data2D<T>(sourcefile));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    model_in->readModel();  
    // Write initial model files
    model_in->setVpfile(VP0FILE);
    model_in->setVsfile(VS0FILE);
    model_in->setRfile(RHO0FILE);
    model_in->writeModel();
    // Write linesearch model files
    model_in->setVpfile(VPLSFILE);
    model_in->setVsfile(VSLSFILE);
    model_in->setRfile(RHOLSFILE);
    model_in->writeModel();
    // Write source initial model
    source_in->read();
    source_in->setFile(SOURCE0FILE);
    source_in->write();
    // Write source linesearch file
    source_in->setFile(SOURCELSFILE);
    source_in->write();
    int N, Ns;
    int Npar = 0;
    switch(this->getParamtype()){
        case PAR_GRID:
            N = (model_in->getGeom())->getNtot();
            Ns = source_in->getNt();
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_vs) Npar += N;
            break;
        case PAR_BSPLINE:
             spline = std::make_shared<rockseis::Bspl2D<T>>(model_in->getNx(), model_in->getNz(), model_in->getDx(), model_in->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N=spline->getNc();
            Ns = source_in->getNt();
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_vs) Npar += N;
            break;
        default:
            rs_error("WemvaElastic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void WemvaElastic2D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelElastic2D<T>> lsmodel (new rockseis::ModelElastic2D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcels (new rockseis::Data2D<T>(SOURCELSFILE));
    lsmodel->readModel();
    sourcels->read();
    if(update_vp){
        // Write out model
        name = dir + "/" + VP_UP + "-" + std::to_string(iter);
        lsmodel->setVpfile(name);
        lsmodel->writeVp();
        // Write out image
        std::shared_ptr<rockseis::Image2D<T>> pimage;
        pimage = std::make_shared<rockseis::Image2D<T>>(PIMAGEFILE);
        if(!pimage->getAllocated()) pimage->allocateImage();
        pimage->read();
        name = dir + "/" + PIMAGE_UP + "-" + std::to_string(iter);
        pimage->setImagefile(name);
        pimage->write();

    }
    if(update_vs){
        name = dir + "/" + VS_UP + "-" + std::to_string(iter);
        lsmodel->setVsfile(name);
        lsmodel->writeVs();
        // Write out image
        std::shared_ptr<rockseis::Image2D<T>> simage;
        simage = std::make_shared<rockseis::Image2D<T>>(SIMAGEFILE);
        if(!simage->getAllocated()) simage->allocateImage();
        simage->read();
        name = dir + "/" + SIMAGE_UP + "-" + std::to_string(iter);
        simage->setImagefile(name);
        simage->write();
    }
}

template<typename T>
void WemvaElastic2D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic2D<T>> model0 (new rockseis::ModelElastic2D<T>(VP0FILE, VS0FILE, RHO0FILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> source0 (new rockseis::Data2D<T>(SOURCE0FILE));
    std::shared_ptr<rockseis::ModelElastic2D<T>> lsmodel (new rockseis::ModelElastic2D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> lssource (new rockseis::Data2D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelElastic2D<T>> mute;

    // Write linesearch model
    model0->readModel();
    lsmodel->readModel();
    source0->read();
    lssource->read();
    T *vp0, *vs0, *rho0, *wav0, *vpls, *vsls, *rhols, *wavls;
    T *c, *mod;
    T *vpmutedata;
    T *vsmutedata;
    vp0 = model0->getVp(); 
    vs0 = model0->getVs(); 
    rho0 = model0->getR(); 
    wav0 = source0->getData();
    vpls = lsmodel->getVp(); 
    vsls = lsmodel->getVs(); 
    rhols = lsmodel->getR(); 
    wavls = lssource->getData();
    int i;
    int N=0, Ns=0, Nmod=0, Npar=0;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelElastic2D<T>>(Mutefile, Mutefile, Mutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("WemvaElastic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readModel();
        vpmutedata = mute->getVp();
        vsmutedata = mute->getVs();
        N = (lsmodel->getGeom())->getNtot();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        vsmutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
            vsmutedata[i] = 1.0;
        }
    }

    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            Npar = 0;
            if(update_vp){
                for(i=0; i< N; i++)
                {
                    vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    vpls[i] = vp0[i];
                }
            }
            if(update_vs){
                for(i=0; i< N; i++)
                {
                    vsls[i] = vs0[i] + x[Npar+i]*vsmutedata[i]*kvs;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    vsls[i] = vs0[i];
                }
            }
            for(i=0; i< N; i++)
            {
                rhols[i] = rho0[i];
            }

            /*Ensure vp/vs boundary */
            for(i=0; i< N; i++)
            {
                if(vpls[i] < 1.2*vsls[i])
                {
                    vsls[i] = vpls[i]/1.2;
                }
            }
            lsmodel->writeModel();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            Npar = 0;
            spline = std::make_shared<rockseis::Bspl2D<T>>(model0->getNx(), model0->getNz(), model0->getDx(), model0->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            c = spline->getSpline();
            if(update_vp){
                for(i=0; i< N; i++)
                {
                    c[i] = x[i];
                }
                spline->bisp();
                mod = spline->getMod();
                for(i=0; i< Nmod; i++)
                {
                    vpls[i] = vp0[i] + mod[i]*vpmutedata[i]*kvp;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    vpls[i] = vp0[i];
                }
            }
            if(update_vs){
                for(i=0; i< N; i++)
                {
                    c[i] = x[Npar+i];
                }
                spline->bisp();
                mod = spline->getMod();

                for(i=0; i< Nmod; i++)
                {
                    vsls[i] = vs0[i] + mod[i]*vsmutedata[i]*kvs;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    vsls[i] = vs0[i];
                }
            }
            for(i=0; i< Nmod; i++)
            {
                rhols[i] = rho0[i];
            }


            /*Ensure vp/vs boundary */
            for(i=0; i< Nmod; i++){
                if(vpls[i] < 1.2*vsls[i])
                {
                    vsls[i] = vpls[i]/1.2;
                }
            }

            lsmodel->writeModel();
            break;
        default:
            rs_error("WemvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }

    // Source wavelet
    Ns = lssource->getNt();
    for(i=0; i< Ns; i++)
    {
        wavls[i] = wav0[i];
    }
    lssource->write();

    // Free allocated arrays
    if(Mutefile.empty()){
        free(vpmutedata);
        free(vsmutedata);
    }
}

template<typename T>
void WemvaElastic2D<T>::readMisfit(double *f)
{
    // Data misfit
    *f = 0.0;
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->input(MISFITFILE);
    T val;
    for(int i=0; i<Fmisfit->getN(1); i++){
        Fmisfit->read(&val, 1); 
        *f += val;
    }
    Fmisfit->close();

    // Regularisation misfit
    Fmisfit->input(VPREGMISFITFILE);
    Fmisfit->read(&val, 1, 0); 
    *f += reg_alpha[0]*val;
    Fmisfit->close();

    Fmisfit->input(VSREGMISFITFILE);
    Fmisfit->read(&val, 1, 0); 
    *f += reg_alpha[1]*val;
    Fmisfit->close();
}

template<typename T>
void WemvaElastic2D<T>::readGrad(double *g)
{
    int i;
    int N,Npar=0;
    float *g_in;
    T *gvp, *gvs;
    std::string vpgradfile;
    std::string vsgradfile;
    if(Mutefile.empty()){
        vpgradfile = VPGRADFILE;
        vsgradfile = VSGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelElastic2D<T>> modelgrad (new rockseis::ModelElastic2D<T>(vpgradfile, vsgradfile, vpgradfile, 1 ,0));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readModel();
            N = (modelgrad->getGeom())->getNtot();
            Npar = 0;
            gvp = modelgrad->getVp(); 
            gvs = modelgrad->getVs(); 
            if(update_vp){
                for(i=0; i < N; i++)
                {
                    g[i] = gvp[i]*kvp;
                }
                Npar += N;
            }
            if(update_vs){
                for(i=0; i < N; i++)
                {
                    g[Npar+i] = gvs[i]*kvs;
                }
                Npar += N;
            }
            break;
        case PAR_BSPLINE:
           spline = std::make_shared<rockseis::Bspl2D<T>>(modelgrad->getNx(), modelgrad->getNz(), modelgrad->getDx(), modelgrad->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            Npar = 0;
            g_in = (float *) calloc(N, sizeof(float));
            Fgrad = std::make_shared<rockseis::File>();

            if(update_vp){
                Fgrad->input(VPPROJGRADFILE);
                Fgrad->read(&g_in[0], N, 0);
                Fgrad->close();
                for(i=0; i< N; i++)
                {
                    g[i] = g_in[i]*kvp;
                }
                Npar += N;
            }
            if(update_vs){
                Fgrad->input(VSPROJGRADFILE);
                Fgrad->read(&g_in[0], N, 0);
                Fgrad->close();
                for(i=0; i< N; i++)
                {
                    g[Npar+i] = g_in[i]*kvs;
                }
                Npar += N;
            }
            // Free temporary array
            free(g_in);
            break;
        default:
            rs_error("WemvaElastic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void WemvaElastic2D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelElastic2D<T>> grad;
    std::shared_ptr<rockseis::ModelElastic2D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADFILE, VSGRADFILE, VPGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelElastic2D<T>>(VPREGGRADFILE, VSREGGRADFILE, VPREGGRADFILE, 1 ,0);

    // Read gradients
    grad->readModel();
    reggrad->readModel();
    T *vp, *vs, *regvp, *regvs;
    vp = grad->getVp(); 
    vs = grad->getVs(); 
    regvp = reggrad->getVp(); 
    regvs = reggrad->getVs(); 
    int i;
    int N;

    N = (grad->getGeom())->getNtot();
    // Compute 
    for(i=0; i< N; i++)
    {
        if(update_vp){
            vp[i] = vp[i] + reg_alpha[0]*regvp[i];
        }
        if(update_vs){
            vs[i] = vs[i] + reg_alpha[1]*regvs[i];
        }
    }
    grad->setVpfile(VPGRADCOMBFILE);
    grad->setVsfile(VSGRADCOMBFILE);
    grad->writeVp();
    grad->writeVs();
}


template<typename T>
void WemvaElastic2D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelElastic2D<T>> model;
        model = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADCOMBFILE, VSGRADCOMBFILE, VPGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelElastic2D<T>> mute (new rockseis::ModelElastic2D<T>(Mutefile, Mutefile, Mutefile, 1 ,0));

        // Mute model and write
        model->readModel();
        mute->readModel();
        T *vp, *vs, *vpmute, *vsmute;
        vp = model->getVp(); 
        vs = model->getVs(); 
        vpmute = mute->getVp(); 
        vsmute = mute->getVs(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
            vs[i] = vs[i]*vsmute[i];
        }
        model->setVpfile(VPGRADMUTEFILE);
        model->setVsfile(VSGRADMUTEFILE);
        model->writeVp();
        model->writeVs();
    }
}

template<typename T>
void WemvaElastic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic2D<T>> model (new rockseis::ModelElastic2D<T>(VPLSFILE, VSLSFILE, VPLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), 1, model->getNz(), model->getDx(), 1.0, model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx, *dvpdz;
    double *dvsdx, *dvsdz;
    T *vpgrad, *vsgrad;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Npar=0,Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    dvsdx = (double *) calloc(Nmod, sizeof(double));
    dvsdz = (double *) calloc(Nmod, sizeof(double));
    gwrk = (double *) calloc(Nmod, sizeof(double));
    model->readModel();
    model->setVpfile(VPREGGRADFILE);
    model->setVsfile(VSREGGRADFILE);
    model->setRfile(VPREGGRADFILE);
    vpgrad = model->getVp();
    vsgrad = model->getVs();
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (model->getGeom())->getNtot();
            Npar=0;
            if(update_vp){
                der->ddx_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdx[i] = df[i]*kvp;
                }
                der->ddz_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdz[i] = df[i]*kvp;
                }
                Npar +=N;
            }
            if(update_vs){
                der->ddx_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    dvsdx[i] = df[i]*kvs;
                }
                der->ddz_fw(x);
                for(i=0; i< N; i++)
                {
                    dvsdz[i] = df[i]*kvs;
                }
                Npar +=N;
            }
            break;
        case PAR_BSPLINE:
            Nmod = (model->getGeom())->getNtot();
            Npar=0;
            spline = std::make_shared<rockseis::Bspl2D<double>>(model->getNx(), model->getNz(), model->getDx(), model->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N = spline->getNc();
            c = spline->getSpline();
            if(update_vp){
                for(i=0; i< N; i++)
                {
                    c[i] = x[i];
                }
                spline->bisp();
                mod = spline->getMod();
                der->ddx_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvpdx[i] = df[i]*kvp;
                }
                der->ddz_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvpdz[i] = df[i]*kvp;
                }
                Npar += N;
            }
            if(update_vs){
                for(i=0; i< N; i++)
                {
                    c[i] = x[Npar+i];
                }
                spline->bisp();
                mod = spline->getMod();
                der->ddx_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvsdx[i] = df[i]*kvs;
                }
                der->ddz_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvsdz[i] = df[i]*kvs;
                }
                Npar += N;
            }

            break;
        default:
            rs_error("WemvaElastic2D<T>::computeRegularization(): Unknown parameterisation."); 
            break;
    }
    // Computing misfit
    double M; 
    T fvp = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i];
        M = sqrt(M);
        fvp += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdx[i]/M;
    }
    der->ddx_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] -= df[i];
    }

    // Computing misfit
    T fvs = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = dvsdx[i]*dvsdx[i] + dvsdz[i]*dvsdz[i];
        M = sqrt(M);
        fvs += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = dvsdx[i]*dvsdx[i] + dvsdz[i]*dvsdz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = dvsdx[i]/M;
    }
    der->ddx_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vsgrad[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvsdx[i]*dvsdx[i] + dvsdz[i]*dvsdz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = dvsdz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vsgrad[i] -= df[i];
    }

    /* Write out misfit */
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->output(VPREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&fvp,1,0);
    Fmisfit->close();

    Fmisfit->output(VSREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&fvs,1,0);
    Fmisfit->close();

    /* Write out gradient */
    model->writeVp();
    model->writeVs();

    // Free variables
    free(dvpdx);
    free(dvpdz);
    free(dvsdx);
    free(dvsdz);
    free(gwrk);
}

template<typename T>
void WemvaElastic2D<T>::computeMisfit(std::shared_ptr<rockseis::Image2D<T>> image, std::string imageresfile)
{
    T f=0;
    if(!image->getAllocated()) image->allocateImage();
    int ix, iz, ihx, ihz;
    // Read image data
    image->read();
    T *imagedata = image->getImagedata();
    int nhx = image->getNhx();
    int nhz = image->getNhz();
    int nx = image->getNx();
    int nz = image->getNz();
    T *wrk = (T *) calloc(nz, sizeof(T));
    int hx, hz;
    T G1 = 0.;
    T G2 = 0.;
    T f1 = 0.;
    T f2 = 0.;

    switch(this->getMisfit_type()){
        case SI:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        // Misfit
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f -= 0.5*G2*wrk[iz]*wrk[iz];
                        }
                        // Residual
                        for (iz=1; iz<nz-1; iz++){
						    wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = G2*wrk[iz];
                        }
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        case DS:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = (hx*hx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1 + (hz*hz);
                    for (ix=0; ix<nx; ix++){
                        // Misfit
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f += 0.5*G2*wrk[iz]*wrk[iz];
                        }
                        // Residual
                        for (iz=1; iz<nz-1; iz++){
						    wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = -1.0*G2*wrk[iz];
                        }
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        case DS_PLUS_SI:
            // Misfit
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        for (iz=0; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                        for (iz=0; iz<nz; iz++){
                            f1 += 0.5*((hx*hx) +  (hz*hz))*wrk[iz]*wrk[iz];
                            f2 += G2*wrk[iz]*wrk[iz];
                        }
                    }
                }
            }
            if(f2 != 0) {
                f = f1/f2;
            }else {
                f = f1;
            }
            //Residual
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = GAUSS(hx, 0.1*nhx);
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*GAUSS(hz, 0.1*nhz);
                    for (ix=0; ix<nx; ix++){
                        for (iz=1; iz<nz-1; iz++){
                            wrk[iz] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }	
                        for (iz=0; iz<nz-1; iz++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = (f1/(f2*f2))*G2*wrk[iz] - (((hx*hx) +  (hz*hz))*1.0/f2)*wrk[iz]; 
                        }	
                        imagedata[ki2D(ix,nz-1,ihx,ihz)] = 0.0;
                    }
                }
            }
            break;
        default:
            f = 0;
            break;
    }

    // Free work array
    free(wrk);

    image->setImagefile(imageresfile);
    image->write();

    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->output(MISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&f, 1, 0);
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Wemva<float>;
template class Wemva<double>;

template class WemvaAcoustic2D<float>;
template class WemvaAcoustic2D<double>;

template class WemvaElastic2D<float>;
template class WemvaElastic2D<double>;

}
