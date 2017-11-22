#include "inversion.h"


namespace rockseis {

// =============== ABSTRACT INVERSION CLASS =============== //
template<typename T>
Inversion<T>::Inversion() {
    //Set default parameters
    fs = false;
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    snapmethod = OPTIMAL; 
    misfit_type = DIFFERENCE;
    paramtype = PAR_GRID;
    dtx = -1;
    dty = -1;
    dtz = -1;
    fnorm = 0.0;
    if(createLog() == INV_ERR)
    {
        rs_error("Inversion<T>::Inversion(): Error creating logfile for writting.");
    }
    if(createProglog() == INV_ERR)
    {
        rs_error("Inversion<T>::Inversion(): Error creating progress logfile for writting.");
    }

}

template<typename T>
Inversion<T>::Inversion(MPImodeling *_mpi) {
    mpi = _mpi;

    //Set default parameters
    fs = false;
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    snapmethod = OPTIMAL; 
    misfit_type = DIFFERENCE;
    paramtype = PAR_GRID;
    dtx = -1;
    dty = -1;
    dtz = -1;
    fnorm = 0.0;
    if(createLog() == INV_ERR)
    {
        rs_error("Inversion<T>::Inversion(): Error creating logfile for writting.");
    }
    if(createProglog() == INV_ERR)
    {
        rs_error("Inversion<T>::Inversion(): Error creating progress logfile for writting.");
    }

}

template<typename T>
void Inversion<T>::normalize(double *v, double *f, int n){
	int i;
    *f /= abs(fnorm);
	for(i=0; i<n; i++) {
		v[i] /= abs(fnorm);
	}
}

template<typename T>
double Inversion<T>::vector_norm(double *v, const int type, const int n){
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
			if(abs(v[i]) >= norm) {
				norm = abs(v[i]);
			}
		}
	}

	return norm;
}

template<typename T>
bool Inversion<T>::createLog(){
	logfile = LOGFILE;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return INV_ERR;
	}else{
		Flog.close();
		return INV_OK;
	}
}

template<typename T>
bool Inversion<T>::createProglog(){
	progresslogfile = PROGLOGFILE;
	Flog.open(progresslogfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return INV_ERR;
	}else{
		Flog.close();
		return INV_OK;
	}
}

template<typename T>
void Inversion<T>::writeLog(std::string text){
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
void Inversion<T>::writeProgress(std::string text){
    if(!progresslogfile.empty()){
        Flog.open(progresslogfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Inversion<T>::createResult(){
    struct stat s;
    // Checking if result folder is present, and creates it if not
    if(stat(RESULTDIR,&s) != 0) {
        int mkdir_return = mkdir(RESULTDIR,0777);
        if(mkdir_return != 0) rs_error("Not able to create result directory: ", RESULTDIR);
    }
}




template<typename T>
Inversion<T>::~Inversion() {
    //Do nothing
}

// =============== 2D ACOUSTIC INVERSION CLASS =============== //
//
template<typename T>
InversionAcoustic2D<T>::InversionAcoustic2D() {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;
    krho = 1.0;
    ksource = 1.0;

    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;

    update_vp = true;
    update_rho = false;
    update_source = false;
}

template<typename T>
InversionAcoustic2D<T>::InversionAcoustic2D(MPImodeling *mpi): Inversion<T>(mpi) {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;
    krho = 1.0;
    ksource = 1.0;
    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;

    update_vp = true;
    update_rho = false;
    update_source = false;
}

template<typename T>
InversionAcoustic2D<T>::~InversionAcoustic2D() {
    //Do nothing
}

template<typename T>
void InversionAcoustic2D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> shot2D;
    std::shared_ptr<rockseis::Data2D<T>> shot2Di;
    std::shared_ptr<rockseis::Data2D<T>> shotmod2D;
    std::shared_ptr<rockseis::Data2D<T>> shotmod2Di;
    std::shared_ptr<rockseis::Data2D<T>> shotres2D;
    std::shared_ptr<rockseis::Data2D<T>> shotres2Di;
    std::shared_ptr<rockseis::Data2D<T>> shotweight2D;
    std::shared_ptr<rockseis::Data2D<T>> shotweight2Di;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::Image2D<T>> rhograd;
    std::shared_ptr<rockseis::Data2D<T>> wavgrad;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> gmodel (new rockseis::ModelAcoustic2D<T>(Vpfile, Rhofile, this->getLpml() ,this->getFs()));
    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> lmodel (new rockseis::ModelAcoustic2D<T>(Vpfile, Rhofile, this->getLpml() ,this->getFs()));

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

        // Wavelet gradient
        wavgrad = std::make_shared<rockseis::Data2D<T>>(1, source->getNt(), source->getDt(), 0.0);
        wavgrad->setFile(Wavgradfile);
        wavgrad->createEmpty(ngathers);

        // Misfit file creation
        Fmisfit->output(Misfitfile);
        Fmisfit->setN(1,ngathers);
        Fmisfit->setD(1,1.0);
        Fmisfit->setData_format(sizeof(T));
        Fmisfit->createEmpty();
        Fmisfit->close();

        // Create a data class for the recorded data
        std::shared_ptr<rockseis::Data2D<T>> shot2D (new rockseis::Data2D<T>(Precordfile));
        // Create modelling and residual data files
        shotmod2D = std::make_shared<rockseis::Data2D<T>>(1, shot2D->getNt(), shot2D->getDt(), shot2D->getOt());
        shotmod2D->setFile(Pmodelledfile);
        shotmod2D->createEmpty(shot2D->getNtrace());

        shotres2D = std::make_shared<rockseis::Data2D<T>>(1, shot2D->getNt(), shot2D->getDt(), shot2D->getOt());
        shotres2D->setFile(Presidualfile);
        shotres2D->createEmpty(shot2D->getNtrace());
        
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

		// Perform work in parallel
		mpi->performWork();

        // Image
        vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile, gmodel, 1, 1);
        vpgrad->createEmpty();

        rhograd = std::make_shared<rockseis::Image2D<T>>(Rhogradfile, gmodel, 1, 1);
        rhograd->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            remove_file(Vpgradfile + "-" + std::to_string(i));
            rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
            remove_file(Rhogradfile + "-" + std::to_string(i));
        }

		//Clear work vector 
		mpi->clearWork();
	}
    else {
        /* Slave */
        std::shared_ptr<rockseis::FwiAcoustic2D<T>> fwi;
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

                // Get the weight
                if(dataweight){
                    Sort->setDatafile(Dataweightfile);
                    shotweight2D = Sort->get2DGather(work.id);
                }

                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot2D, shot2Di);
                shot2Di->makeMap(lmodel->getGeom(), GMAP);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiAcoustic2D<T>>(lmodel, source, shot2Di, this->getOrder(), this->getSnapinc());

                // Create modelled and residual data objects 
                shotmod2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotmod2D->copyCoords(shot2D);
                shotmod2D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDatamodP(shotmod2D);
                shotres2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotres2D->copyCoords(shot2D);
                shotres2D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDataresP(shotres2D);

                // Interpolate weight
                if(dataweight){
                    shotweight2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(shotweight2D, shotweight2Di);
                    shotweight2Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweight(shotweight2Di);
                }
                
                // Setting misfit type
                fwi->setMisfit_type(this->getMisfit_type());

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                rhograd = std::make_shared<rockseis::Image2D<T>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                // Setting up gradient objects in fwi class
                fwi->setVpgrad(vpgrad);
                fwi->setRhograd(rhograd);

                wavgrad = std::make_shared<rockseis::Data2D<T>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
                wavgrad->setField(rockseis::PRESSURE);
                // Copy geometry
                wavgrad->copyCoords(source);
                wavgrad->makeMap(lmodel->getGeom(), SMAP);
                fwi->setWavgrad(wavgrad);

                // Setting Snapshot file 
                fwi->setSnapfile(Psnapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                fwi->setNcheck(this->getNsnaps());
                fwi->setIncore(this->getIncore());

                // Set logfile
                fwi->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        fwi->run();
                        break;
                    case rockseis::OPTIMAL:
                        fwi->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output gradients
                vpgrad->write();
                rhograd->write();
                wavgrad->putTrace(Wavgradfile, work.id);

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fwi->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                shotmod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, shot2D->getNt(), shot2D->getDt(), shot2D->getOt());
                shotmod2Di->setFile(Pmodelledfile);
                interp->interp(shotmod2D, shotmod2Di);
                Sort->put2DGather(shotmod2Di, work.id);

                shotres2Di = std::make_shared<rockseis::Data2D<T>>(ntr, shot2D->getNt(), shot2D->getDt(), shot2D->getOt());
                shotres2Di->setFile(Presidualfile);
                interp->interp(shotres2D, shotres2Di);
                Sort->put2DGather(shotres2Di, work.id);

                
                // Reset all classes
                shot2D.reset();
                shot2Di.reset();
                shotmod2D.reset();
                shotmod2Di.reset();
                shotres2D.reset();
                shotres2Di.reset();
                lmodel.reset();
                vpgrad.reset();
                rhograd.reset();
                wavgrad.reset();
                fwi.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void InversionAcoustic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	T rhosum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;
    std::string rhogradfile;

    if(Mutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
        rhogradfile = RHOGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad (new rockseis::ModelAcoustic2D<T>(vpgradfile, rhogradfile, this->getLpml() ,this->getFs()));

	// Read model
	grad->readModel();
	
	T *vpgrad, *rhograd;
	vpgrad = grad->getVp();
	rhograd = grad->getR();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl2D<T>> spline (new rockseis::Bspl2D<T>(grad->getNx(), grad->getNz(), grad->getDx(), grad->getDz(), this->getDtx(), this->getDtz(), 3, 3));
    int nc = spline->getNc();

	/* Allocating projection arrays */
	float *vpproj= (float *) calloc(nc, sizeof(float));
	if(vpproj==NULL){
		rs_error("InversionAcoustic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
	}
	float *rhoproj= (float *) calloc(nc, sizeof(float));
	if(rhoproj==NULL){
		rs_error("InversionAcoustic2D<T>::runBsproj(): Not enough memory to allocate projection array (rhoproj)");
	}

    if(mpi->getRank() == 0) {
		// Master

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
			rs_error("InversionAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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

		for(long int i=0; i< nc; i++){
			global_stack[i] = 0.0;
		}
		/* Starting reduce operation */
		MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   

		/* Output spline */
        Fout->output(RHOPROJGRADFILE);
        Fout->setN(1,nc);
        Fout->setD(1,1.0);
        Fout->setData_format(sizeof(float));
        Fout->write(global_stack, nc, 0);
        Fout->close();

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
				rhosum = 0.0;
				for(long int i=0; i<grad->getNx()*grad->getNz(); i++){
						vpsum += wrk[i]*vpgrad[i];
						rhosum += wrk[i]*rhograd[i];
				}
				vpproj[work.id]=vpsum;
				rhoproj[work.id]=rhosum;
				c[work.id]=0.0; // Reset coefficient to 0
			}

			// Send result back
			work.status = WORK_FINISHED;
			mpi->sendResult(work);		
		}

		global_stack= (float *) calloc(nc, sizeof(float));
		if(global_stack==NULL){
			rs_error("InversionAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
        MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        for(long int i=0; i< nc; i++){
            global_stack[i] = 0.0;
        }
        MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

	   }
    // Free variables
    free(global_stack);
    free(vpproj);
    free(rhoproj);
}

template<typename T>
int InversionAcoustic2D<T>::setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile)
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
    int N, Ns;
    int Npar = 0;
    switch(this->getParamtype()){
        case PAR_GRID:
            N = (model_in->getGeom())->getNtot();
            Ns = source_in->getNt();
            Npar = 2*N + Ns;
            break;
        case PAR_BSPLINE:
             spline = std::make_shared<rockseis::Bspl2D<T>>(model_in->getNx(), model_in->getNz(), model_in->getDx(), model_in->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N=spline->getNc();
            Ns = source_in->getNt();
            Npar = 2*N + Ns;
            break;
        default:
            rs_error("InversionAcoustic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void InversionAcoustic2D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcels (new rockseis::Data2D<T>(SOURCELSFILE));
    lsmodel->readModel();
    sourcels->read();
    if(update_vp){
        name = dir + "/" + VP_UP + "-" + std::to_string(iter);
        lsmodel->setVpfile(name);
        lsmodel->writeVp();
    }
    if(update_rho){
        name = dir + "/" + RHO_UP + "-" + std::to_string(iter);
        lsmodel->setRfile(name);
        lsmodel->writeR();
    }
    if(update_source){
        name = dir + "/" + RHO_UP + "-" + std::to_string(iter);
        name = dir + "/" + SOURCE_UP + "-" + std::to_string(iter);
        sourcels->setFile(name);
        sourcels->write();
    }
}

template<typename T>
void InversionAcoustic2D<T>::saveLinesearch(double *x)
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
    T *rhomutedata;
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
        if(N != Nmute) rs_error("InversionAcoustic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readModel();
        vpmutedata = mute->getVp();
        rhomutedata = mute->getR();
        N = (lsmodel->getGeom())->getNtot();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        rhomutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
            rhomutedata[i] = 1.0;
        }
    }

    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
                vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
                rhols[i] = rho0[i] + x[N+i]*rhomutedata[i]*krho;
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = wav0[i] + x[2*N+i]*ksource;
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
            }
            for(i=0; i< N; i++)
            {
                c[i] = x[i+N];
            }
            spline->bisp();
            mod = spline->getMod();

            for(i=0; i< Nmod; i++)
            {
                rhols[i] = rho0[i] + mod[i]*rhomutedata[i]*krho;
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = wav0[i] + x[2*N+i]*ksource;
            }
            lssource->write();
            break;
        default:
            rs_error("InversionAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Mutefile.empty()){
        free(vpmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionAcoustic2D<T>::readMisfit(double *f)
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

    Fmisfit->input(RHOREGMISFITFILE);
    Fmisfit->read(&val, 1, 0); 
    *f += reg_alpha[1]*val;
    Fmisfit->close();
}

template<typename T>
void InversionAcoustic2D<T>::readGrad(double *g)
{
    int i;
    int N,Ns;
    float *g_in;
    T *gvp, *grho, *gwav;
    std::string vpgradfile;
    std::string rhogradfile;
    if(Mutefile.empty()){
        vpgradfile = VPGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelAcoustic2D<T>> modelgrad (new rockseis::ModelAcoustic2D<T>(vpgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcegrad (new rockseis::Data2D<T>(SOURCEGRADFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readModel();
            N = (modelgrad->getGeom())->getNtot();
            gvp = modelgrad->getVp(); 
            grho = modelgrad->getR(); 
            for(i=0; i< N; i++)
            {
                g[i] = gvp[i]*kvp;
                g[N+i] = grho[i]*krho;
            }
            Ns = sourcegrad->getNt();
            sourcegrad->read();
            gwav = sourcegrad->getData();
            for(i=0; i< Ns; i++)
            {
                g[2*N+i] = gwav[i]*ksource;

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
            Fgrad->input(RHOPROJGRADFILE);
            Fgrad->read(&g_in[N], N, 0);
            Fgrad->close();
            for(i=0; i< N; i++)
            {
                g[i] = g_in[i]*kvp;
                g[N+i] = g_in[N+i]*krho;
            }
            Ns = sourcegrad->getNt();
            sourcegrad->read();
            gwav = sourcegrad->getData();
            for(i=0; i< Ns; i++)
            {
                g[2*N+i] = gwav[i]*ksource;

            }
            // Free temporary array
            free(g_in);
            break;
        default:
            rs_error("InversionAcoustic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void InversionAcoustic2D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADFILE, RHOGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPREGGRADFILE, RHOREGGRADFILE, 1 ,0);

    // Read gradients
    grad->readModel();
    reggrad->readModel();
    T *vp, *rho, *regvp, *regrho;
    vp = grad->getVp(); 
    rho = grad->getR(); 
    regvp = reggrad->getVp(); 
    regrho = reggrad->getR(); 
    int i;
    int N;

    N = (grad->getGeom())->getNtot();
    // Compute 
    for(i=0; i< N; i++)
    {
        vp[i] = vp[i] + reg_alpha[0]*regvp[i];
        rho[i] = rho[i] + reg_alpha[1]*regrho[i];
    }
    grad->setVpfile(VPGRADCOMBFILE);
    grad->setRfile(RHOGRADCOMBFILE);
    grad->writeModel();
}


template<typename T>
void InversionAcoustic2D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute (new rockseis::ModelAcoustic2D<T>(Mutefile, Mutefile, 1 ,0));

        // Mute model and write
        model->readModel();
        mute->readModel();
        T *vp, *rho, *vpmute, *rhomute;
        vp = model->getVp(); 
        rho = model->getR(); 
        vpmute = mute->getVp(); 
        rhomute = mute->getR(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
            rho[i] = rho[i]*rhomute[i];
        }
        model->setVpfile(VPGRADMUTEFILE);
        model->setRfile(RHOGRADMUTEFILE);
        model->writeModel();
    }
}

template<typename T>
void InversionAcoustic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), 1, model->getNz(), model->getDx(), 1.0, model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx, *dvpdz;
    double *drhodx, *drhodz;
    T *vpgrad, *rhograd;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    drhodx = (double *) calloc(Nmod, sizeof(double));
    drhodz = (double *) calloc(Nmod, sizeof(double));
    gwrk = (double *) calloc(Nmod, sizeof(double));
    model->readModel();
    model->setVpfile(VPREGGRADFILE);
    model->setRfile(RHOREGGRADFILE);
    vpgrad = model->getVp();
    rhograd = model->getR();
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
            der->ddx_fw(&x[N]);
            for(i=0; i< N; i++)
            {
                drhodx[i] = df[i]*krho;
            }
            der->ddz_fw(&x[N]);
            for(i=0; i< N; i++)
            {
                drhodz[i] = df[i]*krho;
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
            for(i=0; i< N; i++)
            {
                c[i] = x[i+N];
            }
            spline->bisp();
            mod = spline->getMod();
            der->ddx_fw(mod);
            for(i=0; i< Nmod; i++)
            {
                drhodx[i] = df[i]*krho;
            }
            der->ddz_fw(mod);
            for(i=0; i< Nmod; i++)
            {
                drhodz[i] = df[i]*krho;
            }

            break;
        default:
            rs_error("InversionAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
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

    T frho = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i];
        M = sqrt(M);
        frho += M;
    }

    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = drhodx[i]/M;
    }
    der->ddx_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = drhodz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] -= df[i];
    }

    /* Write out misfit */
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->output(VPREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&fvp,1,0);
    Fmisfit->close();

    Fmisfit->output(RHOREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&frho,1,0);
    Fmisfit->close();

    /* Write out gradient */
    model->writeModel();

    // Free variables
    free(dvpdx);
    free(dvpdz);
    free(drhodx);
    free(drhodz);
    free(gwrk);
}

// =============== 2D ELASTIC INVERSION CLASS =============== //
//
template<typename T>
InversionElastic2D<T>::InversionElastic2D() {
    // Set default parameters
    apertx = -1;
    sourcetype = 0;

    kvp = 1.0;
    kvs = 1.0;
    krho = 1.0;
    ksource = 1.0;

    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_alpha[2]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;
    reg_eps[2]=1e-3;

    update_vp = true;
    update_vs = true;
    update_rho = false;
    update_source = false;
}

template<typename T>
InversionElastic2D<T>::InversionElastic2D(MPImodeling *mpi): Inversion<T>(mpi) {
    // Set default parameters
    apertx = -1;
    sourcetype = 0;

    kvp = 1.0;
    kvs = 1.0;
    krho = 1.0;
    ksource = 1.0;
    reg_alpha[0]=0.0;
    reg_alpha[1]=0.0;
    reg_alpha[2]=0.0;
    reg_eps[0]=1e-3;
    reg_eps[1]=1e-3;
    reg_eps[2]=1e-3;

    update_vp = true;
    update_vs = true;
    update_rho = false;
    update_source = false;
}

template<typename T>
InversionElastic2D<T>::~InversionElastic2D() {
    //Do nothing
}

template<typename T>
void InversionElastic2D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> Vxdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Vxdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Vzdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Vzdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Vxdatamod2D;
    std::shared_ptr<rockseis::Data2D<T>> Vxdatamod2Di;
    std::shared_ptr<rockseis::Data2D<T>> Vxdatares2D;
    std::shared_ptr<rockseis::Data2D<T>> Vxdatares2Di;
    std::shared_ptr<rockseis::Data2D<T>> Vzdatamod2D;
    std::shared_ptr<rockseis::Data2D<T>> Vzdatamod2Di;
    std::shared_ptr<rockseis::Data2D<T>> Vzdatares2D;
    std::shared_ptr<rockseis::Data2D<T>> Vzdatares2Di;
    std::shared_ptr<rockseis::Data2D<T>> shotweight2D;
    std::shared_ptr<rockseis::Data2D<T>> shotweight2Di;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::Image2D<T>> vsgrad;
    std::shared_ptr<rockseis::Image2D<T>> rhograd;
    std::shared_ptr<rockseis::Data2D<T>> wavgrad;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Vxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelElastic2D<T>> gmodel (new rockseis::ModelElastic2D<T>(Vpfile, Vsfile, Rhofile, this->getLpml() ,this->getFs()));
    // Create a local model class
	std::shared_ptr<rockseis::ModelElastic2D<T>> lmodel (new rockseis::ModelElastic2D<T>(Vpfile, Vsfile, Rhofile, this->getLpml() ,this->getFs()));

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

		// Wavelet gradient
		wavgrad = std::make_shared<rockseis::Data2D<T>>(1, source->getNt(), source->getDt(), 0.0);
		wavgrad->setFile(Wavgradfile);
		wavgrad->createEmpty(ngathers);

		// Misfit file creation
		Fmisfit->output(Misfitfile);
		Fmisfit->setN(1,ngathers);
		Fmisfit->setD(1,1.0);
		Fmisfit->setData_format(sizeof(T));
		Fmisfit->createEmpty();
		Fmisfit->close();

		// Create a data class for the recorded data in order to get parameters from file
		std::shared_ptr<rockseis::Data2D<T>> Vxdata2D (new rockseis::Data2D<T>(Vxrecordfile));

		// Create a data class for the recorded data
		Vxdatamod2D = std::make_shared<rockseis::Data2D<T>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
		Vxdatamod2D->setFile(Vxmodelledfile);
		Vxdatamod2D->createEmpty(Vxdata2D->getNtrace());
		Vxdatares2D = std::make_shared<rockseis::Data2D<T>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
		Vxdatares2D->setFile(Vxresidualfile);
		Vxdatares2D->createEmpty(Vxdata2D->getNtrace());

		Vzdatamod2D = std::make_shared<rockseis::Data2D<T>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
		Vzdatamod2D->setFile(Vzmodelledfile);
		Vzdatamod2D->createEmpty(Vxdata2D->getNtrace());
		Vzdatares2D = std::make_shared<rockseis::Data2D<T>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
		Vzdatares2D->setFile(Vzresidualfile);
		Vzdatares2D->createEmpty(Vxdata2D->getNtrace());

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

        rhograd = std::make_shared<rockseis::Image2D<T>>(Rhogradfile, gmodel, 1, 1);
        rhograd->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            if(update_vp){
                vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
                remove_file(Vpgradfile + "-" + std::to_string(i));
            }
            if(update_vs){
                vsgrad->stackImage(Vsgradfile + "-" + std::to_string(i));
                remove_file(Vsgradfile + "-" + std::to_string(i));
            }
            if(update_rho){
                rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
                remove_file(Rhogradfile + "-" + std::to_string(i));
            }
        }

		//Clear work vector 
		mpi->clearWork();
	}
    else {
        /* Slave */
        std::shared_ptr<rockseis::FwiElastic2D<T>> fwi;
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
                Sort->setDatafile(Vxrecordfile);
                Vxdata2D = Sort->get2DGather(work.id);
                size_t ntr = Vxdata2D->getNtrace();

                Sort->setDatafile(Vzrecordfile);
                Vzdata2D = Sort->get2DGather(work.id);

                // Get the weight
                if(dataweight){
                    Sort->setDatafile(Dataweightfile);
                    shotweight2D = Sort->get2DGather(work.id);
                }

                lmodel = gmodel->getLocal(Vxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Vxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
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

                // Interpolate shot
                Vxdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vxdata2D, Vxdata2Di);
                Vxdata2Di->makeMap(lmodel->getGeom(), GMAP);
                Vxdata2Di->setField(rockseis::VX);

                Vzdata2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vzdata2D, Vzdata2Di);
                Vzdata2Di->makeMap(lmodel->getGeom(), GMAP);
                Vzdata2Di->setField(rockseis::VZ);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic2D<T>>(lmodel, source, Vxdata2Di, Vzdata2Di, this->getOrder(), this->getSnapinc());

                // Create modelled and residual data objects 
                Vxdatamod2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatamod2D->copyCoords(Vxdata2D);
                Vxdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatamod2D->setField(rockseis::VX);
                fwi->setDatamodVx(Vxdatamod2D);
                Vxdatares2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatares2D->copyCoords(Vxdata2D);
                Vxdatares2D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatares2D->setField(rockseis::VX);
                fwi->setDataresVx(Vxdatares2D);

                Vzdatamod2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatamod2D->copyCoords(Vzdata2D);
                Vzdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatamod2D->setField(rockseis::VZ);
                fwi->setDatamodVz(Vzdatamod2D);
                Vzdatares2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatares2D->copyCoords(Vzdata2D);
                Vzdatares2D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatares2D->setField(rockseis::VZ);
                fwi->setDataresVz(Vzdatares2D);

                // Interpolate weight
                if(dataweight){
                    shotweight2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(shotweight2D, shotweight2Di);
                    shotweight2Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweight(shotweight2Di);
                }

                // Setting misfit type
                fwi->setMisfit_type(this->getMisfit_type());

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                fwi->setVpgrad(vpgrad);
                if(update_vs || update_rho){
                    vsgrad = std::make_shared<rockseis::Image2D<T>>(Vsgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                    fwi->setVsgrad(vsgrad);
                }
                if(update_rho){
                    rhograd = std::make_shared<rockseis::Image2D<T>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                    fwi->setRhograd(rhograd);
                }

                if(update_source){
                    wavgrad = std::make_shared<rockseis::Data2D<T>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
                    // Copy geometry
                    wavgrad->copyCoords(source);
                    wavgrad->makeMap(lmodel->getGeom(), SMAP);
                    switch(this->getSourcetype()){
                        case 0:
                            wavgrad->setField(PRESSURE);
                            break;
                        case 1:
                            wavgrad->setField(VX);
                            break;
                        case 3:
                            wavgrad->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown wavgrad type: ", std::to_string(this->getSourcetype()));
                            break;
                    }

                    fwi->setWavgrad(wavgrad);
                }

                // Setting Snapshot file 
                fwi->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                fwi->setNcheck(this->getNsnaps());
                fwi->setIncore(this->getIncore());

                // Set logfile
                fwi->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(this->getSnapmethod()){
                    case rockseis::FULL:
                        fwi->run();
                        break;
                    case rockseis::OPTIMAL:
                        fwi->run_optimal();
                        break;
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output gradients
                if(update_vp){
                    vpgrad->write();
                }
                if(update_vs){
                    vsgrad->write();
                }
                if(update_rho){
                    rhograd->write();
                }
                if(update_source){
                    wavgrad->putTrace(Wavgradfile, work.id);
                }

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fwi->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                Vxdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vxdatamod2Di->setFile(Vxmodelledfile);
                interp->interp(Vxdatamod2D, Vxdatamod2Di);
                Sort->put2DGather(Vxdatamod2Di, work.id);

                Vxdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vxdatares2Di->setFile(Vxresidualfile);
                interp->interp(Vxdatares2D, Vxdatares2Di);
                Sort->put2DGather(Vxdatares2Di, work.id);

                Vzdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vzdatamod2Di->setFile(Vzmodelledfile);
                interp->interp(Vzdatamod2D, Vzdatamod2Di);
                Sort->put2DGather(Vzdatamod2Di, work.id);

                Vzdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vzdatares2Di->setFile(Vzresidualfile);
                interp->interp(Vzdatares2D, Vzdatares2Di);
                Sort->put2DGather(Vzdatares2Di, work.id);

                // Reset all classes
                Vxdata2D.reset();
                Vxdata2Di.reset();
                Vxdatamod2D.reset();
                Vxdatamod2Di.reset();
                Vxdatares2D.reset();
                Vxdatares2Di.reset();
                Vzdata2D.reset();
                Vzdata2Di.reset();
                Vzdatamod2D.reset();
                Vzdatamod2Di.reset();
                Vzdatares2D.reset();
                Vzdatares2Di.reset();
                lmodel.reset();
                vpgrad.reset();
                vsgrad.reset();
                rhograd.reset();
                wavgrad.reset();
                fwi.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void InversionElastic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	T vssum = 0.0; // Sum over splines
	T rhosum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;
    std::string vsgradfile;
    std::string rhogradfile;

    if(Mutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
        vsgradfile = VSGRADCOMBFILE;
        rhogradfile = RHOGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelElastic2D<T>> grad (new rockseis::ModelElastic2D<T>(vpgradfile, vsgradfile, rhogradfile, this->getLpml() ,this->getFs()));
    
    //Read gradients
    grad->readModel();

	T *vpgrad, *vsgrad, *rhograd;
	vpgrad = grad->getVp();
	vsgrad = grad->getVs();
	rhograd = grad->getR();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl2D<T>> spline (new rockseis::Bspl2D<T>(grad->getNx(), grad->getNz(), grad->getDx(), grad->getDz(), this->getDtx(), this->getDtz(), 3, 3));
    int nc = spline->getNc();

    float *vpproj, *vsproj, *rhoproj;
    /* Allocating projection arrays */
    vpproj= (float *) calloc(nc, sizeof(float));
    if(vpproj==NULL){
        rs_error("InversionElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
    }
    vsproj= (float *) calloc(nc, sizeof(float));
    if(vsproj==NULL){
        rs_error("InversionElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vsproj)");
    }
    rhoproj= (float *) calloc(nc, sizeof(float));
    if(rhoproj==NULL){
        rs_error("InversionElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (rhoproj)");
    }

    if(mpi->getRank() == 0) {
		// Master

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
			rs_error("InversionElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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

        if(update_rho){
            /* Starting reduce operation */
            MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   
            /* Output spline */
            Fout->output(RHOPROJGRADFILE);
            Fout->setN(1,nc);
            Fout->setD(1,1.0);
            Fout->setData_format(sizeof(float));
            Fout->write(global_stack, nc, 0);
            Fout->close();
        }

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
                rhosum = 0.0;
                for(long int i=0; i<grad->getNx()*grad->getNz(); i++){
                    if(update_vp){
                        vpsum += wrk[i]*vpgrad[i];
                    }
                    if(update_vs){
                        vssum += wrk[i]*vsgrad[i];
                    }
                    if(update_rho){
                        rhosum += wrk[i]*rhograd[i];
                    }
                }
                vpproj[work.id]=vpsum;
                vsproj[work.id]=vssum;
                rhoproj[work.id]=rhosum;
                c[work.id]=0.0; // Reset coefficient to 0
            }

            // Send result back
            work.status = WORK_FINISHED;
            mpi->sendResult(work);		
        }

        global_stack= (float *) calloc(nc, sizeof(float));
        if(global_stack==NULL){
            rs_error("InversionElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
        }

        /* Starting reduce operation */
        if(update_vp){
            MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
        if(update_vs){
            MPI_Reduce(vsproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
        if(update_rho){
            MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
       }
    // Free allocated variables
    free(vpproj);
    free(vsproj);
    free(rhoproj);
    free(global_stack);
}
		
template<typename T>
int InversionElastic2D<T>::setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile)
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
            if(update_rho) Npar += N;
            if(update_source) Npar += Ns;
            break;
        case PAR_BSPLINE:
             spline = std::make_shared<rockseis::Bspl2D<T>>(model_in->getNx(), model_in->getNz(), model_in->getDx(), model_in->getDz(), this->getDtx(), this->getDtz(), 3, 3);
            N=spline->getNc();
            Ns = source_in->getNt();
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_vs) Npar += N;
            if(update_rho) Npar += N;
            if(update_source) Npar += Ns;
            break;
        default:
            rs_error("InversionElastic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void InversionElastic2D<T>::saveResults(int iter)
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
        name = dir + "/" + VP_UP + "-" + std::to_string(iter);
        lsmodel->setVpfile(name);
        lsmodel->writeVp();
    }
    if(update_vs){
        name = dir + "/" + VS_UP + "-" + std::to_string(iter);
        lsmodel->setVsfile(name);
        lsmodel->writeVs();
    }
    if(update_rho){
        name = dir + "/" + RHO_UP + "-" + std::to_string(iter);
        lsmodel->setRfile(name);
        lsmodel->writeR();
    }
    if(update_source){
        name = dir + "/" + RHO_UP + "-" + std::to_string(iter);
        name = dir + "/" + SOURCE_UP + "-" + std::to_string(iter);
        sourcels->setFile(name);
        sourcels->write();
    }
}

template<typename T>
void InversionElastic2D<T>::saveLinesearch(double *x)
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
    T *rhomutedata;
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
        if(N != Nmute) rs_error("InversionElastic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readModel();
        vpmutedata = mute->getVp();
        vsmutedata = mute->getVs();
        rhomutedata = mute->getR();
        N = (lsmodel->getGeom())->getNtot();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        vsmutedata = (T *) calloc(N, sizeof(T)); 
        rhomutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
            vsmutedata[i] = 1.0;
            rhomutedata[i] = 1.0;
        }
    }

    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            Npar = 0;
            for(i=0; i< N; i++)
            {
                if(update_vp){
                    vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
                }else{
                    vpls[i] = vp0[i];
                }
            }
            Npar += N;
            for(i=0; i< N; i++)
            {
                if(update_vs){
                    vsls[i] = vs0[i] + x[Npar+i]*vsmutedata[i]*kvs;
                }else{
                    vsls[i] = vs0[i];
                }
            }
            Npar += N;
            for(i=0; i< N; i++)
            {
                if(update_rho){
                    rhols[i] = rho0[i] + x[Npar+i]*rhomutedata[i]*krho;
                }else{
                    rhols[i] = rho0[i];
                }
            }
            Npar += N;
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
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    c[i] = x[Npar+i];
                }
                spline->bisp();
                mod = spline->getMod();

                for(i=0; i< Nmod; i++)
                {
                    rhols[i] = rho0[i] + mod[i]*rhomutedata[i]*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    rhols[i] = rho0[i];
                }
            }
            lsmodel->writeModel();
            break;
        default:
            rs_error("InversionAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }

    // Source wavelet
    Ns = lssource->getNt();
    for(i=0; i< Ns; i++)
    {
        if(update_source){
            wavls[i] = wav0[i] + x[Npar+i]*ksource;
        }else{
            wavls[i] = wav0[i];
        }
    }
    lssource->write();

    // Free allocated arrays
    if(Mutefile.empty()){
        free(vpmutedata);
        free(vsmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionElastic2D<T>::readMisfit(double *f)
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

    Fmisfit->input(RHOREGMISFITFILE);
    Fmisfit->read(&val, 1, 0); 
    *f += reg_alpha[2]*val;
    Fmisfit->close();
}

template<typename T>
void InversionElastic2D<T>::readGrad(double *g)
{
    int i;
    int N,Ns,Npar=0;
    float *g_in;
    T *gvp, *gvs, *grho, *gwav;
    std::string vpgradfile;
    std::string vsgradfile;
    std::string rhogradfile;
    if(Mutefile.empty()){
        vpgradfile = VPGRADFILE;
        vsgradfile = VSGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelElastic2D<T>> modelgrad (new rockseis::ModelElastic2D<T>(vpgradfile, vsgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcegrad (new rockseis::Data2D<T>(SOURCEGRADFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readModel();
            N = (modelgrad->getGeom())->getNtot();
            Npar = 0;
            gvp = modelgrad->getVp(); 
            gvs = modelgrad->getVs(); 
            grho = modelgrad->getR(); 
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
            if(update_rho){
                for(i=0; i < N; i++)
                {
                    g[Npar+i] = grho[i]*krho;
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
            if(update_rho){
                Fgrad->input(RHOPROJGRADFILE);
                Fgrad->read(&g_in[0], N, 0);
                Fgrad->close();
                for(i=0; i< N; i++)
                {
                    g[Npar+i] = g_in[i]*krho;
                }
                Npar += N;
            }

            // Free temporary array
            free(g_in);
            break;
        default:
            rs_error("InversionElastic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
    if(update_source){
        Ns = sourcegrad->getNt();
        sourcegrad->read();
        gwav = sourcegrad->getData();
        for(i=0; i< Ns; i++)
        {
            g[Npar+i] = gwav[i]*ksource;

        }
    }
}

template<typename T>
void InversionElastic2D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelElastic2D<T>> grad;
    std::shared_ptr<rockseis::ModelElastic2D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADFILE, VSGRADFILE, RHOGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelElastic2D<T>>(VPREGGRADFILE, VSREGGRADFILE, RHOREGGRADFILE, 1 ,0);

    // Read gradients
    grad->readModel();
    reggrad->readModel();
    T *vp, *vs, *rho, *regvp, *regvs, *regrho;
    vp = grad->getVp(); 
    vs = grad->getVs(); 
    rho = grad->getR(); 
    regvp = reggrad->getVp(); 
    regvs = reggrad->getVs(); 
    regrho = reggrad->getR(); 
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
        if(update_rho){
            rho[i] = rho[i] + reg_alpha[2]*regrho[i];
        }
    }
    grad->setVpfile(VPGRADCOMBFILE);
    grad->setVsfile(VSGRADCOMBFILE);
    grad->setRfile(RHOGRADCOMBFILE);
    grad->writeModel();
}


template<typename T>
void InversionElastic2D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelElastic2D<T>> model;
        model = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADCOMBFILE, VSGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelElastic2D<T>> mute (new rockseis::ModelElastic2D<T>(Mutefile, Mutefile, Mutefile, 1 ,0));

        // Mute model and write
        model->readModel();
        mute->readModel();
        T *vp, *vs, *rho, *vpmute, *vsmute, *rhomute;
        vp = model->getVp(); 
        vs = model->getVs(); 
        rho = model->getR(); 
        vpmute = mute->getVp(); 
        vsmute = mute->getVs(); 
        rhomute = mute->getR(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
            vs[i] = vs[i]*vsmute[i];
            rho[i] = rho[i]*rhomute[i];
        }
        model->setVpfile(VPGRADMUTEFILE);
        model->setVsfile(VSGRADMUTEFILE);
        model->setRfile(RHOGRADMUTEFILE);
        model->writeModel();
    }
}

template<typename T>
void InversionElastic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic2D<T>> model (new rockseis::ModelElastic2D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), 1, model->getNz(), model->getDx(), 1.0, model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx, *dvpdz;
    double *dvsdx, *dvsdz;
    double *drhodx, *drhodz;
    T *vpgrad, *vsgrad, *rhograd;
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
    drhodx = (double *) calloc(Nmod, sizeof(double));
    drhodz = (double *) calloc(Nmod, sizeof(double));
    gwrk = (double *) calloc(Nmod, sizeof(double));
    model->readModel();
    model->setVpfile(VPREGGRADFILE);
    model->setVsfile(VSREGGRADFILE);
    model->setRfile(RHOREGGRADFILE);
    vpgrad = model->getVp();
    vsgrad = model->getVs();
    rhograd = model->getR();
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
                der->ddx_fw(x);
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
            if(update_rho){
                der->ddx_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    drhodx[i] = df[i]*krho;
                }
                der->ddz_fw(&x[N]);
                for(i=0; i< N; i++)
                {
                    drhodz[i] = df[i]*krho;
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
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    c[i] = x[Npar+i];
                }
                spline->bisp();
                mod = spline->getMod();
                der->ddx_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    drhodx[i] = df[i]*krho;
                }
                der->ddz_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    drhodz[i] = df[i]*krho;
                }
                Npar += N;
            }

            break;
        default:
            rs_error("InversionElastic2D<T>::computeRegularization(): Unknown parameterisation."); 
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


    // Computing misfit
    T frho = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i];
        M = sqrt(M);
        frho += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i] + reg_eps[2]*reg_eps[2];
        M = sqrt(M);
        gwrk[i] = drhodx[i]/M;
    }
    der->ddx_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhodz[i]*drhodz[i] + reg_eps[2]*reg_eps[2];
        M = sqrt(M);
        gwrk[i] = drhodz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] -= df[i];
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

    Fmisfit->output(RHOREGMISFITFILE);
    Fmisfit->setN(1,1);
    Fmisfit->setD(1,1.0);
    Fmisfit->setData_format(sizeof(T));
    Fmisfit->write(&frho,1,0);
    Fmisfit->close();

    /* Write out gradient */
    model->writeModel();

    // Free variables
    free(dvpdx);
    free(dvpdz);
    free(dvsdx);
    free(dvsdz);
    free(drhodx);
    free(drhodz);
    free(gwrk);
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Inversion<float>;
template class Inversion<double>;

template class InversionAcoustic2D<float>;
template class InversionAcoustic2D<double>;

template class InversionElastic2D<float>;
template class InversionElastic2D<double>;

}


