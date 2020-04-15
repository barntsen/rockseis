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
    noreverse = false;
    filter = false;
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
    noreverse = false;
    filter = false;
}

template<typename T>
void Inversion<T>::normalize(double *v, double *f, int n){
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
void Inversion<T>::un_normalize(double *v, double f, int n){
	int i;
	for(i=0; i<n; i++) {
		v[i] *= fabs(fnorm);
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
			if(fabs(v[i]) >= norm) {
				norm = fabs(v[i]);
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
    std::shared_ptr<rockseis::Image2D<T>> srcilum;
    std::shared_ptr<rockseis::Data2D<T>> wavgrad;
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

        if(this->srcilumset){
            srcilum = std::make_shared<rockseis::Image2D<T>>(Srcilumfile, gmodel, 1, 1);
            srcilum->createEmpty();
        }

        for(long int i=0; i<ngathers; i++) {
            vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            remove_file(Vpgradfile + "-" + std::to_string(i));
            rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
            remove_file(Rhogradfile + "-" + std::to_string(i));

            if(this->srcilumset){
                srcilum->stackImage(Srcilumfile + "-" + std::to_string(i));
                remove_file(Srcilumfile + "-" + std::to_string(i));
            }
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
                fwi->setFilter(this->getFilter());
                fwi->setAllfreqs(this->getFreqs());

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                rhograd = std::make_shared<rockseis::Image2D<T>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                if(this->srcilumset){
                    srcilum = std::make_shared<rockseis::Image2D<T>>(Srcilumfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                    fwi->setSrcilum(srcilum);
                }

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

                // Output gradients
                if(this->getFilter()){
                    wavgrad->apply_filter(this->getFreqs());
                }
                wavgrad->putTrace(Wavgradfile, work.id);

                vpgrad->write();
                rhograd->write();

                // Output ilumination
                if(this->srcilumset){
                    srcilum->write();
                }
                
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
                srcilum.reset();
                fwi.reset();

                // Send result back
                work.status = WORK_FINISHED;
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

    if(Modmutefile.empty()){
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


    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> sourcemute (new rockseis::Data2D<T>(Srcmutefile));
        Ns = sourcemute->getNt();
        if(Ns != source0->getNt()){
            rs_error("InversionAcoustic2D<T>::saveLinesearch(): Geometry in Srcmutefile does not match geometry in the Source file.");
        }
    }

    // If mute
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionAcoustic2D<T>::saveLinesearch(): Geometry in Modmutefile does not match geometry in the model.");
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
    if(Modmutefile.empty()){
        free(vpmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionAcoustic2D<T>::saveHessian(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> lssource (new rockseis::Data2D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute;

    // Write linesearch model
    lssource->read();
    lsmodel->readModel();
    T *vpls, *rhols, *wavls;
    T *c, *mod;
    T *vpmutedata;
    T *rhomutedata;
    vpls = lsmodel->getVp(); 
    rhols = lsmodel->getR(); 
    wavls = lssource->getData();
    lsmodel->setVpfile(VPHESSFILE);
    lsmodel->setRfile(RHOHESSFILE);
    lssource->setFile(SOURCEHESSFILE);
    int i;
    int N, Ns, Nmod;


    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> sourcemute (new rockseis::Data2D<T>(Srcmutefile));
        Ns = sourcemute->getNt();
        if(Ns != lssource->getNt()){
            rs_error("InversionAcoustic2D<T>::saveHessian(): Geometry in Srcmutefile does not match geometry in the Source file.");
        }
    }

    // If mute
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionAcoustic2D<T>::saveHessian(): Geometry in Modmutefile does not match geometry in the model.");
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
                vpls[i] = x[i]*vpmutedata[i]*kvp*kvp;
                rhols[i] = x[N+i]*rhomutedata[i]*krho*krho;
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = x[2*N+i]*ksource*ksource;
            }
            lssource->write();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            spline = std::make_shared<rockseis::Bspl2D<T>>(lsmodel->getNx(), lsmodel->getNz(), lsmodel->getDx(), lsmodel->getDz(), this->getDtx(), this->getDtz(), 3, 3);
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
                vpls[i] = mod[i]*vpmutedata[i]*kvp*kvp;
            }
            for(i=0; i< N; i++)
            {
                c[i] = x[i+N];
            }
            spline->bisp();
            mod = spline->getMod();

            for(i=0; i< Nmod; i++)
            {
                rhols[i] = mod[i]*rhomutedata[i]*krho*krho;
            }
            lsmodel->writeModel();
            Ns = lssource->getNt();
            for(i=0; i< Ns; i++)
            {
                wavls[i] = x[2*N+i]*ksource*ksource;
            }
            lssource->write();
            break;
        default:
            rs_error("InversionAcoustic2D<T>::saveHessian(): Unknown parameterisation."); 
            break;
    }
    if(Modmutefile.empty()){
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
    int i,j;
    int N,Ns,Nsrc;
    float *g_in;
    T *gvp, *grho, *gwav;
    std::string vpgradfile;
    std::string rhogradfile;
    std::string srcgradfile;
    if(Modmutefile.empty()){
        vpgradfile = VPGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    if(Srcmutefile.empty()){
        srcgradfile = SOURCEGRADFILE;
    }else{
        srcgradfile = SOURCEGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelAcoustic2D<T>> modelgrad (new rockseis::ModelAcoustic2D<T>(vpgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcegrad (new rockseis::Data2D<T>(srcgradfile));
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
            Nsrc = sourcegrad->getNtrace();
            sourcegrad->read();
            gwav = sourcegrad->getData();
            // Zero out gradient
            for(i=0; i< Ns; i++)
            {
                g[2*N+i] = 0;

            }
            // Stack src gradient
            for(j=0; j< Nsrc; j++)
            {
                for(i=0; i< Ns; i++)
                {
                    g[2*N+i] += gwav[Ns*j + i]*ksource;

                }
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
            Nsrc = sourcegrad->getNtrace();
            sourcegrad->read();
            gwav = sourcegrad->getData();
            // Zero out gradient
            for(i=0; i< Ns; i++)
            {
                g[2*N+i] = 0;

            }
            // Stack src gradient
            for(j=0; j< Nsrc; j++)
            {
                for(i=0; i< Ns; i++)
                {
                    g[2*N+i] += gwav[Ns*j + i]*ksource;

                }
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
void InversionAcoustic2D<T>::applySrcilum()
{
    if(this->getSrcilum()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADFILE, RHOGRADFILE, 1 ,0);
        // Ilumination
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> ilum (new rockseis::ModelAcoustic2D<T>(Srcilumfile, Srcilumfile, 1 ,0));

        // Mute model and write
        model->readModel();
        ilum->readModel();
        T *vp, *rho, *ilumodel;
        vp = model->getVp(); 
        rho = model->getR(); 
        ilumodel = ilum->getVp(); 
        int i;
        int N;
        T norm = 0.0;

        N = (model->getGeom())->getNtot();
        // Find norm of ilumination map
        for(i=0; i<N; i++) {
			if(((T) fabs(ilumodel[i])) >= norm) {
				norm = (T) fabs(ilumodel[i]);
			}
		}

        for(i=0; i< N; i++)
        {
            if(ilumodel[i] != 0.0){
                vp[i] = norm*vp[i]/ilumodel[i];
                rho[i] = norm*rho[i]/ilumodel[i];
            }
        }
        model->writeModel();
    }
}

template<typename T>
void InversionAcoustic2D<T>::applyMute()
{
    if(!Modmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute (new rockseis::ModelAcoustic2D<T>(Modmutefile, Modmutefile, 1 ,0));

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

    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> srcgrad (new rockseis::Data2D<T>(SOURCEGRADFILE));
        std::shared_ptr<rockseis::Data2D<T>> srcmute (new rockseis::Data2D<T>(Srcmutefile));
        T *srcgraddata, *srcmutedata;
        srcgrad->read();
        srcmute->read();
        srcgraddata = srcgrad->getData();
        srcmutedata = srcmute->getData();
        int Ns = srcgrad->getNt();
        int Nsrc = srcgrad->getNtrace();
        int i,j;
        for(j=0; j< Nsrc; j++)
        {
            for(i=0; i< Ns; i++)
            {
                srcgraddata[Ns*j+i] *= srcmutedata[i];
            }
        }
        srcgrad->setFile(SOURCEGRADMUTEFILE);
        srcgrad->write();
    }
}

template<typename T>
void InversionAcoustic2D<T>::computeTikhonovRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model (new rockseis::ModelAcoustic2D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readModel();
    T *vpgrad, *rhograd;
    double *c, *mod;
    int i;
    int N, Nmod;

    Nmod = (model->getGeom())->getNtot();
    model->readModel();
    model->setVpfile(VPREGGRADFILE);
    model->setRfile(RHOREGGRADFILE);
    vpgrad = model->getVp();
    rhograd = model->getR();
    T fvp = 0.0;
    T frho = 0.0;
    T fwrk;
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (model->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
                vpgrad[i] = x[i]; // * Covariance matrix
                vpgrad[i] *= kvp; // Diagonal scaling
                rhograd[i] = x[i+N]; // * Covariance matrix
                rhograd[i] *= krho; // Diagonal scaling

                fwrk = x[i]; // * Covariance matrix
                fwrk *= kvp;
                fvp += 0.5*(fwrk)*(fwrk); 

                fwrk = x[i+N]; //* Covariance matrix
                fwrk *= krho;
                fvp += 0.5*(fwrk)*(fwrk); 
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
            for(i=0; i< N; i++)
            {
                vpgrad[i] = mod[i]; // * Covariance matrix
                vpgrad[i] *= kvp; // Diagonal scaling

                fwrk = mod[i]; //* Covariance matrix
                fwrk *= kvp;
                fvp += 0.5*(fwrk)*(fwrk); 
            }
            for(i=0; i< N; i++)
            {
                c[i] = x[i+N];
            }
            spline->bisp();
            mod = spline->getMod();
            for(i=0; i< N; i++)
            {
                rhograd[i] = mod[i]; // * Covariance matrix
                rhograd[i] *= krho; // Diagonal scaling

                fwrk = mod[i]; //* Covariance matrix
                fwrk *= krho;
                fvp += 0.5*(fwrk)*(fwrk); 
            }

            break;
        default:
            rs_error("InversionAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
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

// =============== 3D ACOUSTIC INVERSION CLASS =============== //
//
template<typename T>
InversionAcoustic3D<T>::InversionAcoustic3D() {
    // Set default parameters
    apertx = -1;
    aperty = -1;

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
InversionAcoustic3D<T>::InversionAcoustic3D(MPImodeling *mpi): Inversion<T>(mpi) {
    // Set default parameters
    apertx = -1;
    aperty = -1;

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
InversionAcoustic3D<T>::~InversionAcoustic3D() {
    //Do nothing
}

template<typename T>
void InversionAcoustic3D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data3D<T>> shot3D;
    std::shared_ptr<rockseis::Data3D<T>> shot3Di;
    std::shared_ptr<rockseis::Data3D<T>> shotmod3D;
    std::shared_ptr<rockseis::Data3D<T>> shotmod3Di;
    std::shared_ptr<rockseis::Data3D<T>> shotres3D;
    std::shared_ptr<rockseis::Data3D<T>> shotres3Di;
    std::shared_ptr<rockseis::Data3D<T>> shotweight3D;
    std::shared_ptr<rockseis::Data3D<T>> shotweight3Di;
    std::shared_ptr<rockseis::Image3D<T>> vpgrad;
    std::shared_ptr<rockseis::Image3D<T>> rhograd;
    std::shared_ptr<rockseis::Image3D<T>> srcilum;
    std::shared_ptr<rockseis::Data3D<T>> wavgrad;
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic3D<T>> gmodel (new rockseis::ModelAcoustic3D<T>(Vpfile, Rhofile, this->getLpml() ,this->getFs()));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data3D<T>> source (new rockseis::Data3D<T>(Waveletfile));

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
        wavgrad = std::make_shared<rockseis::Data3D<T>>(1, source->getNt(), source->getDt(), 0.0);
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
        std::shared_ptr<rockseis::Data3D<T>> shot3D (new rockseis::Data3D<T>(Precordfile));
        // Create modelling and residual data files
        shotmod3D = std::make_shared<rockseis::Data3D<T>>(1, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
        shotmod3D->setFile(Pmodelledfile);
        shotmod3D->createEmpty(shot3D->getNtrace());

        shotres3D = std::make_shared<rockseis::Data3D<T>>(1, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
        shotres3D->setFile(Presidualfile);
        shotres3D->createEmpty(shot3D->getNtrace());
        
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

        // Perform work in parallel
        mpi->performWork();

        // Image
        vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile, gmodel, 1, 1, 1);
        vpgrad->createEmpty();

        rhograd = std::make_shared<rockseis::Image3D<T>>(Rhogradfile, gmodel, 1, 1, 1);
        rhograd->createEmpty();

        for(long int i=0; i<ngathers; i++) {
            vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            remove_file(Vpgradfile + "-" + std::to_string(i));
            rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
            remove_file(Rhogradfile + "-" + std::to_string(i));

            if(this->srcilumset){
                srcilum->stackImage(Srcilumfile + "-" + std::to_string(i));
                remove_file(Srcilumfile + "-" + std::to_string(i));
            }
        }

        //Clear work vector 
        mpi->clearWork();
	}
    else {
        /* Slave */
        std::shared_ptr<rockseis::FwiAcoustic3D<T>> fwi;
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
                shot3D = Sort->get3DGather(work.id);
                size_t ntr = shot3D->getNtrace();

                // Get the weight
                if(dataweight){
                    Sort->setDatafile(Dataweightfile);
                    shotweight3D = Sort->get3DGather(work.id);
                }

                lmodel = gmodel->getLocal(shot3D, apertx, aperty, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot3D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot3D, shot3Di);
                shot3Di->makeMap(lmodel->getGeom(), GMAP);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiAcoustic3D<T>>(lmodel, source, shot3Di, this->getOrder(), this->getSnapinc());

                // Create modelled and residual data objects 
                shotmod3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotmod3D->copyCoords(shot3D);
                shotmod3D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDatamodP(shotmod3D);
                shotres3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotres3D->copyCoords(shot3D);
                shotres3D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDataresP(shotres3D);

                // Interpolate weight
                if(dataweight){
                    shotweight3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(shotweight3D, shotweight3Di);
                    shotweight3Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweight(shotweight3Di);
                }
                
                // Setting misfit type
                fwi->setMisfit_type(this->getMisfit_type());
                fwi->setFilter(this->getFilter());
                fwi->setAllfreqs(this->getFreqs());

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1, 1);
                rhograd = std::make_shared<rockseis::Image3D<T>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, 1, 1,1);

//                if(this->srcilumset){
//                    srcilum = std::make_shared<rockseis::Image3D<T>>(Srcilumfile + "-" + std::to_string(work.id), lmodel, 1, 1,1);
//                    fwi->setSrcilum(srcilum);
//                }

                // Setting up gradient objects in fwi class
                fwi->setVpgrad(vpgrad);
                fwi->setRhograd(rhograd);

                wavgrad = std::make_shared<rockseis::Data3D<T>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
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

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fwi->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                shotmod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
                shotmod3Di->setFile(Pmodelledfile);
                interp->interp(shotmod3D, shotmod3Di);
                Sort->put3DGather(shotmod3Di, work.id);

                shotres3Di = std::make_shared<rockseis::Data3D<T>>(ntr, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
                shotres3Di->setFile(Presidualfile);
                interp->interp(shotres3D, shotres3Di);
                Sort->put3DGather(shotres3Di, work.id);

                // Output gradients
                if(this->getFilter()){
                    wavgrad->apply_filter(this->getFreqs());
                }
                wavgrad->putTrace(Wavgradfile, work.id);

                vpgrad->write();
                rhograd->write();

                // Output ilumination
                /*
                if(this->srcilumset){
                    srcilum->write();
                }
                */

                // Reset all classes
                shot3D.reset();
                shot3Di.reset();
                shotmod3D.reset();
                shotmod3Di.reset();
                shotres3D.reset();
                shotres3Di.reset();
                lmodel.reset();
                vpgrad.reset();
                rhograd.reset();
                wavgrad.reset();
                srcilum.reset();
                fwi.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void InversionAcoustic3D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	T rhosum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;
    std::string rhogradfile;

    if(Modmutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
        rhogradfile = RHOGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelAcoustic3D<T>> grad (new rockseis::ModelAcoustic3D<T>(vpgradfile, rhogradfile, this->getLpml() ,this->getFs()));

	// Read model
	grad->readModel();
	
	T *vpgrad, *rhograd;
	vpgrad = grad->getVp();
	rhograd = grad->getR();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl3D<T>> spline (new rockseis::Bspl3D<T>(grad->getNx(), grad->getNy(), grad->getNz(), grad->getDx(), grad->getDy(), grad->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3));
    int nc = spline->getNc();

	/* Allocating projection arrays */
	float *vpproj= (float *) calloc(nc, sizeof(float));
	if(vpproj==NULL){
		rs_error("InversionAcoustic3D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
	}
	float *rhoproj= (float *) calloc(nc, sizeof(float));
	if(rhoproj==NULL){
		rs_error("InversionAcoustic3D<T>::runBsproj(): Not enough memory to allocate projection array (rhoproj)");
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
			rs_error("InversionAcoustic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
				rhosum = 0.0;
                for(long int i=0; i<grad->getNx()*grad->getNy()*grad->getNz(); i++){
                    if(update_vp){
						vpsum += wrk[i]*vpgrad[i];
                    }
                    if(update_rho){
						rhosum += wrk[i]*rhograd[i];
                    }
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
			rs_error("InversionAcoustic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
int InversionAcoustic3D<T>::setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile)
{
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> model_in (new rockseis::ModelAcoustic3D<T>(vpfile, rhofile, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> source_in (new rockseis::Data3D<T>(sourcefile));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
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
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_rho) Npar += N;
            if(update_source) Npar += Ns;
            break;
        case PAR_BSPLINE:
            spline = std::make_shared<rockseis::Bspl3D<T>>(model_in->getNx(), model_in->getNy(), model_in->getNz(), model_in->getDx(), model_in->getDy(), model_in->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
            N=spline->getNc();
            Ns = source_in->getNt();
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_rho) Npar += N;
            if(update_source) Npar += Ns;
            break;
        default:
            rs_error("InversionAcoustic3D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void InversionAcoustic3D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> lsmodel (new rockseis::ModelAcoustic3D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> sourcels (new rockseis::Data3D<T>(SOURCELSFILE));
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
void InversionAcoustic3D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> model0 (new rockseis::ModelAcoustic3D<T>(VP0FILE, RHO0FILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> source0 (new rockseis::Data3D<T>(SOURCE0FILE));
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> lsmodel (new rockseis::ModelAcoustic3D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> lssource (new rockseis::Data3D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> mute;

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
    int N=0, Ns=0, Nmod=0, Npar=0;

    // If mute
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic3D<T>>(Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionAcoustic3D<T>::saveLinesearch(): Geometry in Modmutefile does not match geometry in the model.");
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
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i] + x[Npar+i]*rhomutedata[i]*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i];
                }
            }
            lsmodel->writeModel();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            Npar = 0;
            spline = std::make_shared<rockseis::Bspl3D<T>>(model0->getNx(), model0->getNy(), model0->getNz(), model0->getDx(), model0->getDy(), model0->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
            rs_error("InversionAcoustic3D<T>::saveLinesearch(): Unknown parameterisation."); 
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

    if(Modmutefile.empty()){
        free(vpmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionAcoustic3D<T>::readMisfit(double *f)
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
void InversionAcoustic3D<T>::readGrad(double *g)
{
    int i,j;
    int N,Ns,Nsrc,Npar=0;
    float *g_in;
    T *gvp, *grho, *gwav;
    std::string vpgradfile;
    std::string rhogradfile;
    if(Modmutefile.empty()){
        vpgradfile = VPGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelAcoustic3D<T>> modelgrad (new rockseis::ModelAcoustic3D<T>(vpgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> sourcegrad (new rockseis::Data3D<T>(SOURCEGRADFILE));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readModel();
            N = (modelgrad->getGeom())->getNtot();
            Npar = 0;
            gvp = modelgrad->getVp(); 
            grho = modelgrad->getR(); 
            if(update_vp){
                for(i=0; i < N; i++)
                {
                    g[i] = gvp[i]*kvp;
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
           spline = std::make_shared<rockseis::Bspl3D<T>>(modelgrad->getNx(), modelgrad->getNy(), modelgrad->getNz(), modelgrad->getDx(), modelgrad->getDy(), modelgrad->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
            rs_error("InversionAcoustic3D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
    if(update_source){
        Ns = sourcegrad->getNt();
        Nsrc = sourcegrad->getNtrace();
        sourcegrad->read();
        gwav = sourcegrad->getData();
        for(i=0; i< Ns; i++)
        {
            g[Npar+i] = 0.0;

        }
        for(j=0; j< Nsrc; j++)
        {
            for(i=0; i< Ns; i++)
            {
                g[Npar+i] += gwav[Ns*j + i]*ksource;

            }
        }
    }
}

template<typename T>
void InversionAcoustic3D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> grad;
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelAcoustic3D<T>>(VPGRADFILE, RHOGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelAcoustic3D<T>>(VPREGGRADFILE, RHOREGGRADFILE, 1 ,0);

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
        if(update_vp){
            vp[i] = vp[i] + reg_alpha[0]*regvp[i];
        }
        if(update_rho){
            rho[i] = rho[i] + reg_alpha[1]*regrho[i];
        }
    }
    grad->setVpfile(VPGRADCOMBFILE);
    grad->setRfile(RHOGRADCOMBFILE);
    grad->writeModel();
}

//template<typename T>
//void InversionAcoustic3D<T>::applySrcilum()
//{
//    if(this->getSrcilum()){
//        // Models
//        std::shared_ptr<rockseis::ModelAcoustic3D<T>> model;
//        model = std::make_shared<rockseis::ModelAcoustic3D<T>>(VPGRADFILE, RHOGRADFILE, 1 ,0);
//        // Ilumination
//        std::shared_ptr<rockseis::ModelAcoustic3D<T>> ilum (new rockseis::ModelAcoustic3D<T>(Srcilumfile, Srcilumfile, 1 ,0));
//
//        // Mute model and write
//        model->readModel();
//        ilum->readModel();
//        T *vp, *rho, *ilumodel;
//        vp = model->getVp(); 
//        rho = model->getR(); 
//        ilumodel = ilum->getVp(); 
//        int i;
//        int N;
//        T norm = 0.0;
//
//        N = (model->getGeom())->getNtot();
//        // Find norm of ilumination map
//        for(i=0; i<N; i++) {
//			if(((T) fabs(ilumodel[i])) >= norm) {
//				norm = (T) fabs(ilumodel[i]);
//			}
//		}
//
//        for(i=0; i< N; i++)
//        {
//            if(ilumodel[i] != 0.0){
//                vp[i] = norm*vp[i]/ilumodel[i];
//                rho[i] = norm*rho[i]/ilumodel[i];
//            }
//        }
//        model->writeModel();
//    }
//}

template<typename T>
void InversionAcoustic3D<T>::applyMute()
{
    if(!Modmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic3D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic3D<T>>(VPGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelAcoustic3D<T>> mute (new rockseis::ModelAcoustic3D<T>(Modmutefile, Modmutefile, 1 ,0));

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
void InversionAcoustic3D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> model (new rockseis::ModelAcoustic3D<T>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl3D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx,*dvpdy,*dvpdz;
    double *drhodx,*drhody,*drhodz;
    T *vpgrad, *rhograd;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Npar=0,Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdy = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    drhodx = (double *) calloc(Nmod, sizeof(double));
    drhody = (double *) calloc(Nmod, sizeof(double));
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
            Npar=0;
            if(update_vp){
                der->ddx_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdx[i] = df[i]*kvp;
                }
                der->ddy_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdy[i] = df[i]*kvp;
                }
                der->ddz_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdz[i] = df[i]*kvp;
                }
                Npar +=N;
            }
            if(update_rho){
                der->ddx_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    drhodx[i] = df[i]*krho;
                }
                der->ddy_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    drhody[i] = df[i]*krho;
                }
                der->ddz_fw(&x[Npar]);
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
           spline = std::make_shared<rockseis::Bspl3D<double>>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
                der->ddy_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvpdy[i] = df[i]*kvp;
                }
                der->ddz_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvpdz[i] = df[i]*kvp;
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
                der->ddy_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    drhody[i] = df[i]*krho;
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
            rs_error("InversionAcoustic3D<T>::computeRegularization(): Unknown parameterisation."); 
            break;
    }
    // Computing misfit
    double M; 
    T fvp = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i];
        M = sqrt(M);
        fvp += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
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
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdy[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdz[i]/M;
    }
    der->ddz_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] -= df[i];
    }

    // Computing misfit
    T frho = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i];
        M = sqrt(M);
        frho += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[1]*reg_eps[1];
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
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = drhody[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[1]*reg_eps[1];
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
    free(dvpdy);
    free(dvpdz);
    free(drhodx);
    free(drhody);
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
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uxdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2D;
    std::shared_ptr<rockseis::Data2D<T>> Uzdata2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uxdatamod2D;
    std::shared_ptr<rockseis::Data2D<T>> Uxdatamod2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uxdatares2D;
    std::shared_ptr<rockseis::Data2D<T>> Uxdatares2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uzdatamod2D;
    std::shared_ptr<rockseis::Data2D<T>> Uzdatamod2Di;
    std::shared_ptr<rockseis::Data2D<T>> Uzdatares2D;
    std::shared_ptr<rockseis::Data2D<T>> Uzdatares2Di;
    std::shared_ptr<rockseis::Data2D<T>> xweight2D;
    std::shared_ptr<rockseis::Data2D<T>> xweight2Di;
    std::shared_ptr<rockseis::Data2D<T>> zweight2D;
    std::shared_ptr<rockseis::Data2D<T>> zweight2Di;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::Image2D<T>> vsgrad;
    std::shared_ptr<rockseis::Image2D<T>> rhograd;
    std::shared_ptr<rockseis::Data2D<T>> wavgrad;
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
		std::shared_ptr<rockseis::Data2D<T>> Uxdata2D (new rockseis::Data2D<T>(Uxrecordfile));

		// Create a data class for the recorded data
		Uxdatamod2D = std::make_shared<rockseis::Data2D<T>>(1, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
		Uxdatamod2D->setFile(Uxmodelledfile);
		Uxdatamod2D->createEmpty(Uxdata2D->getNtrace());
		Uxdatares2D = std::make_shared<rockseis::Data2D<T>>(1, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
		Uxdatares2D->setFile(Uxresidualfile);
		Uxdatares2D->createEmpty(Uxdata2D->getNtrace());

		Uzdatamod2D = std::make_shared<rockseis::Data2D<T>>(1, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
		Uzdatamod2D->setFile(Uzmodelledfile);
		Uzdatamod2D->createEmpty(Uxdata2D->getNtrace());
		Uzdatares2D = std::make_shared<rockseis::Data2D<T>>(1, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
		Uzdatares2D->setFile(Uzresidualfile);
		Uzdatares2D->createEmpty(Uxdata2D->getNtrace());

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
                int gatherid;
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/2;
                }
                bool applyweightx = this->getDataweightx();
                bool applyweightz = this->getDataweightz();

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata2D = Sort->get2DGather(gatherid);
                size_t ntr = Uxdata2D->getNtrace();

                Sort->setDatafile(Uzrecordfile);
                Uzdata2D = Sort->get2DGather(gatherid);

                // Get the weight
                if(!Sort->getReciprocity()){
                    if(applyweightx){
                        Sort->setDatafile(Dataweightxfile);
                        xweight2D = Sort->get2DGather(gatherid);
                    }
                    if(applyweightz){
                        Sort->setDatafile(Dataweightzfile);
                        zweight2D = Sort->get2DGather(gatherid);
                    }
                }else{
                    if(work.id % 2 == 0){
                        if(applyweightx){
                            Sort->setDatafile(Dataweightxfile);
                            xweight2D = Sort->get2DGather(gatherid);
                        }
                        applyweightz = true;
                        zweight2D = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                        zweight2D->copyCoords(Uzdata2D);
                    }else{
                        applyweightx = true;
                        xweight2D = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                        xweight2D->copyCoords(Uxdata2D);
                        if(applyweightz){
                            Sort->setDatafile(Dataweightzfile);
                            zweight2D = Sort->get2DGather(gatherid);
                        }
                    }
                }

                // Get local model
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

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic2D<T>>(lmodel, source, Uxdata2Di, Uzdata2Di, this->getOrder(), this->getSnapinc());

                // Create modelled and residual data objects 
                Uxdatamod2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatamod2D->copyCoords(Uxdata2D);
                Uxdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Uxdatares2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatares2D->copyCoords(Uxdata2D);
                Uxdatares2D->makeMap(lmodel->getGeom(), GMAP);

                Uzdatamod2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatamod2D->copyCoords(Uzdata2D);
                Uzdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Uzdatares2D = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatares2D->copyCoords(Uzdata2D);
                Uzdatares2D->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdatamod2D->setField(VX);
                    Uxdatares2D->setField(VX);
                    Uzdatamod2D->setField(VZ);
                    Uzdatares2D->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdatamod2D->setField(PRESSURE);
                            Uxdatares2D->setField(PRESSURE);
                            Uzdatamod2D->setField(PRESSURE);
                            Uzdatares2D->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdatamod2D->setField(VX);
                            Uxdatares2D->setField(VX);
                            Uzdatamod2D->setField(VX);
                            Uzdatares2D->setField(VX);
                            break;
                        case 3:
                            Uxdatamod2D->setField(VZ);
                            Uxdatares2D->setField(VZ);
                            Uzdatamod2D->setField(VZ);
                            Uzdatares2D->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }

                fwi->setDatamodUx(Uxdatamod2D);
                fwi->setDataresUx(Uxdatares2D);
                fwi->setDatamodUz(Uzdatamod2D);
                fwi->setDataresUz(Uzdatares2D);

                // Interpolate weight
                if(applyweightx){
                    xweight2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(xweight2D, xweight2Di);
                    xweight2Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweightx(xweight2Di);
                }
                if(applyweightz){
                    zweight2Di = std::make_shared<rockseis::Data2D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(zweight2D, zweight2Di);
                    zweight2Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweightz(zweight2Di);
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
                    if(!Sort->getReciprocity()){
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
                    }else{
                        if(work.id % 2 == 0){
                            wavgrad->setField(VX);
                        }else{
                            wavgrad->setField(VZ);
                        }
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

                // Set reverse flag (For forward modelling only)
                fwi->setNoreverse(this->getNoreverse());

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


                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fwi->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                if(!Sort->getReciprocity()){
                    Uxdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                    Uxdatamod2Di->setFile(Uxmodelledfile);
                    interp->interp(Uxdatamod2D, Uxdatamod2Di);
                    Sort->put2DGather(Uxdatamod2Di, gatherid);

                    Uxdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                    Uxdatares2Di->setFile(Uxresidualfile);
                    interp->interp(Uxdatares2D, Uxdatares2Di);
                    Sort->put2DGather(Uxdatares2Di, gatherid);

                    Uzdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                    Uzdatamod2Di->setFile(Uzmodelledfile);
                    interp->interp(Uzdatamod2D, Uzdatamod2Di);
                    Sort->put2DGather(Uzdatamod2Di, gatherid);

                    Uzdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                    Uzdatares2Di->setFile(Uzresidualfile);
                    interp->interp(Uzdatares2D, Uzdatares2Di);
                    Sort->put2DGather(Uzdatares2Di, gatherid);
                }else{
                    switch(work.id % 2){
                        case 0:
                            Uxdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                            Uxdatamod2Di->setFile(Uxmodelledfile);
                            interp->interp(Uxdatamod2D, Uxdatamod2Di);
                            Sort->put2DGather(Uxdatamod2Di, gatherid);

                            Uxdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                            Uxdatares2Di->setFile(Uxresidualfile);
                            interp->interp(Uxdatares2D, Uxdatares2Di);
                            Sort->put2DGather(Uxdatares2Di, gatherid);
                            break;
                        case 1:
                            Uzdatamod2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                            Uzdatamod2Di->setFile(Uzmodelledfile);
                            interp->interp(Uzdatamod2D, Uzdatamod2Di);
                            Sort->put2DGather(Uzdatamod2Di, gatherid);

                            Uzdatares2Di = std::make_shared<rockseis::Data2D<T>>(ntr, Uxdata2D->getNt(), Uxdata2D->getDt(), Uxdata2D->getOt());
                            Uzdatares2Di->setFile(Uzresidualfile);
                            interp->interp(Uzdatares2D, Uzdatares2Di);
                            Sort->put2DGather(Uzdatares2Di, gatherid);
                            break;
                        default:
                            break;
                    }
                }

                // Output gradients
                if(update_source){
                    wavgrad->putTrace(Wavgradfile, work.id);
                }

                if(update_vp){
                    vpgrad->write();
                }
                if(update_vs){
                    vsgrad->write();
                }
                if(update_rho){
                    rhograd->write();
                }

                // Reset all classes
                Uxdata2D.reset();
                Uxdata2Di.reset();
                Uxdatamod2D.reset();
                Uxdatamod2Di.reset();
                Uxdatares2D.reset();
                Uxdatares2Di.reset();
                Uzdata2D.reset();
                Uzdata2Di.reset();
                Uzdatamod2D.reset();
                Uzdatamod2Di.reset();
                Uzdatares2D.reset();
                Uzdatares2Di.reset();
                lmodel.reset();
                vpgrad.reset();
                vsgrad.reset();
                rhograd.reset();
                wavgrad.reset();
                fwi.reset();

                // Send result back
                work.status = WORK_FINISHED;
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

    if(Modmutefile.empty()){
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

    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> sourcemute (new rockseis::Data2D<T>(Srcmutefile));
        Ns = sourcemute->getNt();
        if(Ns != source0->getNt()){
            rs_error("InversionElastic2D<T>::saveLinesearch(): Geometry in Srcmutefile does not match geometry in the Source file.");
        }
    }

    // If mute
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelElastic2D<T>>(Modmutefile, Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionElastic2D<T>::saveLinesearch(): Geometry in Modmutefile does not match geometry in the model.");
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
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i] + x[Npar+i]*rhomutedata[i]*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i];
                }
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
    if(Modmutefile.empty()){
        free(vpmutedata);
        free(vsmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionElastic2D<T>::saveHessian(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic2D<T>> lsmodel (new rockseis::ModelElastic2D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> lssource (new rockseis::Data2D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelElastic2D<T>> mute;


    // Write linesearch model
    lsmodel->readModel();
    lssource->read();
    T *vpls, *vsls, *rhols, *wavls;
    T *c, *mod;
    T *vpmutedata;
    T *vsmutedata;
    T *rhomutedata;
    vpls = lsmodel->getVp(); 
    vsls = lsmodel->getVs(); 
    rhols = lsmodel->getR(); 
    wavls = lssource->getData();
    int i;
    int N=0, Ns=0, Nmod=0, Npar=0;


    // Set output files for the Hessian models
    lsmodel->setVpfile(VPHESSFILE);
    lsmodel->setVsfile(VSHESSFILE);
    lsmodel->setRfile(RHOHESSFILE);
    lssource->setFile(SOURCEHESSFILE);

    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> sourcemute (new rockseis::Data2D<T>(Srcmutefile));
        Ns = sourcemute->getNt();
        if(Ns != lssource->getNt()){
            rs_error("InversionElastic2D<T>::saveHessian(): Geometry in Srcmutefile does not match geometry in the Source file.");
        }
    }

    // If mute
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelElastic2D<T>>(Modmutefile, Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionElastic2D<T>::saveHessian(): Geometry in Modmutefile does not match geometry in the model.");
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
            if(update_vp){
                for(i=0; i< N; i++)
                {
                    vpls[i] = x[i]*vpmutedata[i]*kvp*kvp;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    vpls[i] = 0.0;
                }
            }
            if(update_vs){
                for(i=0; i< N; i++)
                {
                    vsls[i] = x[Npar+i]*vsmutedata[i]*kvs*kvs;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    vsls[i] = 0.0;
                }
            }
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    rhols[i] = x[Npar+i]*rhomutedata[i]*krho*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    rhols[i] = 0.0;
                }
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
            spline = std::make_shared<rockseis::Bspl2D<T>>(lsmodel->getNx(), lsmodel->getNz(), lsmodel->getDx(), lsmodel->getDz(), this->getDtx(), this->getDtz(), 3, 3);
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
                    vpls[i] = mod[i]*vpmutedata[i]*kvp*kvp;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    vpls[i] = 0.0;
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
                    vsls[i] = mod[i]*vsmutedata[i]*kvs*kvs;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    vsls[i] = 0.0;
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
                    rhols[i] = mod[i]*rhomutedata[i]*krho*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< Nmod; i++)
                {
                    rhols[i] = 0.0;
                }

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
            rs_error("InversionElastic2D<T>::saveHessian(): Unknown parameterisation."); 
            break;
    }

    // Source wavelet
    Ns = lssource->getNt();
    for(i=0; i< Ns; i++)
    {
        if(update_source){
            wavls[i] = x[Npar+i]*ksource*ksource;
        }else{
            wavls[i] = 0.0;
        }
    }
    lssource->write();

    // Free allocated arrays
    if(Modmutefile.empty()){
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
    int i,j;
    int N,Ns,Nsrc,Npar=0;
    float *g_in;
    T *gvp, *gvs, *grho, *gwav;
    std::string vpgradfile;
    std::string vsgradfile;
    std::string rhogradfile;
    std::string srcgradfile;
    if(Modmutefile.empty()){
        vpgradfile = VPGRADFILE;
        vsgradfile = VSGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    if(Srcmutefile.empty()){
        srcgradfile = SOURCEGRADFILE;
    }else{
        srcgradfile = SOURCEGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelElastic2D<T>> modelgrad (new rockseis::ModelElastic2D<T>(vpgradfile, vsgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<T>> sourcegrad (new rockseis::Data2D<T>(srcgradfile));
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
        Nsrc = sourcegrad->getNtrace();
        sourcegrad->read();
        gwav = sourcegrad->getData();
        for(i=0; i< Ns; i++)
        {
            g[Npar+i] = 0.0;

        }
        for(j=0; j< Nsrc; j++)
        {
            for(i=0; i< Ns; i++)
            {
                g[Npar+i] += gwav[Ns*j + i]*ksource;

            }
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
    if(!Modmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelElastic2D<T>> model;
        model = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADCOMBFILE, VSGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelElastic2D<T>> mute (new rockseis::ModelElastic2D<T>(Modmutefile, Modmutefile, Modmutefile, 1 ,0));

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

    if(!Srcmutefile.empty()){
        std::shared_ptr<rockseis::Data2D<T>> srcgrad (new rockseis::Data2D<T>(SOURCEGRADFILE));
        std::shared_ptr<rockseis::Data2D<T>> srcmute (new rockseis::Data2D<T>(Srcmutefile));
        T *srcgraddata, *srcmutedata;
        srcgrad->read();
        srcmute->read();
        srcgraddata = srcgrad->getData();
        srcmutedata = srcmute->getData();
        int Ns = srcgrad->getNt();
        int Nsrc = srcgrad->getNtrace();
        int i,j;
        for(j=0; j< Nsrc; j++)
        {
            for(i=0; i< Ns; i++)
            {
                srcgraddata[Ns*j+i] *= srcmutedata[i];
            }
        }
        srcgrad->setFile(SOURCEGRADMUTEFILE);
        srcgrad->write();
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

// =============== 3D ELASTIC INVERSION CLASS =============== //
//
template<typename T>
InversionElastic3D<T>::InversionElastic3D() {
    // Set default parameters
    apertx = -1;
    aperty = -1;
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
InversionElastic3D<T>::InversionElastic3D(MPImodeling *mpi): Inversion<T>(mpi) {
    // Set default parameters
    apertx = -1;
    aperty = -1;
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
InversionElastic3D<T>::~InversionElastic3D() {
    //Do nothing
}

template<typename T>
void InversionElastic3D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data3D<T>> Uxdata3D;
    std::shared_ptr<rockseis::Data3D<T>> Uxdata3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uydata3D;
    std::shared_ptr<rockseis::Data3D<T>> Uydata3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uzdata3D;
    std::shared_ptr<rockseis::Data3D<T>> Uzdata3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uxdatamod3D;
    std::shared_ptr<rockseis::Data3D<T>> Uxdatamod3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uxdatares3D;
    std::shared_ptr<rockseis::Data3D<T>> Uxdatares3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uydatamod3D;
    std::shared_ptr<rockseis::Data3D<T>> Uydatamod3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uydatares3D;
    std::shared_ptr<rockseis::Data3D<T>> Uydatares3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uzdatamod3D;
    std::shared_ptr<rockseis::Data3D<T>> Uzdatamod3Di;
    std::shared_ptr<rockseis::Data3D<T>> Uzdatares3D;
    std::shared_ptr<rockseis::Data3D<T>> Uzdatares3Di;
    std::shared_ptr<rockseis::Data3D<T>> xweight3D;
    std::shared_ptr<rockseis::Data3D<T>> xweight3Di;
    std::shared_ptr<rockseis::Data3D<T>> yweight3D;
    std::shared_ptr<rockseis::Data3D<T>> yweight3Di;
    std::shared_ptr<rockseis::Data3D<T>> zweight3D;
    std::shared_ptr<rockseis::Data3D<T>> zweight3Di;
    std::shared_ptr<rockseis::Image3D<T>> vpgrad;
    std::shared_ptr<rockseis::Image3D<T>> vsgrad;
    std::shared_ptr<rockseis::Image3D<T>> rhograd;
    std::shared_ptr<rockseis::Data3D<T>> wavgrad;
    std::shared_ptr<rockseis::ModelElastic3D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Uxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelElastic3D<T>> gmodel (new rockseis::ModelElastic3D<T>(Vpfile, Vsfile, Rhofile, this->getLpml() ,this->getFs()));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data3D<T>> source (new rockseis::Data3D<T>(Waveletfile));

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
            ngathers *= 3;
        }

		// Wavelet gradient
		wavgrad = std::make_shared<rockseis::Data3D<T>>(1, source->getNt(), source->getDt(), 0.0);
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
		std::shared_ptr<rockseis::Data3D<T>> Uxdata3D (new rockseis::Data3D<T>(Uxrecordfile));

		// Create a data class for the recorded data
		Uxdatamod3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uxdatamod3D->setFile(Uxmodelledfile);
		Uxdatamod3D->createEmpty(Uxdata3D->getNtrace());
		Uxdatares3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uxdatares3D->setFile(Uxresidualfile);
		Uxdatares3D->createEmpty(Uxdata3D->getNtrace());

		Uydatamod3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uydatamod3D->setFile(Uymodelledfile);
		Uydatamod3D->createEmpty(Uxdata3D->getNtrace());
		Uydatares3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uydatares3D->setFile(Uyresidualfile);
		Uydatares3D->createEmpty(Uxdata3D->getNtrace());

		Uzdatamod3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uzdatamod3D->setFile(Uzmodelledfile);
		Uzdatamod3D->createEmpty(Uxdata3D->getNtrace());
		Uzdatares3D = std::make_shared<rockseis::Data3D<T>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
		Uzdatares3D->setFile(Uzresidualfile);
		Uzdatares3D->createEmpty(Uxdata3D->getNtrace());

		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,0,0});
			mpi->addWork(work);
		}

        // Perform work in parallel
        mpi->performWork();

        // Images
        vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile, gmodel, 1, 1, 1);
        vpgrad->createEmpty();

        vsgrad = std::make_shared<rockseis::Image3D<T>>(Vsgradfile, gmodel, 1, 1, 1);
        vsgrad->createEmpty();

        rhograd = std::make_shared<rockseis::Image3D<T>>(Rhogradfile, gmodel, 1, 1, 1);
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
        std::shared_ptr<rockseis::FwiElastic3D<T>> fwi;
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

                size_t gatherid;
                if(!Sort->getReciprocity()){
                    gatherid = work.id;
                }else{
                    gatherid = work.id/3;
                }

                bool applyweightx = this->getDataweightx();
                bool applyweighty = this->getDataweighty();
                bool applyweightz = this->getDataweightz();

                // Get the shot
                Sort->setDatafile(Uxrecordfile);
                Uxdata3D = Sort->get3DGather(gatherid);
                size_t ntr = Uxdata3D->getNtrace();

                Sort->setDatafile(Uyrecordfile);
                Uydata3D = Sort->get3DGather(gatherid);

                Sort->setDatafile(Uzrecordfile);
                Uzdata3D = Sort->get3DGather(gatherid);

                // Get the weight
                if(!Sort->getReciprocity()){
                    if(applyweightx){
                        Sort->setDatafile(Dataweightxfile);
                        xweight3D = Sort->get3DGather(gatherid);
                    }
                    if(applyweighty){
                        Sort->setDatafile(Dataweightyfile);
                        yweight3D = Sort->get3DGather(gatherid);
                    }
                    if(applyweightz){
                        Sort->setDatafile(Dataweightzfile);
                        zweight3D = Sort->get3DGather(gatherid);
                    }
                }else{
                    switch(work.id % 3){
                        case 0:
                            if(applyweightx){
                                Sort->setDatafile(Dataweightxfile);
                                xweight3D = Sort->get3DGather(gatherid);
                            }
                            applyweighty = true;
                            yweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            yweight3D->copyCoords(Uydata3D);
                            applyweightz = true;
                            zweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            zweight3D->copyCoords(Uzdata3D);
                            break;
                        case 1:
                            if(applyweighty){
                                Sort->setDatafile(Dataweightyfile);
                                yweight3D = Sort->get3DGather(gatherid);
                            }
                            applyweightx = true;
                            xweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            xweight3D->copyCoords(Uxdata3D);
                            applyweightz = true;
                            zweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            zweight3D->copyCoords(Uzdata3D);
                            break;
                        case 2:
                            if(applyweightz){
                                Sort->setDatafile(Dataweightzfile);
                                zweight3D = Sort->get3DGather(gatherid);
                            }
                            applyweightx = true;
                            xweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            xweight3D->copyCoords(Uxdata3D);
                            applyweighty = true;
                            yweight3D = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            yweight3D->copyCoords(Uydata3D);
                            break;
                        default:
                            break;
                    }
                }

                // Get local model
                lmodel = gmodel->getLocal(Uxdata3D, apertx, aperty, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata3D);
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
                        case 2:
                            source->setField(VY);
                            break;
                        case 3:
                            source->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }else{
                    switch(work.id % 3){
                        case 0:
                            source->setField(VX);
                            break;
                        case 1:
                            source->setField(VY);
                            break;
                        case 2:
                            source->setField(VZ);
                            break;
                        default:
                            break;
                    }
                }

                // Interpolate shot
                Uxdata3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata3D, Uxdata3Di);
                Uxdata3Di->makeMap(lmodel->getGeom(), GMAP);

                Uydata3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uydata3D, Uydata3Di);
                Uydata3Di->makeMap(lmodel->getGeom(), GMAP);

                Uzdata3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata3D, Uzdata3Di);
                Uzdata3Di->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdata3Di->setField(VX);
                    Uydata3Di->setField(VY);
                    Uzdata3Di->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdata3Di->setField(PRESSURE);
                            Uydata3Di->setField(PRESSURE);
                            Uzdata3Di->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdata3Di->setField(VX);
                            Uydata3Di->setField(VX);
                            Uzdata3Di->setField(VX);
                            break;
                        case 2:
                            Uxdata3Di->setField(VX);
                            Uydata3Di->setField(VX);
                            Uzdata3Di->setField(VX);
                            break;
                        case 3:
                            Uxdata3Di->setField(VZ);
                            Uydata3Di->setField(VZ);
                            Uzdata3Di->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }


                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic3D<T>>(lmodel, source, Uxdata3Di, Uydata3Di, Uzdata3Di, this->getOrder(), this->getSnapinc());

                // Create modelled and residual data objects 
                Uxdatamod3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatamod3D->copyCoords(Uxdata3D);
                Uxdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uxdatares3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatares3D->copyCoords(Uxdata3D);
                Uxdatares3D->makeMap(lmodel->getGeom(), GMAP);

                Uydatamod3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uydatamod3D->copyCoords(Uydata3D);
                Uydatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uydatares3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uydatares3D->copyCoords(Uydata3D);
                Uydatares3D->makeMap(lmodel->getGeom(), GMAP);

                Uzdatamod3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatamod3D->copyCoords(Uzdata3D);
                Uzdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uzdatares3D = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatares3D->copyCoords(Uzdata3D);
                Uzdatares3D->makeMap(lmodel->getGeom(), GMAP);

                if(!Sort->getReciprocity()){
                    Uxdatamod3D->setField(VX);
                    Uxdatares3D->setField(VX);
                    Uydatamod3D->setField(VY);
                    Uydatares3D->setField(VY);
                    Uzdatamod3D->setField(VZ);
                    Uzdatares3D->setField(VZ);
                }else{
                    switch(this->getSourcetype()){
                        case 0:
                            Uxdatamod3D->setField(PRESSURE);
                            Uxdatares3D->setField(PRESSURE);
                            Uydatamod3D->setField(PRESSURE);
                            Uydatares3D->setField(PRESSURE);
                            Uzdatamod3D->setField(PRESSURE);
                            Uzdatares3D->setField(PRESSURE);
                            break;
                        case 1:
                            Uxdatamod3D->setField(VX);
                            Uxdatares3D->setField(VX);
                            Uydatamod3D->setField(VX);
                            Uydatares3D->setField(VX);
                            Uzdatamod3D->setField(VX);
                            Uzdatares3D->setField(VX);
                            break;
                        case 2:
                            Uxdatamod3D->setField(VY);
                            Uxdatares3D->setField(VY);
                            Uydatamod3D->setField(VY);
                            Uydatares3D->setField(VY);
                            Uzdatamod3D->setField(VY);
                            Uzdatares3D->setField(VY);
                            break;
                        case 3:
                            Uxdatamod3D->setField(VZ);
                            Uxdatares3D->setField(VZ);
                            Uydatamod3D->setField(VZ);
                            Uydatares3D->setField(VZ);
                            Uzdatamod3D->setField(VZ);
                            Uzdatares3D->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown source type: ", std::to_string(this->getSourcetype()));
                            break;
                    }
                }

                fwi->setDatamodUx(Uxdatamod3D);
                fwi->setDataresUx(Uxdatares3D);
                fwi->setDatamodUy(Uydatamod3D);
                fwi->setDataresUy(Uydatares3D);
                fwi->setDatamodUz(Uzdatamod3D);
                fwi->setDataresUz(Uzdatares3D);

                // Interpolate weight
                if(applyweightx){
                    xweight3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(xweight3D, xweight3Di);
                    xweight3Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweightx(xweight3Di);
                }

                if(applyweighty){
                    yweight3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(yweight3D, yweight3Di);
                    yweight3Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweighty(yweight3Di);
                }
                if(applyweightz){
                    zweight3Di = std::make_shared<rockseis::Data3D<T>>(ntr, source->getNt(), source->getDt(), 0.0);
                    interp->interp(zweight3D, zweight3Di);
                    zweight3Di->makeMap(lmodel->getGeom(), GMAP);
                    fwi->setDataweightz(zweight3Di);
                }

                // Setting misfit type
                fwi->setMisfit_type(this->getMisfit_type());

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1, 1);
                fwi->setVpgrad(vpgrad);
                if(update_vs || update_rho){
                    vsgrad = std::make_shared<rockseis::Image3D<T>>(Vsgradfile + "-" + std::to_string(work.id), lmodel, 1, 1, 1);
                    fwi->setVsgrad(vsgrad);
                }
                if(update_rho){
                    rhograd = std::make_shared<rockseis::Image3D<T>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, 1, 1, 1);
                    fwi->setRhograd(rhograd);
                }

                if(update_source){
                    wavgrad = std::make_shared<rockseis::Data3D<T>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
                    // Copy geometry
                    wavgrad->copyCoords(source);
                    wavgrad->makeMap(lmodel->getGeom(), SMAP);
                    if(!Sort->getReciprocity()){
                        switch(this->getSourcetype()){
                            case 0:
                                wavgrad->setField(PRESSURE);
                                break;
                            case 1:
                                wavgrad->setField(VX);
                                break;
                            case 2:
                                wavgrad->setField(VY);
                                break;
                            case 3:
                                wavgrad->setField(VZ);
                                break;
                            default:
                                rs_error("Unknown wavgrad type: ", std::to_string(this->getSourcetype()));
                                break;
                        }
                    }else{
                        switch(work.id % 3){
                            case 0:
                                wavgrad->setField(VX);
                                break;
                            case 1:
                                wavgrad->setField(VY);
                                break;
                            case 2:
                                wavgrad->setField(VZ);
                                break;
                            default:
                                break;
                        }
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

                // Set reverse flag (For forward modelling only)
                fwi->setNoreverse(this->getNoreverse());

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

                

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fwi->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                if(!Sort->getReciprocity()){
                    Uxdatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                    Uxdatamod3Di->setFile(Uxmodelledfile);
                    interp->interp(Uxdatamod3D, Uxdatamod3Di);
                    Sort->put3DGather(Uxdatamod3Di, gatherid);

                    Uxdatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                    Uxdatares3Di->setFile(Uxresidualfile);
                    interp->interp(Uxdatares3D, Uxdatares3Di);
                    Sort->put3DGather(Uxdatares3Di, gatherid);

                    Uydatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uydata3D->getNt(), Uydata3D->getDt(), Uydata3D->getOt());
                    Uydatamod3Di->setFile(Uymodelledfile);
                    interp->interp(Uydatamod3D, Uydatamod3Di);
                    Sort->put3DGather(Uydatamod3Di, gatherid);

                    Uydatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uydata3D->getNt(), Uydata3D->getDt(), Uydata3D->getOt());
                    Uydatares3Di->setFile(Uyresidualfile);
                    interp->interp(Uydatares3D, Uydatares3Di);
                    Sort->put3DGather(Uydatares3Di, gatherid);

                    Uzdatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                    Uzdatamod3Di->setFile(Uzmodelledfile);
                    interp->interp(Uzdatamod3D, Uzdatamod3Di);
                    Sort->put3DGather(Uzdatamod3Di, gatherid);

                    Uzdatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                    Uzdatares3Di->setFile(Uzresidualfile);
                    interp->interp(Uzdatares3D, Uzdatares3Di);
                    Sort->put3DGather(Uzdatares3Di, gatherid);
                }else{
                    switch(work.id % 3){
                        case 0:
                            Uxdatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            Uxdatamod3Di->setFile(Uxmodelledfile);
                            interp->interp(Uxdatamod3D, Uxdatamod3Di);
                            Sort->put3DGather(Uxdatamod3Di, gatherid);

                            Uxdatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            Uxdatares3Di->setFile(Uxresidualfile);
                            interp->interp(Uxdatares3D, Uxdatares3Di);
                            Sort->put3DGather(Uxdatares3Di, gatherid);
                            break;
                        case 1:
                            Uydatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uydata3D->getNt(), Uydata3D->getDt(), Uydata3D->getOt());
                            Uydatamod3Di->setFile(Uymodelledfile);
                            interp->interp(Uydatamod3D, Uydatamod3Di);
                            Sort->put3DGather(Uydatamod3Di, gatherid);

                            Uydatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uydata3D->getNt(), Uydata3D->getDt(), Uydata3D->getOt());
                            Uydatares3Di->setFile(Uyresidualfile);
                            interp->interp(Uydatares3D, Uydatares3Di);
                            Sort->put3DGather(Uydatares3Di, gatherid);
                            break;
                        case 2:
                            Uzdatamod3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            Uzdatamod3Di->setFile(Uzmodelledfile);
                            interp->interp(Uzdatamod3D, Uzdatamod3Di);
                            Sort->put3DGather(Uzdatamod3Di, gatherid);

                            Uzdatares3Di = std::make_shared<rockseis::Data3D<T>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                            Uzdatares3Di->setFile(Uzresidualfile);
                            interp->interp(Uzdatares3D, Uzdatares3Di);
                            Sort->put3DGather(Uzdatares3Di, gatherid);
                            break;
                        default:
                            break;
                    }
                }

                // Output gradients
                if(update_source){
                    wavgrad->putTrace(Wavgradfile, work.id);
                }

                if(update_vp){
                    vpgrad->write();
                }
                if(update_vs){
                    vsgrad->write();
                }
                if(update_rho){
                    rhograd->write();
                }

                // Reset all classes
                Uxdata3D.reset();
                Uxdata3Di.reset();
                Uxdatamod3D.reset();
                Uxdatamod3Di.reset();
                Uxdatares3D.reset();
                Uxdatares3Di.reset();
                Uydata3D.reset();
                Uydata3Di.reset();
                Uydatamod3D.reset();
                Uydatamod3Di.reset();
                Uydatares3D.reset();
                Uydatares3Di.reset();
                Uzdata3D.reset();
                Uzdata3Di.reset();
                Uzdatamod3D.reset();
                Uzdatamod3Di.reset();
                Uzdatares3D.reset();
                Uzdatares3Di.reset();
                lmodel.reset();
                vpgrad.reset();
                vsgrad.reset();
                rhograd.reset();
                wavgrad.reset();
                fwi.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void InversionElastic3D<T>::runBsproj() {
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

    if(Modmutefile.empty()){
        vpgradfile = VPGRADCOMBFILE;
        vsgradfile = VSGRADCOMBFILE;
        rhogradfile = RHOGRADCOMBFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

	// Get gradients
	std::shared_ptr<rockseis::ModelElastic3D<T>> grad (new rockseis::ModelElastic3D<T>(vpgradfile, vsgradfile, rhogradfile, this->getLpml() ,this->getFs()));
    
    //Read gradients
    grad->readModel();

	T *vpgrad, *vsgrad, *rhograd;
	vpgrad = grad->getVp();
	vsgrad = grad->getVs();
	rhograd = grad->getR();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl3D<T>> spline (new rockseis::Bspl3D<T>(grad->getNx(), grad->getNy(), grad->getNz(), grad->getDx(), grad->getDy(), grad->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3));
    int nc = spline->getNc();

    float *vpproj, *vsproj, *rhoproj;
    /* Allocating projection arrays */
    vpproj= (float *) calloc(nc, sizeof(float));
    if(vpproj==NULL){
        rs_error("InversionElastic3D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
    }
    vsproj= (float *) calloc(nc, sizeof(float));
    if(vsproj==NULL){
        rs_error("InversionElastic3D<T>::runBsproj(): Not enough memory to allocate projection array (vsproj)");
    }
    rhoproj= (float *) calloc(nc, sizeof(float));
    if(rhoproj==NULL){
        rs_error("InversionElastic3D<T>::runBsproj(): Not enough memory to allocate projection array (rhoproj)");
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
			rs_error("InversionElastic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
                rhosum = 0.0;
                for(long int i=0; i<grad->getNx()*grad->getNy()*grad->getNz(); i++){
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
            rs_error("InversionElastic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
int InversionElastic3D<T>::setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile)
{
    std::shared_ptr<rockseis::ModelElastic3D<T>> model_in (new rockseis::ModelElastic3D<T>(vpfile, vsfile, rhofile, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> source_in (new rockseis::Data3D<T>(sourcefile));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
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
             spline = std::make_shared<rockseis::Bspl3D<T>>(model_in->getNx(), model_in->getNy(), model_in->getNz(), model_in->getDx(), model_in->getDy(), model_in->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
            N=spline->getNc();
            Ns = source_in->getNt();
            Npar = 0;
            if(update_vp) Npar += N;
            if(update_vs) Npar += N;
            if(update_rho) Npar += N;
            if(update_source) Npar += Ns;
            break;
        default:
            rs_error("InversionElastic3D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void InversionElastic3D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelElastic3D<T>> lsmodel (new rockseis::ModelElastic3D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> sourcels (new rockseis::Data3D<T>(SOURCELSFILE));
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
void InversionElastic3D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic3D<T>> model0 (new rockseis::ModelElastic3D<T>(VP0FILE, VS0FILE, RHO0FILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> source0 (new rockseis::Data3D<T>(SOURCE0FILE));
    std::shared_ptr<rockseis::ModelElastic3D<T>> lsmodel (new rockseis::ModelElastic3D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> lssource (new rockseis::Data3D<T>(SOURCELSFILE));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
    std::shared_ptr<rockseis::ModelElastic3D<T>> mute;

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
    if(!Modmutefile.empty()){
        mute = std::make_shared <rockseis::ModelElastic3D<T>>(Modmutefile, Modmutefile, Modmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("InversionElastic3D<T>::saveLinesearch(): Geometry in Modmutefile does not match geometry in the model.");
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
            if(update_rho){
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i] + x[Npar+i]*rhomutedata[i]*krho;
                }
                Npar += N;
            }else{
                for(i=0; i< N; i++)
                {
                    rhols[i] = rho0[i];
                }
            }
            lsmodel->writeModel();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            Npar = 0;
            spline = std::make_shared<rockseis::Bspl3D<T>>(model0->getNx(), model0->getNy(), model0->getNz(), model0->getDx(), model0->getDy(), model0->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
            rs_error("InversionAcoustic3D<T>::saveLinesearch(): Unknown parameterisation."); 
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
    if(Modmutefile.empty()){
        free(vpmutedata);
        free(vsmutedata);
        free(rhomutedata);
    }
}

template<typename T>
void InversionElastic3D<T>::readMisfit(double *f)
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
void InversionElastic3D<T>::readGrad(double *g)
{
    int i,j;
    int N,Ns,Nsrc,Npar=0;
    float *g_in;
    T *gvp, *gvs, *grho, *gwav;
    std::string vpgradfile;
    std::string vsgradfile;
    std::string rhogradfile;
    if(Modmutefile.empty()){
        vpgradfile = VPGRADFILE;
        vsgradfile = VSGRADFILE;
        rhogradfile = RHOGRADFILE;
    }else{
        vpgradfile = VPGRADMUTEFILE;
        vsgradfile = VSGRADMUTEFILE;
        rhogradfile = RHOGRADMUTEFILE;
    }

    std::shared_ptr<rockseis::ModelElastic3D<T>> modelgrad (new rockseis::ModelElastic3D<T>(vpgradfile, vsgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data3D<T>> sourcegrad (new rockseis::Data3D<T>(SOURCEGRADFILE));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
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
           spline = std::make_shared<rockseis::Bspl3D<T>>(modelgrad->getNx(), modelgrad->getNy(), modelgrad->getNz(), modelgrad->getDx(), modelgrad->getDy(), modelgrad->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
            rs_error("InversionElastic3D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
    if(update_source){
        Ns = sourcegrad->getNt();
        Nsrc = sourcegrad->getNtrace();
        sourcegrad->read();
        gwav = sourcegrad->getData();
        for(i=0; i< Ns; i++)
        {
            g[Npar+i] = 0.0;

        }
        for(j=0; j< Nsrc; j++)
        {
            for(i=0; i< Ns; i++)
            {
                g[Npar+i] += gwav[Ns*j + i]*ksource;

            }
        }

    }
}

template<typename T>
void InversionElastic3D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelElastic3D<T>> grad;
    std::shared_ptr<rockseis::ModelElastic3D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelElastic3D<T>>(VPGRADFILE, VSGRADFILE, RHOGRADFILE, 1 ,0);
    reggrad = std::make_shared<rockseis::ModelElastic3D<T>>(VPREGGRADFILE, VSREGGRADFILE, RHOREGGRADFILE, 1 ,0);

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
void InversionElastic3D<T>::applyMute()
{
    if(!Modmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelElastic3D<T>> model;
        model = std::make_shared<rockseis::ModelElastic3D<T>>(VPGRADCOMBFILE, VSGRADCOMBFILE, RHOGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelElastic3D<T>> mute (new rockseis::ModelElastic3D<T>(Modmutefile, Modmutefile, Modmutefile, 1 ,0));

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

    if(!Srcmutefile.empty()){
    }
}

template<typename T>
void InversionElastic3D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelElastic3D<T>> model (new rockseis::ModelElastic3D<T>(VPLSFILE, VSLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl3D<double>> spline;

    // Write linesearch model
    model->readModel();
    double *dvpdx,*dvpdy,*dvpdz;
    double *dvsdx,*dvsdy,*dvsdz;
    double *drhodx,*drhody,*drhodz;
    T *vpgrad, *vsgrad, *rhograd;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Npar=0,Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdy = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    dvsdx = (double *) calloc(Nmod, sizeof(double));
    dvsdy = (double *) calloc(Nmod, sizeof(double));
    dvsdz = (double *) calloc(Nmod, sizeof(double));
    drhodx = (double *) calloc(Nmod, sizeof(double));
    drhody = (double *) calloc(Nmod, sizeof(double));
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
                der->ddy_fw(x);
                for(i=0; i< N; i++)
                {
                    dvpdy[i] = df[i]*kvp;
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
                der->ddy_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    dvsdy[i] = df[i]*kvs;
                }
                der->ddz_fw(&x[Npar]);
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
                der->ddy_fw(&x[Npar]);
                for(i=0; i< N; i++)
                {
                    drhody[i] = df[i]*krho;
                }
                der->ddz_fw(&x[Npar]);
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
           spline = std::make_shared<rockseis::Bspl3D<double>>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), this->getDtx(), this->getDty(), this->getDtz(), 3, 3, 3);
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
                der->ddy_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvpdy[i] = df[i]*kvp;
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
                der->ddy_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    dvsdy[i] = df[i]*kvs;
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
                der->ddy_fw(mod);
                for(i=0; i< Nmod; i++)
                {
                    drhody[i] = df[i]*krho;
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
            rs_error("InversionElastic3D<T>::computeRegularization(): Unknown parameterisation."); 
            break;
    }
    // Computing misfit
    double M; 
    T fvp = 0.0;
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i];
        M = sqrt(M);
        fvp += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
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
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdy[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i] + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
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
        M = dvsdx[i]*dvsdx[i] + dvsdy[i]*dvsdy[i] + dvsdz[i]*dvsdz[i];
        M = sqrt(M);
        fvs += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = dvsdx[i]*dvsdx[i] + dvsdy[i]*dvsdy[i] + dvsdz[i]*dvsdz[i] + reg_eps[1]*reg_eps[1];
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
        M = dvsdx[i]*dvsdx[i] + dvsdy[i]*dvsdy[i] + dvsdz[i]*dvsdz[i] + reg_eps[1]*reg_eps[1];
        M = sqrt(M);
        gwrk[i] = dvsdy[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vsgrad[i] -= df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvsdx[i]*dvsdx[i] + dvsdy[i]*dvsdy[i] + dvsdz[i]*dvsdz[i] + reg_eps[1]*reg_eps[1];
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
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i];
        M = sqrt(M);
        frho += M;
    }

    // Computing gradient
    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[2]*reg_eps[2];
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
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[2]*reg_eps[2];
        M = sqrt(M);
        gwrk[i] = drhody[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        rhograd[i] = -1.0*df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = drhodx[i]*drhodx[i] + drhody[i]*drhody[i] + drhodz[i]*drhodz[i] + reg_eps[2]*reg_eps[2];
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
    free(dvpdy);
    free(dvpdz);
    free(dvsdx);
    free(dvsdy);
    free(dvsdz);
    free(drhodx);
    free(drhody);
    free(drhodz);
    free(gwrk);
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Inversion<float>;
template class Inversion<double>;

template class InversionAcoustic2D<float>;
template class InversionAcoustic2D<double>;

template class InversionAcoustic3D<float>;
template class InversionAcoustic3D<double>;

template class InversionElastic2D<float>;
template class InversionElastic2D<double>;

template class InversionElastic3D<float>;
template class InversionElastic3D<double>;

}


