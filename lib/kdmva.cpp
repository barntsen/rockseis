#include "kdmva.h"


namespace rockseis {

// =============== ABSTRACT KDMVA CLASS =============== //
template<typename T>
Kdmva<T>::Kdmva() {
    //Set default parameters
    fs = false;
    lpml = 0;
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
    if(createLog() == KVA_ERR)
    {
        rs_error("Kdmva<T>::Kdmva(): Error creating logfile for writting.");
    }
    if(createProglog() == KVA_ERR)
    {
        rs_error("Kdmva<T>::Kdmva(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
Kdmva<T>::Kdmva(MPImodeling *_mpi) {
    mpi = _mpi;

    //Set default parameters
    fs = false;
    lpml = 0;
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
    if(createLog() == KVA_ERR)
    {
        rs_error("Kdmva<T>::Kdmva(): Error creating logfile for writting.");
    }
    if(createProglog() == KVA_ERR)
    {
        rs_error("Kdmva<T>::Kdmva(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
void Kdmva<T>::normalize(double *v, double *f, int n){
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
double Kdmva<T>::vector_norm(double *v, const int type, const int n){
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
void Kdmva<T>::find_max(T *v, T* max, int *imax, int n){
	// Variables
	int i;

    *max = v[0];
    *imax = 0;
    for(i=1; i < n; i++)
    {
        if (v[i] > *max){
            *max = v[i];
            *imax = i;
        }
    }
}

template<typename T>
bool Kdmva<T>::createLog(){
	logfile = LOGFILE;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return KVA_ERR;
	}else{
		Flog.close();
		return KVA_OK;
	}
}

template<typename T>
bool Kdmva<T>::createProglog(){
	progresslogfile = PROGLOGFILE;
	Flog.open(progresslogfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return KVA_ERR;
	}else{
		Flog.close();
		return KVA_OK;
	}
}

template<typename T>
void Kdmva<T>::writeLog(std::string text){
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
void Kdmva<T>::writeProgress(std::string text){
    if(!progresslogfile.empty()){
        Flog.open(progresslogfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Kdmva<T>::createResult(){
    struct stat s;
    // Checking if result folder is present, and creates it if not
    if(stat(RESULTDIR,&s) != 0) {
        int mkdir_return = mkdir(RESULTDIR,0777);
        if(mkdir_return != 0) rs_error("Not able to create result directory: ", RESULTDIR);
    }
}

template<typename T>
Kdmva<T>::~Kdmva() {
    //Do nothing
}

// =============== 2D ACOUSTIC KDMVA CLASS =============== //
//
template<typename T>
KdmvaAcoustic2D<T>::KdmvaAcoustic2D() {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;

    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
KdmvaAcoustic2D<T>::KdmvaAcoustic2D(MPImodeling *mpi): Kdmva<T>(mpi) {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;
    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
KdmvaAcoustic2D<T>::~KdmvaAcoustic2D() {
    //Do nothing
}

template<typename T>
void KdmvaAcoustic2D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> source;
    std::shared_ptr<rockseis::Data2D<T>> shot2D;
    std::shared_ptr<rockseis::RaysAcoustic2D<T>> rays;
    std::shared_ptr<rockseis::Ttable<T>> Ttable;
    std::shared_ptr<rockseis::Image2D<T>> pimage;
    std::shared_ptr<rockseis::Image2D<T>> lpimage;
	std::shared_ptr<rockseis::ModelEikonal2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> gmodel (new rockseis::ModelAcoustic2D<T>(Vpfile, Rhofile, this->getLpml(), this->getFs()));

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());

    // Test for problematic model sampling
    if(gmodel->getDx() != gmodel->getDz()){
        rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dz).");
    }

    // Read and expand global model
    gmodel->readVelocity();
    gmodel->Expand();

	if(mpi->getRank() == 0) {
		// Master

    // -------------------------------------Calculate traveltimes
        // Get number of receivers
        Sort->createReceivermap(Surveyfile); 
        size_t nrecgath =  Sort->getNensemb();

        // Get number of shots
        Sort->createShotmap(Surveyfile); 
        size_t nsougath =  Sort->getNensemb();

        Sort->writeKeymap();
        Sort->writeSortmap();

        nsoufin = nsougath/souinc + 1;
        if(nsoufin > nsougath) nsoufin = nsougath;

        nrecfin = nrecgath/recinc + 1;
        if(nrecfin > nrecgath) nrecfin = nrecgath;

        // -------------------------------------Create a travel time table class
        size_t ngathers = nsoufin + nrecfin;
        Ttable = std::make_shared<rockseis::Ttable<T>> (gmodel, ngathers);
        Ttable->setFilename(Ttablefile);
        Ttable->createEmptyttable();

        /******************      Creating source side traveltime   ***********************/
		// Create work queue
		for(long int i=0; i<nsoufin; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
			mpi.addWork(work);
		}

        // Broadcast nsoufin
        MPI_Bcast(&nsoufin, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Perform work in parallel
		mpi.performWork();

        // Reset mpi 
        mpi.clearWork();

        /******************       Creating receiver side traveltime    ***********************/
        // Create new list of positions
        Sort->createReceivermap(Precordfile); 
        Sort->setReciprocity(true);
        Sort->writeKeymap();
        Sort->writeSortmap();

        for(long int i=0; i<nrecfin; i++) {
            // Work struct
            std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
            mpi.addWork(work);
        }

        // Perform work in parallel
        mpi.performWork();

		//Clear work vector 
		mpi->clearWork();

    // -------------------------------------Migrate data
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort->getNensemb();

        // Image
        pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile, gmodel, this->getNhx(), this->getNhz());
        pimage->createEmpty();

        // Create work queue
        for(long int i=0; i<ngathers; i++) {
            // Work struct
            std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
            mpi.addWork(work);
        }

        // Perform work in parallel
        mpi.performWork();

        //Calculate and output misfit and residual
        this->computeMisfit(pimage);

		//Clear work vector 
		mpi->clearWork();

// -------------------------------------Compute gradient
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Gradient
        vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile, gmodel, 1, 1);
        vpgrad->createEmpty();

        // Create work queue
        for(long int i=0; i<ngathers; i++) {
            // Work struct
            std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
            mpi.addWork(work);
        }

        // Perform work in parallel
        mpi.performWork();

		//Clear work vector 
		mpi->clearWork();
    }
    else {
        /* Slave */
        size_t number;
        std::shared_ptr<rockseis::Data2D<T>> Shotgeom;

        // Receive nsoufin from master
        MPI_Bcast(&nsoufin, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /******************             Creating source side     ***********************/
        while(1) {
            workModeling_t work = mpi.receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi.sendNoWork(0);
            }
            else {
                // Calculate traveltime
                Sort->readKeymap();
                Sort->readSortmap();

                number = work.id*souinc;
                if (number > Sort->getNensemb()-1) number = Sort->getNensemb()-1;
                Shotgeom = Sort->get2DGather(number);

                // Set shot coordinates and make a map
	            source = std::make_shared<rockseis::Data2D<T>>(1, 1, 1.0, 0.0);
                source->copyCoords(Shotgeom);

                source->makeMap(gmodel->getGeom(), SMAP);

                // Run modelling 
                rays = std::make_shared<rockseis::RaysAcoustic2D<T>>(gmodel);

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();

                // Create traveltime table
                Ttable = std::make_shared<rockseis::Ttable<T>> (Ttablefile);
                Ttable->fetchTtabledata(rays, source, work.id); //Get traveltime data
                Ttable->writeTtable(work.id);
                Ttable.reset();
                
                // Reset all classes
                source.reset();
                Shotgeom.reset();
                rays.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
        /******************             Creating receiver side     ***********************/
        while(1) {
            workModeling_t work = mpi.receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi.sendNoWork(0);
            }
            else {
                // Do some work
                Sort->readKeymap();
                Sort->readSortmap();

                number = work.id*recinc;
                if (number > Sort->getNensemb()-1) number = Sort->getNensemb()-1;
                Shotgeom = Sort->get2DGather(number);

                // Set shot coordinates and make a map
	            source = std::make_shared<rockseis::Data2D<T>>(1, 1, 1.0, 0.0);
                source->copyCoords(Shotgeom);

                source->makeMap(gmodel->getGeom(), SMAP);

                // Run modelling 
                rays = std::make_shared<rockseis::RaysAcoustic2D<T>>(gmodel);

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();

                // Create traveltime table
                Ttable = std::make_shared<rockseis::Ttable<T>> (Ttablefile);
                Ttable->fetchTtabledata(rays, source, work.id+nsoufin); //Get traveltime data
                Ttable->writeTtable(work.id+nsoufin);
                Ttable.reset();
                
                // Reset all classes
                source.reset();
                Shotgeom.reset();
                rays.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		

                // ------------------------------------Migrate data
            }
        }
        std::shared_ptr<rockseis::KdmigAcoustic2D<T>> kdmig;
        while(1) {
            workModeling_t work = mpi.receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi.sendNoWork(0);
            }
            else {
                // Do migration
                Sort->readKeymap();
                Sort->readSortmap();

                // Get the shot
                shot2D = Sort->get2DGather(work.id);

                // Make local model
                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);
                lmodel->Expand();

                // Make image class
                pimage = std::make_shared<rockseis::Image2D<T>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhz);

                // Create traveltime table class
                ttable = std::make_shared<rockseis::Ttable<T>>(Ttablefile);
                ttable->allocTtable();

                // Create imaging class
                kdmig = std::make_shared<rockseis::KdmigAcoustic2D<T>>(lmodel, ttable, shot2D, pimage);

                // Set frequency decimation 
                kdmig->setFreqinc(1);

                // Set minimum and maximum frequency to migrate
                kdmig->setMinfreq(0.0);
                kdmig->setMaxfreq(125.0);

                // Set radius of interpolation
                kdmig->setRadius(radius);

                // Set logfile
                kdmig->setLogfile("log.txt-" + std::to_string(work.id));

                // Run migration
                kdmig->run();

                // Send result back
                work.status = PARALLEL_IO;
                mpi.sendResult(work);		

                // Stack image
                pimage->stackImage_parallel(Pimagefile);

                // Reset all classes
                shot2D.reset();
                lmodel.reset();
                pimage.reset();
                kdmig.reset();
                ttable.reset();
                pimage.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi.sendResult(work);		
            }
        }

        // Compute Adjoint
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
                Sort->readKeymap();
                Sort->readSortmap();

                // Get the shot
                shot2D = Sort->get2DGather(work.id);

                // Make local model
                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);
                lmodel->Expand();

                // Make image class
                pimage = std::make_shared<rockseis::Image2D<T>>(PIMAGERESFILE);
                lpimage = pimage->getLocal(shot2D, apertx, SMAP);

                // Create traveltime table class
                ttable = std::make_shared<rockseis::Ttable<T>>(Ttablefile);
                ttable->allocTtable();

                // Create imaging class
                kdmig = std::make_shared<rockseis::KdmigAcoustic2D<T>>(lmodel, ttable, shot2D, lpimage);

                // Set frequency decimation 
                kdmig->setFreqinc(1);

                // Set minimum and maximum frequency to migrate
                kdmig->setMinfreq(0.0);
                kdmig->setMaxfreq(125.0);

                // Set radius of interpolation
                kdmig->setRadius(radius);

                // Creating gradient object
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);

                // Setting up gradient objects in wemvafwi class
                mva->setVpgrad(vpgrad);

                // Set logfile
                kdmig->setLogfile("log.txt-" + std::to_string(work.id));

                // Run migration
                kdmig->run_adj();

                // Send result back
                work.status = PARALLEL_IO;
                mpi.sendResult(work);		

                // Stack image
                vpgrad->stackImage_parallel(Vpgradfile);

                // Reset all classes
                shot2D.reset();
                lmodel.reset();
                pimage.reset();
                kdmig.reset();
                ttable.reset();
                vpgrad.reset();
                lpimage.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi.sendResult(work);		
            }
        }
    }

}

template<typename T>
void KdmvaAcoustic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;

    if(Modelmutefile.empty()){
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
		rs_error("KdmvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
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
			rs_error("KdmvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
			rs_error("KdmvaAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
        MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

	   }
    // Free variables
    free(global_stack);
    free(vpproj);
}

template<typename T>
int KdmvaAcoustic2D<T>::setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile)
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
            rs_error("KdmvaAcoustic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return N;
}

template<typename T>
void KdmvaAcoustic2D<T>::saveResults(int iter)
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
void KdmvaAcoustic2D<T>::saveLinesearch(double *x)
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
    if(!Modelmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Modelmutefile, Modelmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("KdmvaAcoustic2D<T>::saveLinesearch(): Geometry in Modelmutefile does not match geometry in the model.");
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
            rs_error("KdmvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Modelmutefile.empty()){
        free(vpmutedata);
    }
}

template<typename T>
void KdmvaAcoustic2D<T>::computeMisfit(std::shared_ptr<rockseis::Image2D<T>> pimage)
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
    T f3 = 0.;


    std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute;
    int N;
    T *mutedata;
    // If mute
    if(!Residualmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Residualmutefile, Residualmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (nx*nz);
        if(N != Nmute) rs_error("KdmvaAcoustic2D<T>::computeMisfit(): Geometry in Residualmutefile does not match geometry in the image.");
        mute->readModel();
        mutedata = mute->getVp();
    }else{
        N = (nx*nz);
        mutedata = (T *) calloc(N, sizeof(T)); 
        for(ix=0; ix < N; ix++){
            mutedata[ix] = 1.0;
        }
    }


    /* Variables for DS_hmax misfit and residual computation*/
    T *cip = (T *) calloc(nz*nhx*nhz, sizeof(T));
    T *env = (T *) calloc(nz*nhx*nhz, sizeof(T));
    T *hmax = (T *) calloc(nz, sizeof(T));
    T *hsort = (T *) calloc(nz, sizeof(T));
    T *hwrk = (T *) calloc(nhx, sizeof(T));
    int *imax = (int *) calloc(nz, sizeof(int));
    std::shared_ptr<Hilbert<T>> hilb_cip1 (new Hilbert<T>(nz, nhx, nhz, nz));
    std::shared_ptr<Hilbert<T>> hilb_cip2 (new Hilbert<T>(nz, nhx, nhz, nz));
    T pclip;
    T *imag;
    int pos;
    T num, den;
    T rms=0;
    size_t fold=0;
    

    switch(this->getMisfit_type()){
        case SI:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
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
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
                    }
                }
            }
            break;
        case DS_PLUS_SI:
            // Misfit
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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

            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
                    }
                }
            }
            break;
        case DS_HMAX:
            // Misfit
            // Get a CIP gather
            for (ix=0; ix<nx; ix++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        // Derivative
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                    }
                }
                // Hilbert transform over first non-singleton axis
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();

                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=0; iz<nz; iz++){
                           env[kres2D(iz,ihx,ihz)] = SQ(cip[kres2D(iz,ihx,ihz)]) + SQ(imag[kres2D(iz,ihx,ihz)]);
                        }
                    }
                }
                // Find maximum and index of maximum along hx axis and store in an array
                for (iz=0; iz<nz; iz++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            hwrk[ihx] = env[kres2D(iz,ihx,ihz)];
                        }
                    }
                    this->find_max(hwrk, &hmax[iz], &imax[iz], nhx);
                }

                // Threshold 
                for (iz=0; iz<nz; iz++){
                    hsort[iz] = hmax[iz];
                }
                std::sort(hsort, hsort+nz); 
                pos = (int) (THRES*nz/100);
                pclip = hsort[pos];
                for (iz=0; iz<nz; iz++){
                    if(hmax[iz] < pclip){
                        hmax[iz] = 0;
                        imax[iz] = (nhx-1)/2;
                    }
                }
                // Calculate delta h, and misfit. 
                for (iz=0; iz<nz; iz++){
                    f1 += 0.5*SQ(imax[iz]-((nhx-1)/2));
                    f2 += 0.5*SQ(cip[kres2D(iz,(nhx-1)/2,(nhz-1)/2)]);
                }

                for (iz=0; iz<nz; iz++){
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        rms += SQ((env[kres2D(iz,imax[iz]+2,0)] - 2.0*env[kres2D(iz,imax[iz],0)] + env[kres2D(iz,imax[iz]-2,0)])/4.0);
                        fold++;
                    }
                }
            }
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
                            f3 += 0.5*G2*wrk[iz]*wrk[iz];
                        }
                                            }
                }
            }

            f = f1*f3/f2;

            if(fold >0){
                rms = (T) std::sqrt(rms)/fold;
            }

            // Residual
            for (ix=0; ix<nx; ix++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        // Derivative
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                    }
                }
                // Hilbert transform over first non-singleton axis
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();

                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=0; iz<nz; iz++){
                           env[kres2D(iz,ihx,ihz)] = SQ(cip[kres2D(iz,ihx,ihz)]) + SQ(imag[kres2D(iz,ihx,ihz)]);
                        }
                    }
                }
                // Find maximum and index of maximum along hx axis and store in an array
                for (iz=0; iz<nz; iz++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            hwrk[ihx] = env[kres2D(iz,ihx,ihz)];
                        }
                    }
                    this->find_max(hwrk, &hmax[iz], &imax[iz], nhx);
                }

                // Threshold 
                for (iz=0; iz<nz; iz++){
                    hsort[iz] = hmax[iz];
                }
                std::sort(hsort, hsort+nz); 
                pos = (int) (THRES*nz/100);
                pclip = hsort[pos];
                for (iz=0; iz<nz; iz++){
                    if(hmax[iz] < pclip){
                        hmax[iz] = 0;
                        imax[iz] = (nhx-1)/2;
                    }
                }

                // Calculate second derivative
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                    }
                }

                // Calculate double Hilbert transform
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();
                hilb_cip2->hilbertx(imag);
                imag = hilb_cip2->getDf();

                //Build residual
                for (iz=0; iz<nz; iz++){
                    den = 0.0;
                    num = 0.0;
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        den = (env[kres2D(iz,imax[iz]+2,0)] - 2.0*env[kres2D(iz,imax[iz],0)] + env[kres2D(iz,imax[iz]-2,0)])/4.0;
                    }
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        num = -2.0*(cip[kres2D(iz,imax[iz]+1,0)] - cip[kres2D(iz,imax[iz]-1,0)])/2.0;
                        num += 2.0*(imag[kres2D(iz,imax[iz]+1,0)] - imag[kres2D(iz,imax[iz]-1,0)])/2.0;
                    }

                    // Zero out residual data not in hmax
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = 0.0;
                        }
                    }

                    if(den != 0.0){
                        imagedata[ki2D(ix,iz,imax[iz],0)] = (f3/f2)*(imax[iz]-((nhx-1)/2))*num*den/(SQ(den) + (rms*1.0e-6));
                    }
                    // SI residual
                    imagedata[ki2D(ix,iz,(nhx-1)/2,(nhz-1)/2)] += (f1*f3/(f2*f2))*cip[kres2D(iz,(nhx-1)/2,(nhz-1)/2)];

                    // DS resudual
                    for (ihx=0; ihx<nhx; ihx++){
                        hx= -(nhx-1)/2 + ihx;
                        G1 = (hx*hx);
                        for (ihz=0; ihz<nhz; ihz++){
                            hz= -(nhz-1)/2 + ihz;
                            G2 = G1 + (hz*hz);
                            imagedata[ki2D(ix,iz,ihx,ihz)] -= (f1/f2)*G2*cip[kres2D(iz,ihx,ihz)];

                        }
                    }
                }
            }
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
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
    free(cip);
    free(env);
    free(hmax);
    free(hsort);
    free(hwrk);
    free(imax);

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
void KdmvaAcoustic2D<T>::readMisfit(double *f)
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
void KdmvaAcoustic2D<T>::readGrad(double *g)
{
    int i;
    int N;
    float *g_in;
    T *gvp;
    std::string vpgradfile;
    if(Modelmutefile.empty()){
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
            rs_error("KdmvaAcoustic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void KdmvaAcoustic2D<T>::combineGradients()
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
void KdmvaAcoustic2D<T>::applyMute()
{
    if(!Modelmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> model;
        model = std::make_shared<rockseis::ModelAcoustic2D<T>>(VPGRADCOMBFILE, VPGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute (new rockseis::ModelAcoustic2D<T>(Modelmutefile, Modelmutefile, 1 ,0));

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
void KdmvaAcoustic2D<T>::computeRegularisation(double *x)
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
            rs_error("KdmvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
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

// =============== 2D ELASTIC KDMVA CLASS =============== //
// 
template<typename T>
KdmvaElastic2D<T>::KdmvaElastic2D() {
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
KdmvaElastic2D<T>::KdmvaElastic2D(MPImodeling *mpi): Kdmva<T>(mpi) {
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
KdmvaElastic2D<T>::~KdmvaElastic2D() {
    //Do nothing
}

template<typename T>
void KdmvaElastic2D<T>::setWavemode(int mode) {
    if(mode != 0 && mode != 1) rs_error("KdmvaElastic2D<T>::setWavemod: Invalid wavemode.");
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
void KdmvaElastic2D<T>::runPPgrad() {
    MPImodeling *mpi = this->getMpi();
    // Calculate traveltimes
    // Migrate data
    // Compute Adjoint
}

template<typename T>
void KdmvaElastic2D<T>::runPSgrad() {
    MPImodeling *mpi = this->getMpi();
    // Calculate traveltimes
    // Migrate data
    // Compute Adjoint
}


template<typename T>
void KdmvaElastic2D<T>::runBsproj() {
    MPImodeling *mpi = this->getMpi();
	T vpsum = 0.0; // Sum over splines
	T vssum = 0.0; // Sum over splines
	float *global_stack;
    T *c;
    T *wrk;

    std::string vpgradfile;
    std::string vsgradfile;

    if(Modelmutefile.empty()){
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
        rs_error("KdmvaElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
    }
    vsproj= (float *) calloc(nc, sizeof(float));
    if(vsproj==NULL){
        rs_error("KdmvaElastic2D<T>::runBsproj(): Not enough memory to allocate projection array (vsproj)");
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
			rs_error("KdmvaElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
            rs_error("KdmvaElastic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
int KdmvaElastic2D<T>::setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile)
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
            rs_error("KdmvaElastic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return Npar;
}

template<typename T>
void KdmvaElastic2D<T>::saveResults(int iter)
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
void KdmvaElastic2D<T>::saveLinesearch(double *x)
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
    if(!Modelmutefile.empty()){
        mute = std::make_shared <rockseis::ModelElastic2D<T>>(Modelmutefile, Modelmutefile, Modelmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("KdmvaElastic2D<T>::saveLinesearch(): Geometry in Modelmutefile does not match geometry in the model.");
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
            rs_error("KdmvaAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
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
    if(Modelmutefile.empty()){
        free(vpmutedata);
        free(vsmutedata);
    }
}

template<typename T>
void KdmvaElastic2D<T>::readMisfit(double *f)
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
void KdmvaElastic2D<T>::readGrad(double *g)
{
    int i;
    int N,Npar=0;
    float *g_in;
    T *gvp, *gvs;
    std::string vpgradfile;
    std::string vsgradfile;
    if(Modelmutefile.empty()){
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
            rs_error("KdmvaElastic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void KdmvaElastic2D<T>::combineGradients()
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
void KdmvaElastic2D<T>::applyMute()
{
    if(!Modelmutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelElastic2D<T>> model;
        model = std::make_shared<rockseis::ModelElastic2D<T>>(VPGRADCOMBFILE, VSGRADCOMBFILE, VPGRADCOMBFILE, 1 ,0);
        // Mute
        std::shared_ptr<rockseis::ModelElastic2D<T>> mute (new rockseis::ModelElastic2D<T>(Modelmutefile, Modelmutefile, Modelmutefile, 1 ,0));

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
void KdmvaElastic2D<T>::computeRegularisation(double *x)
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
            rs_error("KdmvaElastic2D<T>::computeRegularization(): Unknown parameterisation."); 
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
void KdmvaElastic2D<T>::computeMisfit(std::shared_ptr<rockseis::Image2D<T>> image, std::string imageresfile)
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
    T f3 = 0.;

    std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute;
    int N;
    T *mutedata;
    // If mute
    if(!Residualmutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Residualmutefile, Residualmutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (nx*nz);
        if(N != Nmute) rs_error("KdmvaAcoustic2D<T>::computeMisfit(): Geometry in Residualmutefile does not match geometry in the image.");
        mute->readModel();
        mutedata = mute->getVp();
    }else{
        N = (nx*nz);
        mutedata = (T *) calloc(N, sizeof(T)); 
        for(ix=0; ix < N; ix++){
            mutedata[ix] = 1.0;
        }
    }

    /* Variables for DS_hmax misfit and residual computation*/
    T *cip = (T *) calloc(nz*nhx*nhz, sizeof(T));
    T *env = (T *) calloc(nz*nhx*nhz, sizeof(T));
    T *hmax = (T *) calloc(nz, sizeof(T));
    T *hsort = (T *) calloc(nz, sizeof(T));
    T *hwrk = (T *) calloc(nhx, sizeof(T));
    int *imax = (int *) calloc(nz, sizeof(int));
    std::shared_ptr<Hilbert<T>> hilb_cip1 (new Hilbert<T>(nz, nhx, nhz, nz));
    std::shared_ptr<Hilbert<T>> hilb_cip2 (new Hilbert<T>(nz, nhx, nhz, nz));
    T pclip;
    T *imag;
    int pos;
    T num, den;
    T rms=0;
    size_t fold=0;


    switch(this->getMisfit_type()){
        case SI:
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
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
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
                    }
                }
            }
            break;
        case DS_PLUS_SI:
            // Misfit
            for (ihx=0; ihx<nhx; ihx++){
                hx= -(nhx-1)/2 + ihx;
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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
                G1 = WEIGHT((T) hx, (T) ((nhx-1)/2));
                for (ihz=0; ihz<nhz; ihz++){
                    hz= -(nhz-1)/2 + ihz;
                    G2 = G1*WEIGHT((T) hz, (T) ((nhz-1)/2));
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
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
                    }
                }
            }
            break;
        case DS_HMAX:
            // Misfit
            // Get a CIP gather
            for (ix=0; ix<nx; ix++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        // Derivative
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                    }
                }
                // Hilbert transform over first non-singleton axis
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();

                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=0; iz<nz; iz++){
                           env[kres2D(iz,ihx,ihz)] = SQ(cip[kres2D(iz,ihx,ihz)]) + SQ(imag[kres2D(iz,ihx,ihz)]);
                        }
                    }
                }
                // Find maximum and index of maximum along hx axis and store in an array
                for (iz=0; iz<nz; iz++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            hwrk[ihx] = env[kres2D(iz,ihx,ihz)];
                        }
                    }
                    this->find_max(hwrk, &hmax[iz], &imax[iz], nhx);
                }

                // Threshold 
                for (iz=0; iz<nz; iz++){
                    hsort[iz] = hmax[iz];
                }
                std::sort(hsort, hsort+nz); 
                pos = (int) (THRES*nz/100);
                pclip = hsort[pos];
                for (iz=0; iz<nz; iz++){
                    if(hmax[iz] < pclip){
                        hmax[iz] = 0;
                        imax[iz] = (nhx-1)/2;
                    }
                }
                // Calculate delta h, and misfit. 
                for (iz=0; iz<nz; iz++){
                    f1 += 0.5*SQ(imax[iz]-((nhx-1)/2));
                    f2 += 0.5*SQ(cip[kres2D(iz,(nhx-1)/2,(nhz-1)/2)]);
                }

                for (iz=0; iz<nz; iz++){
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        rms += SQ((env[kres2D(iz,imax[iz]+2,0)] - 2.0*env[kres2D(iz,imax[iz],0)] + env[kres2D(iz,imax[iz]-2,0)])/4.0);
                        fold++;
                    }
                }
            }
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
                            f3 += 0.5*G2*wrk[iz]*wrk[iz];
                        }
                                            }
                }
            }

            f = f1*f3/f2;

            if(fold >0){
                rms = (T) std::sqrt(rms)/fold;
            }

            // Residual
            for (ix=0; ix<nx; ix++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        // Derivative
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - imagedata[ki2D(ix,iz,ihx,ihz)];
                        }
                    }
                }
                // Hilbert transform over first non-singleton axis
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();

                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=0; iz<nz; iz++){
                           env[kres2D(iz,ihx,ihz)] = SQ(cip[kres2D(iz,ihx,ihz)]) + SQ(imag[kres2D(iz,ihx,ihz)]);
                        }
                    }
                }
                // Find maximum and index of maximum along hx axis and store in an array
                for (iz=0; iz<nz; iz++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            hwrk[ihx] = env[kres2D(iz,ihx,ihz)];
                        }
                    }
                    this->find_max(hwrk, &hmax[iz], &imax[iz], nhx);
                }

                // Threshold 
                for (iz=0; iz<nz; iz++){
                    hsort[iz] = hmax[iz];
                }
                std::sort(hsort, hsort+nz); 
                pos = (int) (THRES*nz/100);
                pclip = hsort[pos];
                for (iz=0; iz<nz; iz++){
                    if(hmax[iz] < pclip){
                        hmax[iz] = 0;
                        imax[iz] = (nhx-1)/2;
                    }
                }

                // Calculate second derivative
                for (ihx=0; ihx<nhx; ihx++){
                    for (ihz=0; ihz<nhz; ihz++){
                        for (iz=1; iz<nz-1; iz++){
                            cip[kres2D(iz,ihx,ihz)] = imagedata[ki2D(ix,iz+1,ihx,ihz)] - 2.0*imagedata[ki2D(ix,iz,ihx,ihz)] + imagedata[ki2D(ix,iz-1,ihx,ihz)];
                        }
                    }
                }

                // Calculate double Hilbert transform
                hilb_cip1->hilbertx(cip);
                imag = hilb_cip1->getDf();
                hilb_cip2->hilbertx(imag);
                imag = hilb_cip2->getDf();

                //Build residual
                for (iz=0; iz<nz; iz++){
                    den = 0.0;
                    num = 0.0;
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        den = (env[kres2D(iz,imax[iz]+2,0)] - 2.0*env[kres2D(iz,imax[iz],0)] + env[kres2D(iz,imax[iz]-2,0)])/4.0;
                    }
                    if(imax[iz] >1 && imax[iz] < nhx-1){
                        num = -2.0*(cip[kres2D(iz,imax[iz]+1,0)] - cip[kres2D(iz,imax[iz]-1,0)])/2.0;
                        num += 2.0*(imag[kres2D(iz,imax[iz]+1,0)] - imag[kres2D(iz,imax[iz]-1,0)])/2.0;
                    }

                    // Zero out residual data not in hmax
                    for (ihz=0; ihz<nhz; ihz++){
                        for (ihx=0; ihx<nhx; ihx++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] = 0.0;
                        }
                    }

                    if(den != 0.0){
                        imagedata[ki2D(ix,iz,imax[iz],0)] = (f3/f2)*(imax[iz]-((nhx-1)/2))*num*den/(SQ(den) + (rms*1.0e-6));
                    }
                    // SI residual
                    imagedata[ki2D(ix,iz,(nhx-1)/2,(nhz-1)/2)] += (f1*f3/(f2*f2))*cip[kres2D(iz,(nhx-1)/2,(nhz-1)/2)];

                    // DS resudual
                    for (ihx=0; ihx<nhx; ihx++){
                        hx= -(nhx-1)/2 + ihx;
                        G1 = (hx*hx);
                        for (ihz=0; ihz<nhz; ihz++){
                            hz= -(nhz-1)/2 + ihz;
                            G2 = G1 + (hz*hz);
                            imagedata[ki2D(ix,iz,ihx,ihz)] -= (f1/f2)*G2*cip[kres2D(iz,ihx,ihz)];

                        }
                    }
                }
            }
            // Apply mute
            for (ihz=0; ihz<nhz; ihz++){
                for (ihx=0; ihx<nhx; ihx++){
                    for (iz=0; iz<nz; iz++){
                        for (ix=0; ix<nx; ix++){
                            imagedata[ki2D(ix,iz,ihx,ihz)] *= mutedata[km2D(ix,iz)];
                        }
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
template class Kdmva<float>;
template class Kdmva<double>;

template class KdmvaAcoustic2D<float>;
template class KdmvaAcoustic2D<double>;

template class KdmvaElastic2D<float>;
template class KdmvaElastic2D<double>;

}
