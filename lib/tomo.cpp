#include "tomo.h"


namespace rockseis {

// =============== ABSTRACT TOMO CLASS =============== //
template<typename T>
Tomo<T>::Tomo() {
    //Set default parameters
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    paramtype = PAR_GRID;
    dtx = -1;
    dty = -1;
    dtz = -1;
    fnorm = 0.0;
    if(createLog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating logfile for writting.");
    }
    if(createProglog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
Tomo<T>::Tomo(MPImodeling *_mpi) {
    mpi = _mpi;

    //Set default parameters
    lpml = 10;
    incore = true;
    order = 4;
    snapinc=4;
    nsnaps=0;
    paramtype = PAR_GRID;
    dtx = -1;
    dty = -1;
    dtz = -1;
    fnorm = 0.0;
    if(createLog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating logfile for writting.");
    }
    if(createProglog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating progress logfile for writting.");
    }
    noreverse = false;
}

template<typename T>
void Tomo<T>::normalize(double *v, double *f, int n){
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
double Tomo<T>::vector_norm(double *v, const int type, const int n){
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
bool Tomo<T>::createLog(){
	logfile = LOGFILE;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return TOM_ERR;
	}else{
		Flog.close();
		return TOM_OK;
	}
}

template<typename T>
bool Tomo<T>::createProglog(){
	progresslogfile = PROGLOGFILE;
	Flog.open(progresslogfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return TOM_ERR;
	}else{
		Flog.close();
		return TOM_OK;
	}
}

template<typename T>
void Tomo<T>::writeLog(std::string text){
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
void Tomo<T>::writeProgress(std::string text){
    if(!progresslogfile.empty()){
        Flog.open(progresslogfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Tomo<T>::createResult(){
    struct stat s;
    // Checking if result folder is present, and creates it if not
    if(stat(RESULTDIR,&s) != 0) {
        int mkdir_return = mkdir(RESULTDIR,0777);
        if(mkdir_return != 0) rs_error("Not able to create result directory: ", RESULTDIR);
    }
}

template<typename T>
Tomo<T>::~Tomo() {
    //Do nothing
}

// =============== 2D ACOUSTIC TOMO CLASS =============== //
//
template<typename T>
TomoAcoustic2D<T>::TomoAcoustic2D() {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;

    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
TomoAcoustic2D<T>::TomoAcoustic2D(MPImodeling *mpi): Tomo<T>(mpi) {
    // Set default parameters
    apertx = -1;

    kvp = 1.0;
    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
TomoAcoustic2D<T>::~TomoAcoustic2D() {
    //Do nothing
}


template<typename T>
void TomoAcoustic2D<T>::clipGrad(std::shared_ptr<rockseis::Image2D<T>> grad)
{
    int nx, nz;
    nx = grad->getNx();
    nz = grad->getNz();
	if(!grad->getAllocated()) grad->allocateImage();
    grad->read();
    T * data = grad->getImagedata();
    T *thres = (T *) calloc(nx*nz, sizeof(T));
    for (int i=0; i<nx*nz; i++){
        thres[i] = ABS(data[i]);
    }
    std::sort(thres, thres+(nx*nz)); 
    int pos = (int) (97.0*(nx*nz)/100.0);
    T pclip = thres[pos];

    for (int i=0; i<nx*nz; i++){
        if(ABS(data[i]) >= pclip){
            data[i] = 0.0;
        }
    }
    grad->write();
    free(thres);
}


template<typename T>
void TomoAcoustic2D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data2D<T>> Tobs2D;
    std::shared_ptr<rockseis::Data2D<T>> Tmod2D;
    std::shared_ptr<rockseis::Data2D<T>> Tres2D;
    std::shared_ptr<rockseis::Data2D<T>> weight2D;
    std::shared_ptr<rockseis::Image2D<T>> vpgrad;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Trecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> gmodel (new rockseis::ModelAcoustic2D<T>(Vpfile, Vpfile, 0 ,false));

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());


	if(mpi->getRank() == 0) {
		// Master

        // Get shot map
        Sort->readKeymap();
        Sort->readSortmap();
        size_t ngathers =  Sort->getNensemb();

        // Create a data class for the recorded data
        std::shared_ptr<rockseis::Data2D<T>> Tobs2D (new rockseis::Data2D<T>(Trecordfile));
        // Create modelling and residual data files
        Tmod2D = std::make_shared<rockseis::Data2D<T>>(1, Tobs2D->getNt(), Tobs2D->getDt(), Tobs2D->getOt());
        Tmod2D->setFile(Tmodelledfile);
        Tmod2D->createEmpty(Tobs2D->getNtrace());

        Tres2D = std::make_shared<rockseis::Data2D<T>>(1, Tobs2D->getNt(), Tobs2D->getDt(), Tobs2D->getOt());
        Tres2D->setFile(Tresidualfile);
        Tres2D->createEmpty(Tobs2D->getNtrace());

        // Misfit file creation
        Fmisfit->output(Misfitfile);
        Fmisfit->setN(1,ngathers);
        Fmisfit->setD(1,1.0);
        Fmisfit->setData_format(sizeof(T));
        Fmisfit->createEmpty();
        Fmisfit->close();

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

        // Clip extreme values from gradient
        //this->clipGrad(vpgrad);

		//Clear work vector 
		mpi->clearWork();

    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::FatAcoustic2D<T>> fat;
        while(1) {
            workModeling_t work = mpi->receiveWork();

            if(work.MPItag == MPI_TAG_DIE) {
                break;
            }

            if(work.MPItag == MPI_TAG_NO_WORK) {
                mpi->sendNoWork(0);
            }
            else {
                // Do gradient computation
                Sort->readKeymap();
                Sort->readSortmap();

                // Get the shot
                Sort->setDatafile(Trecordfile);
                Tobs2D = Sort->get2DGather(work.id);
                size_t ntr = Tobs2D->getNtrace();

                lmodel = gmodel->getLocal(Tobs2D, -3.0*gmodel->getDx(), SMAP);

                //Create source 
                std::shared_ptr<rockseis::Data2D<T>> source (new rockseis::Data2D<T>(ntr, 1, 1.0, 0.0));
                source->copyCoords(Tobs2D);
                source->makeMap(lmodel->getGeom(), SMAP);


                // Create first arrival tomography object
                 fat = std::make_shared<rockseis::FatAcoustic2D<T>>(lmodel, source, Tobs2D);
                // Create modelled and residual data objects 
                Tmod2D = std::make_shared<rockseis::Data2D<T>>(ntr, Tobs2D->getNt(), Tobs2D->getDt(), 0.0);
                Tmod2D->copyCoords(Tobs2D);
                Tmod2D->makeMap(lmodel->getGeom());
                fat->setTmod(Tmod2D);
                Tres2D = std::make_shared<rockseis::Data2D<T>>(ntr, Tobs2D->getNt(), Tobs2D->getDt(), 0.0);
                Tres2D->copyCoords(Tobs2D);
                Tres2D->makeMap(lmodel->getGeom());
                fat->setTres(Tres2D);

                if(dataweight){
                    Sort->setDatafile(Dataweightfile);
                    weight2D = Sort->get2DGather(work.id);
                    fat->setTweight(weight2D);
                }

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image2D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1);
                fat->setVpgrad(vpgrad);

                // Set logfile
                fat->setLogfile("log.txt-" + std::to_string(work.id));

                /////// Run eikonal solver
                fat->run();

                // Output gradient
                vpgrad->write();

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fat->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                Tmod2D->setFile(Tmodelledfile);
                Sort->put2DGather(Tmod2D, work.id);
                Tres2D->setFile(Tresidualfile);
                Sort->put2DGather(Tres2D, work.id);

                
                // Reset all classes
                Tobs2D.reset();
                Tmod2D.reset();
                Tres2D.reset();
                if(dataweight){
                    weight2D.reset();
                }
                source.reset();
                lmodel.reset();
                vpgrad.reset();
                fat.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void TomoAcoustic2D<T>::runBsproj() {
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
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad (new rockseis::ModelAcoustic2D<T>(vpgradfile, vpgradfile, this->getLpml() ,false));

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
		rs_error("TomoAcoustic2D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
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
			rs_error("TomoAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
			rs_error("TomoAcoustic2D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
        MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

	   }
    // Free variables
    free(global_stack);
    free(vpproj);
}

template<typename T>
int TomoAcoustic2D<T>::setInitial(double *x, std::string vpfile)
{
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model_in (new rockseis::ModelAcoustic2D<T>(vpfile, vpfile, 1 ,0));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    model_in->readModel();  
    // Write initial model files
    model_in->setVpfile(VP0FILE);
    model_in->writeVp();
    // Write linesearch model files
    model_in->setVpfile(VPLSFILE);
    model_in->writeVp();
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
            rs_error("TomoAcoustic2D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return N;
}

template<typename T>
void TomoAcoustic2D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, VPLSFILE, 1 ,0));
    lsmodel->readModel();
    name = dir + "/" + VP_UP + "-" + std::to_string(iter);
    lsmodel->setVpfile(name);
    lsmodel->writeVp();
}

template<typename T>
void TomoAcoustic2D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model0 (new rockseis::ModelAcoustic2D<T>(VP0FILE, VP0FILE, 1 ,0));
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> lsmodel (new rockseis::ModelAcoustic2D<T>(VPLSFILE, VPLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> mute;

    // Write linesearch model
    model0->readModel();
    lsmodel->readModel();
    T *vp0, *vpls;
    T *c, *mod;
    T *vpmutedata;
    vp0 = model0->getVp(); 
    vpls = lsmodel->getVp(); 
    int i;
    int N, Nmod;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelAcoustic2D<T>>(Mutefile, Mutefile, 1 ,0);
        int Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("TomoAcoustic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
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
            }
            lsmodel->writeVp();
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
            lsmodel->writeVp();
            break;
        default:
            rs_error("TomoAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Mutefile.empty()){
        free(vpmutedata);
    }
}

template<typename T>
void TomoAcoustic2D<T>::readMisfit(double *f)
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
void TomoAcoustic2D<T>::readGrad(double *g)
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
            rs_error("TomoAcoustic2D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void TomoAcoustic2D<T>::combineGradients()
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
void TomoAcoustic2D<T>::applyMute()
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
void TomoAcoustic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> model (new rockseis::ModelAcoustic2D<T>(VPLSFILE, VPLSFILE, 1 ,0));
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
            rs_error("TomoAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
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





// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Tomo<float>;
template class Tomo<double>;

template class TomoAcoustic2D<float>;
template class TomoAcoustic2D<double>;

}