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
    constrain = false;
    if(createLog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating logfile for writting.");
    }
    if(createProglog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating progress logfile for writting.");
    }
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
    constrain = false;
    if(createLog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating logfile for writting.");
    }
    if(createProglog() == TOM_ERR)
    {
        rs_error("Tomo<T>::Tomo(): Error creating progress logfile for writting.");
    }
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
void Tomo<T>::un_normalize(double *v, double f, int n){
	int i;
	for(i=0; i<n; i++) {
		v[i] *= fabs(f);
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
    std::shared_ptr<rockseis::ModelEikonal2D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Trecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelEikonal2D<T>> gmodel (new rockseis::ModelEikonal2D<T>(Vpfile, this->getLpml()));

     // Test for problematic model sampling
    if(gmodel->getDx() != gmodel->getDz()){
        rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dz).");
    }

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

        // Stack gradient
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

                lmodel = gmodel->getLocal(Tobs2D, -1.0*gmodel->getDx(), SMAP);
                lmodel->Expand();

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

                // Stack gradient
                vpgrad->write();
                
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

                // Send result back
                work.status = WORK_FINISHED;
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
   std::shared_ptr<rockseis::ModelEikonal2D<T>> grad (new rockseis::ModelEikonal2D<T>(vpgradfile, this->getLpml()));

   // Read model
   grad->readVelocity();

   T *vpgrad;
   vpgrad = grad->getVelocity();

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
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> bound;
    std::shared_ptr<rockseis::ModelEikonal2D<T>> model_in (new rockseis::ModelEikonal2D<T>(vpfile, 1));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;

    // Read initial model
    model_in->readVelocity();  

    // Check bounds 
    if(this->getConstrain()){
       bound = std::make_shared <rockseis::ModelAcoustic2D<T>>(Lboundfile, Uboundfile, 1 ,0);
       T *lbounddata;
       T *ubounddata;
       T *vp0;
       bound->readModel();
       long Nbound = (bound->getGeom())->getNtot();
       long N = (long) (model_in->getGeom())->getNtot();
       if(N != Nbound) rs_error("TomoAcoustic2D<T>::setInitial(): Geometry in Boundary files does not match geometry in the model.");
       lbounddata = bound->getVp();
       ubounddata = bound->getR();
       vp0 = model_in->getVelocity();
       for (long i=0; i< N; i++){
          if(lbounddata[i] == ubounddata[i]) rs_error("TomoAcoustic2D<T>::setInitial():The lower bound cannot be equal to the upper bound, use a mute function insted.");
          if(vp0[i] <= lbounddata[i]) vp0[i] = lbounddata[i] + 1e-2;
          if(vp0[i] >= ubounddata[i]) vp0[i] = ubounddata[i] - 1e-2;
       }
    }

    // Write initial model files
    model_in->setVelocityfile(VP0FILE);
    model_in->writeVelocity();
    // Write linesearch model files
    model_in->setVelocityfile(VPLSFILE);
    model_in->writeVelocity();
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
    std::shared_ptr<rockseis::ModelEikonal2D<T>> lsmodel (new rockseis::ModelEikonal2D<T>(VPLSFILE, 1));
    lsmodel->readVelocity();
    name = dir + "/" + VP_UP + "-" + std::to_string(iter);
    lsmodel->setVelocityfile(name);
    lsmodel->writeVelocity();
}

template<typename T>
void TomoAcoustic2D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelEikonal2D<T>> model0 (new rockseis::ModelEikonal2D<T>(VP0FILE, 1));
    std::shared_ptr<rockseis::ModelEikonal2D<T>> lsmodel (new rockseis::ModelEikonal2D<T>(VPLSFILE, 1));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelEikonal2D<T>> mute;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> bound;

    // Write linesearch model
    model0->readVelocity();
    lsmodel->readVelocity();
    T *vp0, *vpls;
    T *c, *mod;
    T *vpmutedata;
    T *lbounddata;
    T *ubounddata;
    T *x0;
    double log_in, log_out, exp_in, exp_out;

    vp0 = model0->getVelocity(); 
    vpls = lsmodel->getVelocity(); 
    int i;
    int N, Nmod;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelEikonal2D<T>>(Mutefile, 1);
        long Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("TomoAcoustic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readVelocity();
        vpmutedata = mute->getVelocity();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
        }
    }
    if(this->getConstrain()){
        bound = std::make_shared <rockseis::ModelAcoustic2D<T>>(Lboundfile, Uboundfile, 1 ,0);
        bound->readModel();
        long Nbound = (bound->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nbound) rs_error("TomoAcoustic2D<T>::saveLinesearch(): Geometry in Boundary files does not match geometry in the model.");
        lbounddata = bound->getVp();
        ubounddata = bound->getR();
        x0 = (T *) calloc(N, sizeof(T));
    }
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
               if(!this->getConstrain()){
                  vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
               }else{
                  log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
                  log_out = log(log_in);
                  x0[i] =((T) log_out);
                  exp_in = (double) (-(x[i]*vpmutedata[i]*kvp + x0[i]));
                  exp_out = exp(exp_in);
                  vpls[i] = lbounddata[i] + (ubounddata[i]-lbounddata[i])*(1.0/(1.0 + (T) exp_out));
               }
            }
            lsmodel->writeVelocity();
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
               if(!this->getConstrain()){
                  vpls[i] = vp0[i] + mod[i]*vpmutedata[i]*kvp;
               }else{
                  log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
                  log_out = log(log_in);
                  x0[i] =((T) log_out);
                  exp_in = (double) (-(mod[i]*vpmutedata[i]*kvp + x0[i]));
                  exp_out = exp(exp_in);
                  vpls[i] = lbounddata[i] + (ubounddata[i]-lbounddata[i])*(1.0/(1.0 + (T) exp_out));
               }
            }
            lsmodel->writeVelocity();
            break;
        default:
            rs_error("TomoAcoustic2D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Mutefile.empty()){
       free(vpmutedata);
    }
    if(this->getConstrain()){
       free(x0);
    }
}

template<typename T>
void TomoAcoustic2D<T>::saveHessian(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelEikonal2D<T>> lsmodel (new rockseis::ModelEikonal2D<T>(VPLSFILE, 1));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::ModelEikonal2D<T>> mute;
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> bound;

    // Write linesearch model
    lsmodel->readVelocity();
    T *c, *mod;
    T *vpls;
    T *vpmutedata;
    T *lbounddata;
    T *ubounddata;
    vpls = lsmodel->getVelocity(); 
    lsmodel->setVelocityfile(VPHESSFILE);
    int i;
    int N, Nmod;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelEikonal2D<T>>(Mutefile, 1);
        long Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("TomoAcoustic2D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readVelocity();
        vpmutedata = mute->getVelocity();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
        }
    }
    if(this->getConstrain()){
        bound = std::make_shared <rockseis::ModelAcoustic2D<T>>(Lboundfile, Uboundfile, 1, 0);
        bound->readModel();
        long Nbound = (bound->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nbound) rs_error("TomoAcoustic2D<T>::saveLinesearch(): Geometry in Boundary files does not match geometry in the model.");
        lbounddata = bound->getVp();
        ubounddata = bound->getR();
    }
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
                vpls[i] = x[i]*vpmutedata[i]*kvp*kvp;
                if(this->getConstrain()){
                    if(vpls[i] < lbounddata[i]) vpls[i] = lbounddata[i];
                    if(vpls[i] > ubounddata[i]) vpls[i] = ubounddata[i];
                }
            }
            lsmodel->writeVelocity();
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
                if(this->getConstrain()){
                    if(vpls[i] < lbounddata[i]) vpls[i] = lbounddata[i];
                    if(vpls[i] > ubounddata[i]) vpls[i] = ubounddata[i];
                }
            }
            lsmodel->writeVelocity();
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

    std::shared_ptr<rockseis::ModelEikonal2D<T>> modelgrad (new rockseis::ModelEikonal2D<T>(vpgradfile, 1));
    std::shared_ptr<rockseis::Bspl2D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readVelocity();
            N = (modelgrad->getGeom())->getNtot();
            gvp = modelgrad->getVelocity(); 
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
   std::shared_ptr<rockseis::ModelEikonal2D<T>> grad;
   std::shared_ptr<rockseis::ModelEikonal2D<T>> reggrad;
   grad = std::make_shared<rockseis::ModelEikonal2D<T>>(VPGRADFILE, 1);
   reggrad = std::make_shared<rockseis::ModelEikonal2D<T>>(VPREGGRADFILE, 1);

   // Read gradients
   grad->readVelocity();
   reggrad->readVelocity();
   T *vp, *regvp;
   vp = grad->getVelocity(); 
   regvp = reggrad->getVelocity(); 
   int i;
   int N;
   N = (grad->getGeom())->getNtot();
   
   // Compute 
   for(i=0; i< N; i++)
   {
      vp[i] = vp[i] + reg_alpha[0]*regvp[i];
   }
   grad->setVelocityfile(VPGRADCOMBFILE);
   grad->writeVelocity();
}


template<typename T>
void TomoAcoustic2D<T>::applyChainrule(double *x)
{
   // Models
   std::shared_ptr<rockseis::ModelEikonal2D<T>> model0 (new rockseis::ModelEikonal2D<T>(VP0FILE, 1));
   std::shared_ptr<rockseis::ModelEikonal2D<T>> grad (new rockseis::ModelEikonal2D<T>(VPGRADCOMBFILE, 1));
   std::shared_ptr<rockseis::Bspl2D<T>> spline;
   std::shared_ptr<rockseis::ModelEikonal2D<T>> mute;
   std::shared_ptr<rockseis::ModelAcoustic2D<T>> bound;

   // Write linesearch model
   model0->readVelocity();
   grad->readVelocity();
   T *vp0, *vpgrad;
   T *c, *mod;
   T *vpmutedata;
   T *lbounddata;
   T *ubounddata;
   T *x0;
   double log_in, log_out, exp_in, exp_out;

   vp0 = model0->getVelocity(); 
   vpgrad = grad->getVelocity(); 
   int i;
   int N, Nmod;

   // If mute
   if(!Mutefile.empty()){
      mute = std::make_shared <rockseis::ModelEikonal2D<T>>(Mutefile, 1);
      long Nmute = (mute->getGeom())->getNtot();
      N = (grad->getGeom())->getNtot();
      if(N != Nmute) rs_error("TomoAcoustic2D<T>::applyChainrule(): Geometry in Mutefile does not match geometry in the model.");
      mute->readVelocity();
      vpmutedata = mute->getVelocity();
   }else{
      N = (grad->getGeom())->getNtot();
      vpmutedata = (T *) calloc(N, sizeof(T)); 
      for(i=0; i < N; i++){
         vpmutedata[i] = 1.0;
      }
   }
   bound = std::make_shared <rockseis::ModelAcoustic2D<T>>(Lboundfile, Uboundfile, 1 ,0);
   bound->readModel();
   long Nbound = (bound->getGeom())->getNtot();
   N = (grad->getGeom())->getNtot();
   if(N != Nbound) rs_error("TomoAcoustic2D<T>::applyChainrule(): Geometry in Boundary files does not match geometry in the model.");
   lbounddata = bound->getVp();
   ubounddata = bound->getR();
   x0 = (T *) calloc(N, sizeof(T));
   switch (this->getParamtype()){
      case PAR_GRID:
         N = (grad->getGeom())->getNtot();
         for(i=0; i< N; i++)
         {
            log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
            log_out = log(log_in);
            x0[i] =((T) log_out);
            exp_in = (double) (-(x[i]*vpmutedata[i]*kvp + x0[i]));
            exp_out = exp(exp_in);
            vpgrad[i] *= kvp*(ubounddata[i]-lbounddata[i])*((T) exp_out/((1.0 + (T) exp_out)*(1.0 + (T) exp_out)));
         }
         grad->writeVelocity();
         break;
      case PAR_BSPLINE:
         Nmod = (grad->getGeom())->getNtot();
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
            log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
            log_out = log(log_in);
            x0[i] =((T) log_out);
            exp_in = (double) (-(mod[i]*vpmutedata[i]*kvp + x0[i]));
            exp_out = exp(exp_in);
            vpgrad[i] *= kvp*(ubounddata[i]-lbounddata[i])*((T) exp_out/((1.0 + (T) exp_out)*(1.0 + (T) exp_out)));
         }
         grad->writeVelocity();
         break;
      default:
         rs_error("TomoAcoustic2D<T>::applyChainrule(): Unknown parameterisation."); 
         break;
   }
   if(Mutefile.empty()){
      free(vpmutedata);
   }
   free(x0);
}

template<typename T>
void TomoAcoustic2D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelEikonal2D<T>> model;
        model = std::make_shared<rockseis::ModelEikonal2D<T>>(VPGRADCOMBFILE, 1);
        // Mute
        std::shared_ptr<rockseis::ModelEikonal2D<T>> mute (new rockseis::ModelEikonal2D<T>(Mutefile, 1));

        // Mute model and write
        model->readVelocity();
        mute->readVelocity();
        T *vp, *vpmute;
        vp = model->getVelocity(); 
        vpmute = mute->getVelocity(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
        }
        model->setVelocityfile(VPGRADMUTEFILE);
        model->writeVelocity();
    }
}

template<typename T>
void TomoAcoustic2D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelEikonal2D<T>> model (new rockseis::ModelEikonal2D<T>(VPLSFILE, 1));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), 1, model->getNz(), model->getDx(), 1.0, model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl2D<double>> spline;

    // Write linesearch model
    model->readVelocity();
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
    model->readVelocity();
    model->setVelocityfile(VPREGGRADFILE);
    vpgrad = model->getVelocity();
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
    model->writeVelocity();

    // Free variables
    free(dvpdx);
    free(dvpdz);
    free(gwrk);
}

// =============== 3D ACOUSTIC TOMO CLASS =============== //
//
template<typename T>
TomoAcoustic3D<T>::TomoAcoustic3D() {
    // Set default parameters
    apertx = -1;
    aperty = -1;

    kvp = 1.0;

    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
TomoAcoustic3D<T>::TomoAcoustic3D(MPImodeling *mpi): Tomo<T>(mpi) {
    // Set default parameters
    apertx = -1;
    aperty = -1;

    kvp = 1.0;
    reg_alpha[0]=0.0;
    reg_eps[0]=1e-3;
}

template<typename T>
TomoAcoustic3D<T>::~TomoAcoustic3D() {
    //Do nothing
}


template<typename T>
void TomoAcoustic3D<T>::clipGrad(std::shared_ptr<rockseis::Image3D<T>> grad)
{
    int nx, ny, nz;
    nx = grad->getNx();
    ny = grad->getNy();
    nz = grad->getNz();
	if(!grad->getAllocated()) grad->allocateImage();
    grad->read();
    T * data = grad->getImagedata();
    T *thres = (T *) calloc(nx*ny*nz, sizeof(T));
    for (int i=0; i<nx*ny*nz; i++){
        thres[i] = ABS(data[i]);
    }
    std::sort(thres, thres+(nx*ny*nz)); 
    int pos = (int) (97.0*(nx*ny*nz)/100.0);
    T pclip = thres[pos];

    for (int i=0; i<nx*ny*nz; i++){
        if(ABS(data[i]) >= pclip){
            data[i] = 0.0;
        }
    }
    grad->write();
    free(thres);
}


template<typename T>
void TomoAcoustic3D<T>::runGrad() {
    MPImodeling *mpi = this->getMpi();
    std::shared_ptr<rockseis::Data3D<T>> Tobs3D;
    std::shared_ptr<rockseis::Data3D<T>> Tmod3D;
    std::shared_ptr<rockseis::Data3D<T>> Tres3D;
    std::shared_ptr<rockseis::Data3D<T>> weight3D;
    std::shared_ptr<rockseis::Image3D<T>> vpgrad;
    std::shared_ptr<rockseis::ModelEikonal3D<T>> lmodel;

    // Create a sort class
    std::shared_ptr<rockseis::Sort<T>> Sort (new rockseis::Sort<T>());
    Sort->setDatafile(Trecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelEikonal3D<T>> gmodel (new rockseis::ModelEikonal3D<T>(Vpfile, this->getLpml()));

     // Test for problematic model sampling
    if(gmodel->getDx() != gmodel->getDy() || gmodel->getDx() != gmodel->getDz()){
        rs_error("Input model has different dx, dy and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dy = dz).");
    }

    // Create a file to output data misfit values
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());


	if(mpi->getRank() == 0) {
		// Master

        // Get shot map
        Sort->readKeymap();
        Sort->readSortmap();
        size_t ngathers =  Sort->getNensemb();

        // Create a data class for the recorded data
        std::shared_ptr<rockseis::Data3D<T>> Tobs3D (new rockseis::Data3D<T>(Trecordfile));
        // Create modelling and residual data files
        Tmod3D = std::make_shared<rockseis::Data3D<T>>(1, Tobs3D->getNt(), Tobs3D->getDt(), Tobs3D->getOt());
        Tmod3D->setFile(Tmodelledfile);
        Tmod3D->createEmpty(Tobs3D->getNtrace());

        Tres3D = std::make_shared<rockseis::Data3D<T>>(1, Tobs3D->getNt(), Tobs3D->getDt(), Tobs3D->getOt());
        Tres3D->setFile(Tresidualfile);
        Tres3D->createEmpty(Tobs3D->getNtrace());

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

        // Stack gradient
        vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile, gmodel, 1, 1, 1);
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
        std::shared_ptr<rockseis::FatAcoustic3D<T>> fat;
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
                Tobs3D = Sort->get3DGather(work.id);
                size_t ntr = Tobs3D->getNtrace();

                lmodel = gmodel->getLocal(Tobs3D, -1.0*gmodel->getDx(), -1.0*gmodel->getDx(), SMAP);
                lmodel->Expand();

                //Create source 
                std::shared_ptr<rockseis::Data3D<T>> source (new rockseis::Data3D<T>(ntr, 1, 1.0, 0.0));
                source->copyCoords(Tobs3D);
                source->makeMap(lmodel->getGeom(), SMAP);


                // Create first arrival tomography object
                 fat = std::make_shared<rockseis::FatAcoustic3D<T>>(lmodel, source, Tobs3D);
                // Create modelled and residual data objects 
                Tmod3D = std::make_shared<rockseis::Data3D<T>>(ntr, Tobs3D->getNt(), Tobs3D->getDt(), 0.0);
                Tmod3D->copyCoords(Tobs3D);
                Tmod3D->makeMap(lmodel->getGeom());
                fat->setTmod(Tmod3D);
                Tres3D = std::make_shared<rockseis::Data3D<T>>(ntr, Tobs3D->getNt(), Tobs3D->getDt(), 0.0);
                Tres3D->copyCoords(Tobs3D);
                Tres3D->makeMap(lmodel->getGeom());
                fat->setTres(Tres3D);

                if(dataweight){
                    Sort->setDatafile(Dataweightfile);
                    weight3D = Sort->get3DGather(work.id);
                    fat->setTweight(weight3D);
                }

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image3D<T>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, 1, 1, 1);
                fat->setVpgrad(vpgrad);

                // Set logfile
                fat->setLogfile("log.txt-" + std::to_string(work.id));

                /////// Run eikonal solver
                fat->run();

                // Output misfit
                Fmisfit->append(Misfitfile);
                T val = fat->getMisfit();
                Fmisfit->write(&val, 1, work.id*sizeof(T));
                Fmisfit->close();

                // Output modelled and residual data
                Tmod3D->setFile(Tmodelledfile);
                Sort->put3DGather(Tmod3D, work.id);
                Tres3D->setFile(Tresidualfile);
                Sort->put3DGather(Tres3D, work.id);

                // Stack gradient
                vpgrad->write();
                
                // Reset all classes
                Tobs3D.reset();
                Tmod3D.reset();
                Tres3D.reset();
                if(dataweight){
                    weight3D.reset();
                }
                source.reset();
                lmodel.reset();
                vpgrad.reset();
                fat.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi->sendResult(work);		
            }
        }
    }
}

template<typename T>
void TomoAcoustic3D<T>::runBsproj() {
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
	std::shared_ptr<rockseis::ModelEikonal3D<T>> grad (new rockseis::ModelEikonal3D<T>(vpgradfile, this->getLpml()));

	// Read model
	grad->readVelocity();
	
	T *vpgrad;
	vpgrad = grad->getVelocity();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl3D<T>> spline (new rockseis::Bspl3D<T>(grad->getNx(), grad->getNy(), grad->getNz(), grad->getDx(), grad->getDy(), grad->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3));
    int nc = spline->getNc();

	/* Allocating projection arrays */
	float *vpproj= (float *) calloc(nc, sizeof(float));
	if(vpproj==NULL){
		rs_error("TomoAcoustic3D<T>::runBsproj(): Not enough memory to allocate projection array (vpproj)");
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
			rs_error("TomoAcoustic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
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
				for(long int i=0; i<grad->getNx()*grad->getNy()*grad->getNz(); i++){
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
			rs_error("TomoAcoustic3D<T>::runBsproj(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
        MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

	   }
    // Free variables
    free(global_stack);
    free(vpproj);
}

template<typename T>
int TomoAcoustic3D<T>::setInitial(double *x, std::string vpfile)
{
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> bound;
    std::shared_ptr<rockseis::ModelEikonal3D<T>> model_in (new rockseis::ModelEikonal3D<T>(vpfile, 1));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;

    // Read initial model
    model_in->readVelocity();  

    // Check bounds 
    if(this->getConstrain()){
       bound = std::make_shared <rockseis::ModelAcoustic3D<T>>(Lboundfile, Uboundfile, 1 ,0);
       T *lbounddata;
       T *ubounddata;
       T *vp0;
       bound->readModel();
       long Nbound = (bound->getGeom())->getNtot();
       long N = (long) (model_in->getGeom())->getNtot();
       if(N != Nbound) rs_error("TomoAcoustic2D<T>::setInitial(): Geometry in Boundary files does not match geometry in the model.");
       lbounddata = bound->getVp();
       ubounddata = bound->getR();
       vp0 = model_in->getVelocity();
       for (long i=0; i< N; i++){
          if(lbounddata[i] == ubounddata[i]) rs_error("TomoAcoustic2D<T>::setInitial():The lower bound cannot be equal to the upper bound, use a mute function insted.");
          if(vp0[i] <= lbounddata[i]) vp0[i] = lbounddata[i] + 1e-2;
          if(vp0[i] >= ubounddata[i]) vp0[i] = ubounddata[i] - 1e-2;
       }
    }

    // Write initial model files
    model_in->setVelocityfile(VP0FILE);
    model_in->writeVelocity();
    // Write linesearch model files
    model_in->setVelocityfile(VPLSFILE);
    model_in->writeVelocity();
    int N=0;
    switch(this->getParamtype()){
        case PAR_GRID:
            N = (model_in->getGeom())->getNtot();
            break;
        case PAR_BSPLINE:
            spline = std::make_shared<rockseis::Bspl3D<T>>(model_in->getNx(), model_in->getNy(), model_in->getNz(), model_in->getDx(), model_in->getDy(), model_in->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3);
            N=spline->getNc();
            break;
        default:
            rs_error("TomoAcoustic3D<T>::setInitial(): Unknown parameterisation."); 
            break;
    }

    return N;
}

template<typename T>
void TomoAcoustic3D<T>::saveResults(int iter)
{
    std::string name;
    std::string dir;
    dir = RESULTDIR;
    // Write out new models
    std::shared_ptr<rockseis::ModelEikonal3D<T>> lsmodel (new rockseis::ModelEikonal3D<T>(VPLSFILE, 1));
    lsmodel->readVelocity();
    name = dir + "/" + VP_UP + "-" + std::to_string(iter);
    lsmodel->setVelocityfile(name);
    lsmodel->writeVelocity();
}

template<typename T>
void TomoAcoustic3D<T>::saveLinesearch(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelEikonal3D<T>> model0 (new rockseis::ModelEikonal3D<T>(VP0FILE, 1));
    std::shared_ptr<rockseis::ModelEikonal3D<T>> lsmodel (new rockseis::ModelEikonal3D<T>(VPLSFILE, 1));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
    std::shared_ptr<rockseis::ModelEikonal3D<T>> mute;
    std::shared_ptr<rockseis::ModelAcoustic3D<T>> bound;

    // Write linesearch model
    model0->readVelocity();
    lsmodel->readVelocity();
    T *vp0, *vpls;
    T *c, *mod;
    T *vpmutedata;
    T *lbounddata;
    T *ubounddata;
    T *x0;
    double log_in, log_out, exp_in, exp_out;
    vp0 = model0->getVelocity(); 
    vpls = lsmodel->getVelocity(); 
    int i;
    int N, Nmod;

    // If mute
    if(!Mutefile.empty()){
        mute = std::make_shared <rockseis::ModelEikonal3D<T>>(Mutefile, 1);
        long Nmute = (mute->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nmute) rs_error("TomoAcoustic3D<T>::saveLinesearch(): Geometry in Mutefile does not match geometry in the model.");
        mute->readVelocity();
        vpmutedata = mute->getVelocity();
        N = (lsmodel->getGeom())->getNtot();
    }else{
        N = (lsmodel->getGeom())->getNtot();
        vpmutedata = (T *) calloc(N, sizeof(T)); 
        for(i=0; i < N; i++){
            vpmutedata[i] = 1.0;
        }
    }

    if(this->getConstrain()){
        bound = std::make_shared <rockseis::ModelAcoustic3D<T>>(Lboundfile, Uboundfile, 1 ,0);
        bound->readModel();
        long Nbound = (bound->getGeom())->getNtot();
        N = (lsmodel->getGeom())->getNtot();
        if(N != Nbound) rs_error("TomoAcoustic3D<T>::saveLinesearch(): Geometry in Boundary files does not match geometry in the model.");
        lbounddata = bound->getVp();
        ubounddata = bound->getR();
        x0 = (T *) calloc(N, sizeof(T));
    }

    switch (this->getParamtype()){
        case PAR_GRID:
            N = (lsmodel->getGeom())->getNtot();
            for(i=0; i< N; i++)
            {
                if(!this->getConstrain()){
                    vpls[i] = vp0[i] + x[i]*vpmutedata[i]*kvp;
                }else{
                    log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
                    log_out = log(log_in);
                    x0[i] =((T) log_out);
                    exp_in = (double) (-(x[i]*vpmutedata[i]*kvp + x0[i]));
                    exp_out = exp(exp_in);
                    vpls[i] = lbounddata[i] + (ubounddata[i]-lbounddata[i])*(1.0/(1.0 + (T) exp_out));
                }
            }
            lsmodel->writeVelocity();
            break;
        case PAR_BSPLINE:
            Nmod = (lsmodel->getGeom())->getNtot();
            spline = std::make_shared<rockseis::Bspl3D<T>>(model0->getNx(), model0->getNy(), model0->getNz(), model0->getDx(), model0->getDy(), model0->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3);
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
                if(!this->getConstrain()){
                    vpls[i] = vp0[i] + mod[i]*vpmutedata[i]*kvp;
                }else{
                    log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
                    log_out = log(log_in);
                    x0[i] =((T) log_out);
                    exp_in = (double) (-(mod[i]*vpmutedata[i]*kvp + x0[i]));
                    exp_out = exp(exp_in);
                    vpls[i] = lbounddata[i] + (ubounddata[i]-lbounddata[i])*(1.0/(1.0 + (T) exp_out));
                }
            }
            lsmodel->writeVelocity();
            break;
        default:
            rs_error("TomoAcoustic3D<T>::saveLinesearch(): Unknown parameterisation."); 
            break;
    }
    if(Mutefile.empty()){
        free(vpmutedata);
    }

    if(this->getConstrain()){
       free(x0);
    }
}

template<typename T>
void TomoAcoustic3D<T>::readMisfit(double *f)
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
void TomoAcoustic3D<T>::readGrad(double *g)
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

    std::shared_ptr<rockseis::ModelEikonal3D<T>> modelgrad (new rockseis::ModelEikonal3D<T>(vpgradfile, 1));
    std::shared_ptr<rockseis::Bspl3D<T>> spline;
    std::shared_ptr<rockseis::File> Fgrad;
    switch (this->getParamtype()){
        case PAR_GRID:
            modelgrad->readVelocity();
            N = (modelgrad->getGeom())->getNtot();
            gvp = modelgrad->getVelocity(); 
            for(i=0; i< N; i++)
            {
                g[i] = gvp[i]*kvp;
            }
            break;
        case PAR_BSPLINE:
           spline = std::make_shared<rockseis::Bspl3D<T>>(modelgrad->getNx(), modelgrad->getNy(), modelgrad->getNz(), modelgrad->getDx(), modelgrad->getDy(), modelgrad->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3);
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
            rs_error("TomoAcoustic3D<T>::readGrad(): Unknown parameterisation."); 
            break;
    }
}

template<typename T>
void TomoAcoustic3D<T>::combineGradients()
{
    // Gradients
    std::shared_ptr<rockseis::ModelEikonal3D<T>> grad;
    std::shared_ptr<rockseis::ModelEikonal3D<T>> reggrad;
    grad = std::make_shared<rockseis::ModelEikonal3D<T>>(VPGRADFILE, 1);
    reggrad = std::make_shared<rockseis::ModelEikonal3D<T>>(VPREGGRADFILE, 1);

    // Read gradients
    grad->readVelocity();
    reggrad->readVelocity();
    T *vp, *regvp;
    vp = grad->getVelocity(); 
    regvp = reggrad->getVelocity(); 
    int i;
    int N;

    N = (grad->getGeom())->getNtot();
    // Compute 
    for(i=0; i< N; i++)
    {
        vp[i] = vp[i] + reg_alpha[0]*regvp[i];
    }
    grad->setVelocityfile(VPGRADCOMBFILE);
    grad->writeVelocity();
}

template<typename T>
void TomoAcoustic3D<T>::applyChainrule(double *x)
{
   // Models
   std::shared_ptr<rockseis::ModelEikonal3D<T>> model0 (new rockseis::ModelEikonal3D<T>(VP0FILE, 1));
   std::shared_ptr<rockseis::ModelEikonal3D<T>> grad (new rockseis::ModelEikonal3D<T>(VPGRADCOMBFILE, 1));
   std::shared_ptr<rockseis::Bspl3D<T>> spline;
   std::shared_ptr<rockseis::ModelEikonal3D<T>> mute;
   std::shared_ptr<rockseis::ModelAcoustic3D<T>> bound;

   // Write linesearch model
   model0->readVelocity();
   grad->readVelocity();
   T *vp0, *vpgrad;
   T *c, *mod;
   T *vpmutedata;
   T *lbounddata;
   T *ubounddata;
   T *x0;
   double log_in, log_out, exp_in, exp_out;

   vp0 = model0->getVelocity(); 
   vpgrad = grad->getVelocity(); 
   int i;
   int N, Nmod;

   // If mute
   if(!Mutefile.empty()){
      mute = std::make_shared <rockseis::ModelEikonal3D<T>>(Mutefile, 1);
      long Nmute = (mute->getGeom())->getNtot();
      N = (grad->getGeom())->getNtot();
      if(N != Nmute) rs_error("TomoAcoustic3D<T>::applyChainrule(): Geometry in Mutefile does not match geometry in the model.");
      mute->readVelocity();
      vpmutedata = mute->getVelocity();
   }else{
      N = (grad->getGeom())->getNtot();
      vpmutedata = (T *) calloc(N, sizeof(T)); 
      for(i=0; i < N; i++){
         vpmutedata[i] = 1.0;
      }
   }
   bound = std::make_shared <rockseis::ModelAcoustic3D<T>>(Lboundfile, Uboundfile, 1 ,0);
   bound->readModel();
   long Nbound = (bound->getGeom())->getNtot();
   N = (grad->getGeom())->getNtot();
   if(N != Nbound) rs_error("TomoAcoustic3D<T>::applyChainrule(): Geometry in Boundary files does not match geometry in the model.");
   lbounddata = bound->getVp();
   ubounddata = bound->getR();
   x0 = (T *) calloc(N, sizeof(T));
   switch (this->getParamtype()){
      case PAR_GRID:
         N = (grad->getGeom())->getNtot();
         for(i=0; i< N; i++)
         {
            log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
            log_out = log(log_in);
            x0[i] =((T) log_out);
            exp_in = (double) (-(x[i]*vpmutedata[i]*kvp + x0[i]));
            exp_out = exp(exp_in);
            vpgrad[i] *= kvp*(ubounddata[i]-lbounddata[i])*((T) exp_out/((1.0 + (T) exp_out)*(1.0 + (T) exp_out)));
         }
         grad->writeVelocity();
         break;
      case PAR_BSPLINE:
         Nmod = (grad->getGeom())->getNtot();
         spline = std::make_shared<rockseis::Bspl3D<T>>(model0->getNx(), model0->getNy(), model0->getNz(), model0->getDx(), model0->getDy(), model0->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3);
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
            log_in = (double) ((vp0[i] - lbounddata[i])/(ubounddata[i] - vp0[i]));
            log_out = log(log_in);
            x0[i] =((T) log_out);
            exp_in = (double) (-(mod[i]*vpmutedata[i]*kvp + x0[i]));
            exp_out = exp(exp_in);
            vpgrad[i] *= kvp*(ubounddata[i]-lbounddata[i])*((T) exp_out/((1.0 + (T) exp_out)*(1.0 + (T) exp_out)));
         }
         grad->writeVelocity();
         break;
      default:
         rs_error("TomoAcoustic3D<T>::applyChainrule(): Unknown parameterisation."); 
         break;
   }
   if(Mutefile.empty()){
      free(vpmutedata);
   }
   free(x0);
}


template<typename T>
void TomoAcoustic3D<T>::applyMute()
{
    if(!Mutefile.empty()){
        // Models
        std::shared_ptr<rockseis::ModelEikonal3D<T>> model;
        model = std::make_shared<rockseis::ModelEikonal3D<T>>(VPGRADCOMBFILE, 1);
        // Mute
        std::shared_ptr<rockseis::ModelEikonal3D<T>> mute (new rockseis::ModelEikonal3D<T>(Mutefile, 1));

        // Mute model and write
        model->readVelocity();
        mute->readVelocity();
        T *vp, *vpmute;
        vp = model->getVelocity(); 
        vpmute = mute->getVelocity(); 
        int i;
        int N;

        N = (model->getGeom())->getNtot();
        for(i=0; i< N; i++)
        {
            vp[i] = vp[i]*vpmute[i];
        }
        model->setVelocityfile(VPGRADMUTEFILE);
        model->writeVelocity();
    }
}

template<typename T>
void TomoAcoustic3D<T>::computeRegularisation(double *x)
{
    // Models
    std::shared_ptr<rockseis::ModelEikonal3D<T>> model (new rockseis::ModelEikonal3D<T>(VPLSFILE, 1));
    std::shared_ptr<Der<double>> der (new Der<double>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), 8));
    std::shared_ptr<rockseis::Bspl3D<double>> spline;

    // Write linesearch model
    model->readVelocity();
    double *dvpdx, *dvpdy, *dvpdz;
    T *vpgrad;
    double *gwrk;
    double *c, *mod;
    double *df = der->getDf();
    int i;
    int N, Nmod;

    Nmod = (model->getGeom())->getNtot();
    dvpdx = (double *) calloc(Nmod, sizeof(double));
    dvpdy = (double *) calloc(Nmod, sizeof(double));
    dvpdz = (double *) calloc(Nmod, sizeof(double));
    gwrk = (double *) calloc(Nmod, sizeof(double));
    model->readVelocity();
    model->setVelocityfile(VPREGGRADFILE);
    vpgrad = model->getVelocity();
    switch (this->getParamtype()){
        case PAR_GRID:
            N = (model->getGeom())->getNtot();
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
            
            break;
        case PAR_BSPLINE:
            Nmod = (model->getGeom())->getNtot();
            spline = std::make_shared<rockseis::Bspl3D<double>>(model->getNx(), model->getNy(), model->getNz(), model->getDx(), model->getDy(), model->getDz(), this->getDtx(),this->getDty(), this->getDtz(), 3, 3, 3);
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
            der->ddy_fw(mod);
            for(i=0; i< N; i++)
            {
                dvpdy[i] = df[i]*kvp;
            }
            der->ddz_fw(mod);
            for(i=0; i< Nmod; i++)
            {
                dvpdz[i] = df[i]*kvp;
            }

            break;
        default:
            rs_error("TomoAcoustic3D<T>::saveLinesearch(): Unknown parameterisation."); 
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
        M = dvpdx[i]*dvpdx[i]  + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
        M = sqrt(M);
        gwrk[i] = dvpdy[i]/M;
    }
    der->ddy_bw(gwrk);

    for(i=0; i< Nmod; i++)
    {
        vpgrad[i] -= df[i];
    }

    for(i=0; i< Nmod; i++)
    {
        M = dvpdx[i]*dvpdx[i]  + dvpdy[i]*dvpdy[i] + dvpdz[i]*dvpdz[i] + reg_eps[0]*reg_eps[0];
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
    model->writeVelocity();

    // Free variables
    free(dvpdx);
    free(dvpdy);
    free(dvpdz);
    free(gwrk);
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Tomo<float>;
template class Tomo<double>;

template class TomoAcoustic2D<float>;
template class TomoAcoustic2D<double>;

template class TomoAcoustic3D<float>;
template class TomoAcoustic3D<double>;

}
