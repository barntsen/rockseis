#include "inversion.h"


namespace rockseis {

// =============== ABSTRACT INVERSION CLASS =============== //
template<typename T>
Inversion<T>::Inversion() {
    //Do nothing
}

template<typename T>
Inversion<T>::Inversion(MPImodeling *_mpi) {
    mpi = _mpi;
}

template<typename T>
Inversion<T>::~Inversion() {
    //Do nothing
}

// =============== 2D ACOUSTIC INVERSION CLASS =============== //
//
template<typename T>
InversionAcoustic2D<T>::InversionAcoustic2D() {
    // Do nothing
}

template<typename T>
InversionAcoustic2D<T>::InversionAcoustic2D(MPImodeling *mpi): Inversion<T>(mpi) {
    // Do nothing
}

template<typename T>
InversionAcoustic2D<T>::~InversionAcoustic2D() {
    //Do nothing
}

template<typename T>
void InversionAcoustic2D<T>::runAcousticfwigrad2d() {
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
void InversionAcoustic2D<T>::runBsprojection2d() {
    MPImodeling *mpi = this->getMpi();
	float vpsum = 0.0; // Sum over splines
	float rhosum = 0.0; // Sum over splines
	float *global_stack;
    float *c;
    float *wrk;

	// Get gradients
	std::shared_ptr<rockseis::ModelAcoustic2D<T>> grad (new rockseis::ModelAcoustic2D<T>(Vpgradfile, Rhogradfile, this->getLpml() ,this->getFs()));

	// Read model
	grad->readModel();
	
	T *vpgrad, *rhograd;
	vpgrad = grad->getVp();
	rhograd = grad->getR();

    /* Initializing spline */
    std::shared_ptr<rockseis::Bspl> spline (new rockseis::Bspl(grad->getNx(), 1, grad->getNz(), grad->getDx(), 1.0, grad->getDz(), this->getDtx(), 1.0, this->getDtz(), 3, 2));
    int nc = spline->getNc();

	/* Allocating projection arrays */
	float *vpproj= (float *) calloc(grad->getNx()*grad->getNz(), sizeof(float));
	if(vpproj==NULL){
		rs_error("InversionAcoustic2D<T>::runBsprojection2d(): Not enough memory to allocate projection array (vpproj)");
	}
	float *rhoproj= (float *) calloc(grad->getNx()*grad->getNz(), sizeof(float));
	if(rhoproj==NULL){
		rs_error("InversionAcoustic2D<T>::runBsprojection2d(): Not enough memory to allocate projection array (rhoproj)");
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

		global_stack= (float *) calloc(grad->getNx()*grad->getNz(), sizeof(float));
		if(global_stack==NULL){
			rs_error("InversionAcoustic2D<T>::runBsprojection2d(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
		MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   

		/* Output spline */
        std::shared_ptr<File> Fout (new File());
        Fout->output(VPPROJGRADFILE);
        Fout->setN(1,nc);
        Fout->setData_format(sizeof(float));
        Fout->write(global_stack, nc);
        Fout->close();

		for(long int i=0; i< nc; i++){
			global_stack[i] = 0.0;
		}
		/* Starting reduce operation */
		MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);   

		/* Output spline */
        Fout->output(RHOPROJGRADFILE);
        Fout->setN(1,nc);
        Fout->setData_format(sizeof(float));
        Fout->write(global_stack, nc);
        Fout->close();

       }else {
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
                // Do work
                c = spline->getSpline();
                wrk = spline->getMod();
				c[work.id]=1.0; // Projection point
				spline->bisp(); // Evaluate spline for this coefficient
				vpsum = 0.0;
				rhosum = 0.0;
				for(long int i=0; i<nc; i++){
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

		global_stack= (float *) calloc(grad->getNx()*grad->getNz(), sizeof(float));
		if(global_stack==NULL){
			rs_error("InversionAcoustic2D<T>::runBsprojection2d(): Not enough memory to allocate global stack array");
		}

		/* Starting reduce operation */
		MPI_Reduce(vpproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
		MPI_Reduce(rhoproj, global_stack, nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
	   }
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Inversion<float>;
template class Inversion<double>;

template class InversionAcoustic2D<float>;
template class InversionAcoustic2D<double>;

}


