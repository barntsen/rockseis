#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "fwi.h"
#include "pml.h"
#include "waves.h"
#include "utils.h"
#include "der.h"
#include "sort.h"
#include "data.h"
#include "file.h"
#include "interp.h"
#include "parallel.h"

using namespace rockseis;

int main(int argc, char** argv) {
    // Initializing MPI
    MPImodeling mpi = MPImodeling(&argc,&argv);
    if(mpi.getNrank() < 2){
        rs_error("This is a parallel program, it must run with at least 2 processors, use mpirun.");
    }

    if(argc < 2){
        if(mpi.getRank() == 0){
			PRINT_DOC(# MPI 3d acoustic fullwaveform inversion gradient computation configuration file);
			PRINT_DOC();
			PRINT_DOC(# Modelling parameters);
			PRINT_DOC(freesurface = "false";  # True if free surface should be on);
			PRINT_DOC(order = "8";  # Order of finite difference stencil);
			PRINT_DOC(lpml = "10"; # Size of pml absorbing boundary );
			PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
			PRINT_DOC(apertx = "1000"; # Aperture for local model (source is in the middle));
			PRINT_DOC(aperty = "1000"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Checkpointing parameters);
			PRINT_DOC(snapmethod = "0"; # 0 - fullcheckpointing; 1 - optimal checkpointing );
			PRINT_DOC(misfit_type= "0"; # 0 - Difference; 1 - Correlation );
			PRINT_DOC(nsnaps = "11"; # Number of checkpoints to store);
			PRINT_DOC(incore = "true"; # Do checkpointing in memory);
			PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp3d.rss";);
			PRINT_DOC(Rho = "Rho3d.rss";);
			PRINT_DOC(Wavelet = "Wav3d.rss";);
			PRINT_DOC(Precordfile = "Prec3d.rss"; # Input observed data);
			PRINT_DOC(Pmodelledfile = "Pmod3d.rss"; # File to output modelled data);
			PRINT_DOC(Presidualfile = "Pres3d.rss"; # File to output residuals);
			PRINT_DOC(Vpgradfile = "Vpgrad3d.rss"; # File to output gradient with respect to Vp);
			PRINT_DOC(Rhogradfile = "Rhograd3d.rss"; # File to output gradient with respect to Rho);
			PRINT_DOC(Wavgradfile = "Wavgrad3d.rss"; # File to output gradient with respect to Wav);
			PRINT_DOC(Psnapfile = "Psnaps3d.rss"; # File to output temporary snapshots);
			PRINT_DOC();
		}
        exit(1);
    }
    bool status;
	/* General input parameters */
	int lpml;
	bool fs;
    bool incore = false;
	int order;
	int snapinc;
	int nsnaps = 0;
	int snapmethod;
	int misfit_type;
    float apertx, aperty;
    int nhx=1, nhy=1, nhz=1;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Psnapfile;
    std::string Precordfile;
    std::string Pmodelledfile;
    std::string Presidualfile;
    std::shared_ptr<rockseis::Data3D<float>> shot3D;
    std::shared_ptr<rockseis::Data3D<float>> shot3Di;
    std::shared_ptr<rockseis::Data3D<float>> shotmod3D;
    std::shared_ptr<rockseis::Data3D<float>> shotmod3Di;
    std::shared_ptr<rockseis::Data3D<float>> shotres3D;
    std::shared_ptr<rockseis::Data3D<float>> shotres3Di;
    std::shared_ptr<rockseis::Image3D<float>> vpgrad;
    std::shared_ptr<rockseis::Image3D<float>> rhograd;
    std::shared_ptr<rockseis::Data3D<float>> wavgrad;

    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("order", &order) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vpgradfile", &Vpgradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rhogradfile", &Rhogradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavgradfile", &Wavgradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Presidualfile", &Presidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Pmodelledfile", &Pmodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapmethod", &snapmethod) == INPARSE_ERR) status = true;
    rockseis::rs_snapmethod checkpoint = static_cast<rockseis::rs_snapmethod>(snapmethod);
    switch(checkpoint){
        case rockseis::FULL:
            break;
        case rockseis::OPTIMAL:
            if(Inpar->getPar("nsnaps", &nsnaps) == INPARSE_ERR) status = true;
            if(Inpar->getPar("incore", &incore) == INPARSE_ERR) status = true;
            break;
        default:
            rockseis::rs_error("Invalid option of snapshot saving (snapmethod)."); 
    }
    if(Inpar->getPar("misfit_type", &misfit_type) == INPARSE_ERR) status = true;
    rockseis::rs_fwimisfit fwimisfit = static_cast<rockseis::rs_fwimisfit>(misfit_type);

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> gmodel (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));
    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> lmodel (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort->getNensemb();

        // Wavelet gradient
        wavgrad = std::make_shared<rockseis::Data3D<float>>(1, source->getNt(), source->getDt(), 0.0);
        wavgrad->setFile(Wavgradfile);
        wavgrad->createEmpty(ngathers);

        // Create a data class for the recorded data
        std::shared_ptr<rockseis::Data3D<float>> shot3D (new rockseis::Data3D<float>(Precordfile));
        // Create modelling and residual data files
        shotmod3D = std::make_shared<rockseis::Data3D<float>>(1, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
        shotmod3D->setFile(Pmodelledfile);
        shotmod3D->createEmpty(shot3D->getNtrace());

        shotres3D = std::make_shared<rockseis::Data3D<float>>(1, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
        shotres3D->setFile(Presidualfile);
        shotres3D->createEmpty(shot3D->getNtrace());
        
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();

        // Image
        vpgrad = std::make_shared<rockseis::Image3D<float>>(Vpgradfile, gmodel, nhx, nhy, nhz);
        vpgrad->createEmpty();

        rhograd = std::make_shared<rockseis::Image3D<float>>(Rhogradfile, gmodel, nhx, nhy, nhz);
        rhograd->createEmpty();

		for(long int i=0; i<ngathers; i++) {
            vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            remove_file(Vpgradfile + "-" + std::to_string(i));
            rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
            remove_file(Rhogradfile + "-" + std::to_string(i));
        }
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::FwiAcoustic3D<float>> fwi;
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
                shot3D = Sort->get3DGather(work.id);
                size_t ntr = shot3D->getNtrace();

                lmodel = gmodel->getLocal(shot3D, apertx, aperty, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot3D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot3D, shot3Di);
                shot3Di->makeMap(lmodel->getGeom(), GMAP);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiAcoustic3D<float>>(lmodel, source, shot3Di, order, snapinc);

                // Create modelled and residual data objects 
                shotmod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotmod3D->copyCoords(shot3D);
                shotmod3D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDatamodP(shotmod3D);
                shotres3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                shotres3D->copyCoords(shot3D);
                shotres3D->makeMap(lmodel->getGeom(), GMAP);
                fwi->setDataresP(shotres3D);
                
                // Setting misfit type
                fwi->setMisfit_type(fwimisfit);

                // Creating gradient objects
                vpgrad = std::make_shared<rockseis::Image3D<float>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                rhograd = std::make_shared<rockseis::Image3D<float>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);

                // Setting up gradient objects in fwi class
                fwi->setVpgrad(vpgrad);
                fwi->setRhograd(rhograd);

                wavgrad = std::make_shared<rockseis::Data3D<float>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
                wavgrad->setField(rockseis::PRESSURE);
                // Copy geometry
                wavgrad->copyCoords(source);
                wavgrad->makeMap(lmodel->getGeom(), SMAP);
                fwi->setWavgrad(wavgrad);

                // Setting Snapshot file 
                fwi->setSnapfile(Psnapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                fwi->setNcheck(nsnaps);
                fwi->setIncore(incore);

                // Set logfile
                fwi->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                // Run simulation
                switch(checkpoint){
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

                // Output modelled and residual data
                shotmod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
                shotmod3Di->setFile(Pmodelledfile);
                interp->interp(shotmod3D, shotmod3Di);
                Sort->put3DGather(shotmod3Di, work.id);

                shotres3Di = std::make_shared<rockseis::Data3D<float>>(ntr, shot3D->getNt(), shot3D->getDt(), shot3D->getOt());
                shotres3Di->setFile(Presidualfile);
                interp->interp(shotres3D, shotres3Di);
                Sort->put3DGather(shotres3Di, work.id);

                
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
                fwi.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

