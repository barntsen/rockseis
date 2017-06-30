#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "rtm.h"
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
			PRINT_DOC(# MPI 2d acoustic reverse-time migration configuration file);
			PRINT_DOC();
			PRINT_DOC(# Modelling parameters);
			PRINT_DOC(freesurface = "false";  # True if free surface should be on);
			PRINT_DOC(order = "8";  # Order of finite difference stencil);
			PRINT_DOC(lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
			PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
			PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Checkpointing parameters);
			PRINT_DOC(snapmethod = "0";  );
			PRINT_DOC(nsnaps = "11";);
			PRINT_DOC(incore = "true";  );
			PRINT_DOC();
			PRINT_DOC(# Migration parameters);
			PRINT_DOC(nhx = "1";);
			PRINT_DOC(nhz = "1";);
			PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp2d.rss";);
			PRINT_DOC(Rho = "Rho2d.rss";);
			PRINT_DOC(Wavelet = "Wav2d.rss";);
			PRINT_DOC(Precordfile = "Pshots2d.rss";);
			PRINT_DOC(Pimagefile = "Pimage2d.rss";);
			PRINT_DOC(Psnapfile = "Psnaps2d.rss";);
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
    float apertx;
    int nhx=1, nhz=1;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Pimagefile;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data2D<float>> shot2D;
    std::shared_ptr<rockseis::Data2D<float>> shot2Di;
    std::shared_ptr<rockseis::ImageAcoustic2D<float>> pimage;

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
    if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhz", &nhz) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapmethod", &snapmethod) == INPARSE_ERR) status = true;
    rockseis::rs_snapmethod checkpoint = static_cast<rockseis::rs_snapmethod>(snapmethod);
    switch(checkpoint){
        case rockseis::FULL:
            break;
        case rockseis::OPTIMAL:
            if(Inpar->getPar("nsnaps", &nsnaps) == INPARSE_ERR) status = true;
            if(Inpar->getPar("incore", &incore) == INPARSE_ERR) status = true;
            break;
        case rockseis::EDGES:
            break;
        default:
            rockseis::rs_error("Invalid option of snapshot saving (snapmethod)."); 
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> gmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Rhofile, lpml ,fs));
    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> lmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort->getNensemb();
        
		// Create work queue
		for(unsigned long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,{'\0'}});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();

        // Image
        pimage = std::make_shared<rockseis::ImageAcoustic2D<float>>(Pimagefile, gmodel, nhx, nhz);
        pimage->createEmpty();

		for(unsigned long int i=0; i<ngathers; i++) {
            pimage->stackImage(Pimagefile + "-" + std::to_string(i));
        }
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::RtmAcoustic2D<float>> rtm;
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
                size_t ntr = shot2D->getNtrace();

                lmodel = gmodel->getLocal(shot2D, apertx, SMAP);
                pimage = std::make_shared<rockseis::ImageAcoustic2D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhz);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(shot2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                // Interpolate shot
                shot2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(shot2D, shot2Di);
                shot2Di->makeMap(lmodel->getGeom(), GMAP);

                rtm = std::make_shared<rockseis::RtmAcoustic2D<float>>(lmodel, pimage, source, shot2Di, order, snapinc);
                // Setting Snapshot file 
                rtm->setSnapfile(Psnapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                rtm->setNcheck(nsnaps);
                rtm->setIncore(incore);

                // Set logfile
                rtm->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

                switch(checkpoint){
                    case rockseis::FULL:
                        rtm->run();
                        break;
                    case rockseis::OPTIMAL:
                        rtm->run_optimal();
                        break;
                    case rockseis::EDGES:
                        rtm->run_edge();
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
                mpi.sendResult(work);		
            }
        }
    }
}

