#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "kdmig.h"
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
			PRINT_DOC(lpml = "10"; # Size of pml absorbing boundary (should be larger than order + 5 ));
			PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Migration parameters);
			PRINT_DOC(freqinc = "4"; # Integer frequency interval to sum over);
			PRINT_DOC(nhx = "1";);
			PRINT_DOC(nhz = "1";);
			PRINT_DOC();
			PRINT_DOC(# Booleans);
			PRINT_DOC(Gather = "false"; # If surface gathers are to be output);
			PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp2d.rss";);
			PRINT_DOC(Precordfile = "Pshots2d.rss";);
			PRINT_DOC(Pimagefile = "Pimage2d.rss";);
			PRINT_DOC(Pgatherfile = "Pgather2d.rss";);
			PRINT_DOC();
		}
        exit(1);
    }
    bool status;
	/* General input parameters */
	int lpml;
    float apertx;
    int nhx=1, nhz=1;
	int freqinc;
    std::string Vpfile;
    std::string Pimagefile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data2D<float>> shot2D;
    std::shared_ptr<rockseis::Image2D<float>> pimage;
    bool Gather;
    std::string Pgatherfile;
    std::shared_ptr<rockseis::Data2D<float>> pgather;
	std::shared_ptr<rockseis::ModelEikonal2D<float>> lmodel;

    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freqinc", &freqinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhz", &nhz) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Gather", &Gather) == INPARSE_ERR) status = true;
    if(Gather){
        if(Inpar->getPar("Pgatherfile", &Pgatherfile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Precordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelEikonal2D<float>> gmodel (new rockseis::ModelEikonal2D<float>(Vpfile, lpml));

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort->getNensemb();
        
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();

        // Image gathers
        if(Gather){
            std::shared_ptr<rockseis::File> Fimg (new rockseis::File());
            Fimg->input(Pimagefile + "-" + std::to_string(0));
            pgather = std::make_shared<rockseis::Data2D<float>>(Fimg->getN(1),Fimg->getN(3),Fimg->getD(3),Fimg->getO(3));
            pgather->setFile(Pgatherfile);
            pgather->open("o");
            for(long int i=0; i<ngathers; i++) {
                pgather->putImage(Pimagefile + "-" + std::to_string(i));
            }
            pgather->close();
            Fimg->close();
        }
        // Image
        pimage = std::make_shared<rockseis::Image2D<float>>(Pimagefile, gmodel, nhx, nhz);
        pimage->createEmpty();
		for(long int i=0; i<ngathers; i++) {
            pimage->stackImage(Pimagefile + "-" + std::to_string(i));
            remove_file(Pimagefile + "-" + std::to_string(i));
        }
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::KdmigAcoustic2D<float>> kdmig;
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
                pimage = std::make_shared<rockseis::Image2D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhz);

                // Map coordinates to model
                shot2D->makeMap(lmodel->getGeom(), SMAP);
                shot2D->makeMap(lmodel->getGeom(), GMAP);

                // Create imaging class
                kdmig = std::make_shared<rockseis::KdmigAcoustic2D<float>>(lmodel, shot2D, pimage);

                // Set frequency decimation 
                kdmig->setFreqinc(freqinc);

                // Set logfile
                kdmig->setLogfile("log.txt-" + std::to_string(work.id));

                // Run migration
                kdmig->run();

                // Output image
                pimage->write();
                
                // Reset all classes
                shot2D.reset();
                lmodel.reset();
                pimage.reset();
                kdmig.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}
