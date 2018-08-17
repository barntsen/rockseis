#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "rays.h"
#include "utils.h"
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
            PRINT_DOC(# MPI 2d acoustic modelling default configuration file);
            PRINT_DOC();
            PRINT_DOC(# Booleans);
            PRINT_DOC(        TTrecord = "true";  # Set these to true if recording or snapshoting is to be made.);
            PRINT_DOC(        TTsnap = "false";);
            PRINT_DOC();
            PRINT_DOC(# Input files);
            PRINT_DOC(        Vp = "Vp2d.rss";);
            PRINT_DOC(        Survey = "2DSurvey.rss";);
            PRINT_DOC();
            PRINT_DOC(# Output files);
            PRINT_DOC(        TTrecordfile = "TTrecord.rss";);
            PRINT_DOC(        TTsnapfile = "TTsnaps.rss";);
        }
        exit(1);
    }
    bool status;
	/* General input parameters */
    std::string Surveyfile;
    std::string Vpfile;

    bool TTsnap=0, TTrecord=0;
    std::string TTsnapfile;
    std::string TTrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> TTdata2D;

    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> lmodel;

    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("TTrecord", &TTrecord) == INPARSE_ERR) status = true;
    if(TTrecord){
        if(Inpar->getPar("TTrecordfile", &TTrecordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("TTsnap", &TTsnap) == INPARSE_ERR) status = true;
    if(TTsnap){
        if(Inpar->getPar("TTsnapfile", &TTsnapfile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Surveyfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> gmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Vpfile, 0 ,false));

    // Test for problematic model sampling
    if(gmodel->getDx() != gmodel->getDz()){
        rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dz).");
    }

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Surveyfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort->getNensemb();
        
        if(TTrecord){
            // create an empty data file
            Sort->createEmptydataset(TTrecordfile, 1, 1.0, 0.0);
        }

		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();
	
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::Data2D<float>> Shotgeom;
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

                Shotgeom = Sort->get2DGather(work.id);
                size_t ntr = Shotgeom->getNtrace();
                lmodel = gmodel->getLocal(Shotgeom, -6.0*gmodel->getDx(), SMAP);

                // Set shot coordinates and make a map
	            std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(ntr, 1, 1.0, 0.0));
                source->copyCoords(Shotgeom);
                source->makeMap(lmodel->getGeom(), SMAP);

                                // Run modelling 
                std::shared_ptr<rockseis::RaysAcoustic2D<float>> rays (new rockseis::RaysAcoustic2D<float>(lmodel));

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();


                // Output record
                if(TTrecord){
                    TTdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                    // Copy geometry to Data
                    TTdata2D->copyCoords(Shotgeom);
                    TTdata2D->makeMap(lmodel->getGeom());
                    rays->recordData(TTdata2D, GMAP);

                    TTdata2D->setFile(TTrecordfile);
                    Sort->put2DGather(TTdata2D, work.id);
                }

                // Create snapshots
                if(TTsnap){ 
                    std::shared_ptr<WavesAcoustic2D<float>> waves (new WavesAcoustic2D<float>(lmodel, 1, 1.0, 0.0));
                    std::shared_ptr<rockseis::Snapshot2D<float>> Ttimes;
                    Ttimes = std::make_shared<rockseis::Snapshot2D<float>>(waves, 1);
                    Ttimes->openSnap(TTsnapfile + "-" + std::to_string(work.id), 'w'); // Create a new snapshot file
                    Ttimes->setData(rays->getTT(), 0); //Set field to snap
                    Ttimes->writeSnap(0);
                }
                
                // Reset all classes
                Shotgeom.reset();
                lmodel.reset();
                TTdata2D.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

