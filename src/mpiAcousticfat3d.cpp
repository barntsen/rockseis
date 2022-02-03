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
            PRINT_DOC(# MPI 3d acoustic modelling default configuration file);
            PRINT_DOC();
            PRINT_DOC(# Booleans);
            PRINT_DOC(        TTrecord = "true";  # Set these to true if recording or snapshoting is to be made.);
            PRINT_DOC(        TTsnap = "false";);
            PRINT_DOC();
            PRINT_DOC(# Input files);
            PRINT_DOC(        Vp = "Vp3d.rss";);
            PRINT_DOC(        Survey = "3DSurvey.rss";);
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
    std::shared_ptr<rockseis::Data3D<float>> TTdata3D;

    // Create a local model class
	std::shared_ptr<rockseis::ModelEikonal3D<float>> lmodel;

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
	std::shared_ptr<rockseis::ModelEikonal3D<float>> gmodel (new rockseis::ModelEikonal3D<float>(Vpfile, 10));

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
        std::shared_ptr<rockseis::Data3D<float>> Shotgeom;
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

                Shotgeom = Sort->get3DGather(work.id);
                size_t ntr = Shotgeom->getNtrace();
                lmodel = gmodel->getLocal(Shotgeom, -1.0*gmodel->getDx(),-1.0*gmodel->getDy(), SMAP);
                lmodel->Expand();

                // Set shot coordinates and make a map
	            std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(1, 1, 1.0, 0.0));
                source->copyCoords(Shotgeom);
                source->makeMap(lmodel->getGeom(), SMAP);

                                // Run modelling 
                std::shared_ptr<rockseis::RaysAcoustic3D<float>> rays (new rockseis::RaysAcoustic3D<float>(lmodel));

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();


                // Output record
                if(TTrecord){
                    TTdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                    // Copy geometry to Data
                    TTdata3D->copyCoords(Shotgeom);
                    TTdata3D->makeMap(lmodel->getGeom());
                    rays->recordData(TTdata3D, GMAP);

                    TTdata3D->setFile(TTrecordfile);
                    Sort->put3DGather(TTdata3D, work.id);
                    TTdata3D.reset();
                }

                // Create snapshots
                if(TTsnap){ 
                    std::shared_ptr<rockseis::Snapshot3D<float>> Ttimes;
                    Ttimes = std::make_shared<rockseis::Snapshot3D<float>>(lmodel, 1);
                    Ttimes->openSnap(TTsnapfile + "-" + std::to_string(work.id), 'w'); // Create a new snapshot file
                    Ttimes->setData(rays->getTT(), 0); //Set field to snap
                    Ttimes->writeSnap(0);
                    Ttimes.reset();
                }
                
                // Reset all classes
                source.reset();
                rays.reset();
                Shotgeom.reset();
                lmodel.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

