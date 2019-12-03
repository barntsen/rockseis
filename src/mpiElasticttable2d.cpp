#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "rays.h"
#include "ttable.h"
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
            PRINT_DOC(# Padding);
            PRINT_DOC(        Lpml = "3";);
            PRINT_DOC();
            PRINT_DOC(# Sampling);
            PRINT_DOC(        Souinc = "1";);
            PRINT_DOC(        Recinc = "1";);
            PRINT_DOC();
            PRINT_DOC(# Input files);
            PRINT_DOC(        Vp = "Vp2d.rss";);
            PRINT_DOC(        Vs = "Vs2d.rss";);
            PRINT_DOC(        Survey = "2DSurvey.rss";);
            PRINT_DOC();
            PRINT_DOC(# Output files);
            PRINT_DOC(        Sou_ttablefile = "Sou_ttable2d.rss";);
            PRINT_DOC(        Rec_ttablefile = "Rec_ttable2d.rss";);
        }
        exit(1);
    }
    bool status;
	/* General input parameters */
    int Lpml;
    int souinc, recinc;
    int nsoufin, nrecfin;
    std::string Surveyfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Sou_ttablefile;
    std::string Rec_ttablefile;
    std::shared_ptr<rockseis::Data2D<float>> source;
    std::shared_ptr<rockseis::RaysAcoustic2D<float>> rays;
    std::shared_ptr<rockseis::Ttable2D<float>> Sou_ttable;
    std::shared_ptr<rockseis::Ttable2D<float>> Rec_ttable;


    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Lpml", &Lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Souinc", &souinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Recinc", &recinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Sou_ttablefile", &Sou_ttablefile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rec_ttablefile", &Rec_ttablefile) == INPARSE_ERR) status = true;

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Surveyfile);

    // Create a global model class
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vpgmodel (new rockseis::ModelEikonal2D<float>(Vpfile, Lpml));
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vsgmodel (new rockseis::ModelEikonal2D<float>(Vsfile, Lpml));

    // Test for problematic model sampling
    if(vpgmodel->getDx() != vpgmodel->getDz()){
        rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dz).");
    }

    // Test for the compatibility of the models
    if((vpgmodel->getGeom())->compare(vsgmodel->getGeom()) == true){
        rs_error("Geometries of Vp and Vs model are not the same.");
    }

    // Read and expand global model
    vpgmodel->readVelocity();
    vpgmodel->Expand();

    vsgmodel->readVelocity();
    vsgmodel->Expand();

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

	if(mpi.getRank() == 0) {
		// Master
    
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

        // Create a travel time table class
        Sou_ttable = std::make_shared<rockseis::Ttable2D<float>> (vpgmodel, nsoufin);
        Sou_ttable->setFilename(Sou_ttablefile);
        Sou_ttable->createEmptyttable();

        Rec_ttable = std::make_shared<rockseis::Ttable2D<float>> (vsgmodel, nrecfin);
        Rec_ttable->setFilename(Rec_ttablefile);
        Rec_ttable->createEmptyttable();
        
        /******************             Creating source side     ***********************/
		// Create work queue
		for(long int i=0; i<nsoufin; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();

        // Reset mpi 
        mpi.clearWork();
        /******************             Creating receiver side     ***********************/

        // Create new list of positions
        Sort->createReceivermap(Surveyfile); 
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

    }
    else {
        /* Slave */
        size_t number;
        std::shared_ptr<rockseis::Data2D<float>> Shotgeom;

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
                // Do some work
                Sort->readKeymap();
                Sort->readSortmap();

                number = work.id*souinc;
                if (number > Sort->getNensemb()-1) number = Sort->getNensemb()-1;
                Shotgeom = Sort->get2DGather(number);

                // Set shot coordinates and make a map
	            source = std::make_shared<rockseis::Data2D<float>>(1, 1, 1.0, 0.0);
                source->copyCoords(Shotgeom);

                source->makeMap(vpgmodel->getGeom(), SMAP);

                // Run modelling 
                rays = std::make_shared<rockseis::RaysAcoustic2D<float>>(vpgmodel);

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();

                // Create traveltime table
                Sou_ttable = std::make_shared<rockseis::Ttable2D<float>> (Sou_ttablefile);
                Sou_ttable->fetchTtabledata(rays, source, work.id); //Get traveltime data
                Sou_ttable->writeTtable(work.id);
                Sou_ttable.reset();
                
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
	            source = std::make_shared<rockseis::Data2D<float>>(1, 1, 1.0, 0.0);
                source->copyCoords(Shotgeom);

                source->makeMap(vsgmodel->getGeom(), SMAP);

                // Run modelling 
                rays = std::make_shared<rockseis::RaysAcoustic2D<float>>(vsgmodel);

                /* initialize traveltime field at source positions */
                rays->insertSource(source, SMAP);
                rays->solve();

                // Create traveltime table
                Rec_ttable = std::make_shared<rockseis::Ttable2D<float>> (Rec_ttablefile);
                Rec_ttable->fetchTtabledata(rays, source, work.id); //Get traveltime data
                Rec_ttable->writeTtable(work.id);
                Rec_ttable.reset();
                
                // Reset all classes
                source.reset();
                Shotgeom.reset();
                rays.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

