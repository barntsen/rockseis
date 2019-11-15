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
			PRINT_DOC(Vs = "Vs2d.rss";);
			PRINT_DOC(Vxrecordfile = "Pshots2d.rss";);
			PRINT_DOC(Simagefile = "Simage2d.rss";);
			PRINT_DOC(Sgatherfile = "Sgather2d.rss";);
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
    std::string Vsfile;
    std::string Simagefile;
    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> shot2D;
    std::shared_ptr<rockseis::Image2D<float>> simage;
    bool Gather;
    std::string Sgatherfile;
    std::shared_ptr<rockseis::Data2D<float>> sgather;
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> vplmodel;
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> vslmodel;

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
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Simagefile", &Simagefile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhz", &nhz) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Gather", &Gather) == INPARSE_ERR) status = true;
    if(Gather){
        if(Inpar->getPar("Sgatherfile", &Sgatherfile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Vxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> vpgmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Vpfile, lpml ,false));
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> vsgmodel (new rockseis::ModelAcoustic2D<float>(Vsfile, Vsfile, lpml ,false));

    // Test for the compatibility of the models
    if((vpgmodel->getGeom())->compare(vsgmodel->getGeom()) == true){
        rs_error("Geometries of Vp and Vs model are not the same.");
    }

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Vxrecordfile); 
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
            Fimg->input(Simagefile + "-" + std::to_string(0));
            sgather = std::make_shared<rockseis::Data2D<float>>(Fimg->getN(1),Fimg->getN(3),Fimg->getD(3),Fimg->getO(3));
            sgather->setFile(Sgatherfile);
            sgather->open("o");
            for(long int i=0; i<ngathers; i++) {
                sgather->putImage(Simagefile + "-" + std::to_string(i));
            }
            sgather->close();
            Fimg->close();
        }
        // Image
        simage = std::make_shared<rockseis::Image2D<float>>(Simagefile, vpgmodel, nhx, nhz);
        simage->createEmpty();
		for(long int i=0; i<ngathers; i++) {
            simage->stackImage(Simagefile + "-" + std::to_string(i));
            remove_file(Simagefile + "-" + std::to_string(i));
        }
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::KdmigElastic2D<float>> kdmig;
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
                vplmodel = vpgmodel->getLocal(shot2D, apertx, SMAP);
                vplmodel->staggerModels_Eikonal();

                vslmodel = vsgmodel->getLocal(shot2D, apertx, SMAP);
                vslmodel->staggerModels_Eikonal();

                // Make image class
                simage = std::make_shared<rockseis::Image2D<float>>(Simagefile + "-" + std::to_string(work.id), vplmodel, nhx, nhz);

                // Map coordinates to model
                shot2D->makeMap(vplmodel->getGeom(), SMAP);
                shot2D->makeMap(vplmodel->getGeom(), GMAP);

                // Create imaging class
                kdmig = std::make_shared<rockseis::KdmigElastic2D<float>>(vplmodel, vslmodel, shot2D, simage);

                // Set frequency decimation 
                kdmig->setFreqinc(freqinc);

                // Set logfile
                kdmig->setLogfile("log.txt-" + std::to_string(work.id));

                // Run migration
                kdmig->run();

                // Output image
                simage->write();
                
                // Reset all classes
                shot2D.reset();
                vplmodel.reset();
                vslmodel.reset();
                simage.reset();
                kdmig.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

