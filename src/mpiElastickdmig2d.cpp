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
			PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Migration parameters);
			PRINT_DOC(freqinc = "4"; # Integer frequency interval to sum over);
			PRINT_DOC(minfreq = "100.0"; # Minimum frequency to migrate);
			PRINT_DOC(maxfreq = "100.0"; # Maximum frequency to migrate);
			PRINT_DOC(radius = "100.0"; # Radius of traveltime interpolation);
			PRINT_DOC(nhx = "1";);
			PRINT_DOC(nhz = "1";);
			PRINT_DOC();
			PRINT_DOC(# Booleans);
			PRINT_DOC(Gather = "false"; # If surface gathers are to be output);
			PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp2d.rss";);
			PRINT_DOC(Vs = "Vs2d.rss";);
			PRINT_DOC(Sou_ttable = "Sou_ttable2d.rss";);
			PRINT_DOC(Rec_ttable = "Rec_ttable2d.rss";);
			PRINT_DOC(Vxrecordfile = "Pshots2d.rss";);
			PRINT_DOC(Simagefile = "Simage2d.rss";);
			PRINT_DOC(Sgatherfile = "Sgather2d.rss";);
			PRINT_DOC();
		}
        exit(1);
    }
    bool status;
	/* General input parameters */
	int lpml=3;
    float apertx;
    int nhx=1, nhz=1;
	int freqinc;
    float maxfreq;
    float minfreq;
    float radius;
    std::string Vpfile;
    std::string Vsfile;

    std::string sou_ttablefile;
    std::string rec_ttablefile;
    std::string Simagefile;
    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> shot2D;
    std::shared_ptr<rockseis::Image2D<float>> simage;
    bool Gather;
    std::string Sgatherfile;
    std::shared_ptr<rockseis::Data2D<float>> sgather;
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vplmodel;
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vslmodel;

	std::shared_ptr<rockseis::Ttable2D<float>> sou_ttable;
	std::shared_ptr<rockseis::Ttable2D<float>> rec_ttable;

    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("freqinc", &freqinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("minfreq", &minfreq) == INPARSE_ERR) status = true;
    if(Inpar->getPar("maxfreq", &maxfreq) == INPARSE_ERR) status = true;
    if(Inpar->getPar("radius", &radius) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Sou_ttable", &sou_ttablefile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rec_ttable", &rec_ttablefile) == INPARSE_ERR) status = true;
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
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vpgmodel (new rockseis::ModelEikonal2D<float>(Vpfile, lpml));
	std::shared_ptr<rockseis::ModelEikonal2D<float>> vsgmodel (new rockseis::ModelEikonal2D<float>(Vsfile, lpml));

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

        // Image
        simage = std::make_shared<rockseis::Image2D<float>>(Simagefile, vpgmodel, nhx, nhz);
        simage->createEmpty();
        
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
                remove_file(Simagefile + "-" + std::to_string(i));
            }
            sgather->close();
            Fimg->close();
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
                vplmodel->Expand();

                vslmodel = vsgmodel->getLocal(shot2D, apertx, SMAP);
                vslmodel->Expand();

                // Create traveltime table classes
                sou_ttable = std::make_shared<rockseis::Ttable2D<float>>(sou_ttablefile);
                sou_ttable->allocTtable();

                rec_ttable = std::make_shared<rockseis::Ttable2D<float>>(rec_ttablefile);
                rec_ttable->allocTtable();

                // Make image class
                simage = std::make_shared<rockseis::Image2D<float>>(Simagefile + "-" + std::to_string(work.id), vplmodel, nhx, nhz);

                // Create imaging class
                kdmig = std::make_shared<rockseis::KdmigElastic2D<float>>(vplmodel, vslmodel, sou_ttable, rec_ttable, shot2D, simage);

                // Set frequency decimation 
                kdmig->setFreqinc(freqinc);

                // Set minimum and maximum frequency to migrate
                kdmig->setMinfreq(minfreq);
                kdmig->setMaxfreq(maxfreq);

                // Set radius of interpolation
                kdmig->setRadius(radius);

                // Set logfile
                kdmig->setLogfile("log.txt-" + std::to_string(work.id));

                // Run migration
                kdmig->run();

                // Output image
                if(Gather){
                    simage->write();
                }

                // Send result back
                work.status = PARALLEL_IO;
                mpi.sendResult(work);		

                // Stack image
                simage->stackImage_parallel(Simagefile);

                // Reset all classes
                shot2D.reset();
                vplmodel.reset();
                vslmodel.reset();
                simage.reset();
                kdmig.reset();
                sou_ttable.reset();
                rec_ttable.reset();

                // Send result back
                work.status = WORK_FINISHED;
                mpi.sendResult(work);		
            }
        }
    }
}

