#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "modelling.h"
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
            PRINT_DOC(# MPI 2d acoustic modelling default configuration file);
            PRINT_DOC();
            PRINT_DOC(# Modelling parameters);
            PRINT_DOC(        freesurface = "true";  # True if free surface should be on);
            PRINT_DOC(            order = "8";  # Order of finite difference stencil 2-8);
            PRINT_DOC(            lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(            snapinc = "10"; # Snap interval in multiples of modelling interval);
            PRINT_DOC(            dtrec = "4e-3"; # Recording interval in seconds);
            PRINT_DOC(            apertx = "900"; # Aperture for local model (source is in the middle));
            PRINT_DOC(            aperty = "900"; # Aperture for local model (source is in the middle));
            PRINT_DOC();
            PRINT_DOC(# Booleans);
            PRINT_DOC(            Precord = "true";  # Set these to true if recording or snapshoting is to be made.);
            PRINT_DOC(            Vxrecord = "false";);
            PRINT_DOC(        Vyrecord = "false";);
            PRINT_DOC(        Vzrecord = "false";);
            PRINT_DOC(        Psnap = "false";);
            PRINT_DOC(        Vxsnap = "false";);
            PRINT_DOC(        Vysnap = "false";);
            PRINT_DOC(        Vzsnap = "false";);
            PRINT_DOC();
            PRINT_DOC(# Files);
            PRINT_DOC(        Vp = "Vp2d.rss";);
            PRINT_DOC(        Rho = "Rho2d.rss";);
            PRINT_DOC(        Wavelet = "Wav2d.rss";);
            PRINT_DOC(        Input_Pdata = "3Ddata.rss";);
            PRINT_DOC(        Datum_survey = "3Dsurvey.rss";);
            PRINT_DOC(        Precordfile = "Pshot.rss";);
            PRINT_DOC(        Vxrecordfile = "Vxshot.rss";);
            PRINT_DOC(        Vyrecordfile = "Vyshot.rss";);
            PRINT_DOC(        Vzrecordfile = "Vzshot.rss";);
            PRINT_DOC(        Psnapfile = "Psnaps.rss";);
            PRINT_DOC(        Vxsnapfile = "Vxsnaps.rss";);
            PRINT_DOC(        Vysnapfile = "Vysnaps.rss";);
            PRINT_DOC(        Vzsnapfile = "Vzsnaps.rss";);
        }
        exit(1);
    }
    bool status;
    /* General input parameters */
    int lpml;
    bool fs;
    int order;
    int snapinc;
    float apertx;
    float aperty;
    float dtrec;
    std::string Indatafile;
    std::string Surveyfile;
    std::string Waveletfile;

    std::shared_ptr<rockseis::Data3D<float>> Indata3D;
    std::shared_ptr<rockseis::Data3D<float>> Indata3Di;
    std::string Vpfile;
    std::string Rhofile;
    bool Psnap=0, Precord=0;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data3D<float>> Pdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Pdata3Di;

    bool Vxsnap=0, Vxrecord=0;
    std::string Vxsnapfile;
    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3Di;

    bool Vysnap=0, Vyrecord=0;
    std::string Vysnapfile;
    std::string Vyrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3Di;

    bool Vzsnap=0, Vzrecord=0;
    std::string Vzsnapfile;
    std::string Vzrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3Di;

    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> lmodel;

    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("dtrec", &dtrec) == INPARSE_ERR) status = true;
    if(Inpar->getPar("order", &order) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Input_Pdata", &Indatafile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Datum_survey", &Surveyfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Psnap", &Psnap) == INPARSE_ERR) status = true;
    if(Psnap){
        if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vxsnap", &Vxsnap) == INPARSE_ERR) status = true;
    if(Vxsnap){
        if(Inpar->getPar("Vxsnapfile", &Vxsnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vysnap", &Vysnap) == INPARSE_ERR) status = true;
    if(Vysnap){
        if(Inpar->getPar("Vysnapfile", &Vysnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vzsnap", &Vzsnap) == INPARSE_ERR) status = true;
    if(Vzsnap){
        if(Inpar->getPar("Vzsnapfile", &Vzsnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Precord", &Precord) == INPARSE_ERR) status = true;
    if(Precord){
        if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vxrecord", &Vxrecord) == INPARSE_ERR) status = true;
    if(Vxrecord){
        if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vzrecord", &Vzrecord) == INPARSE_ERR) status = true;
    if(Vzrecord){
        if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort classes
    std::shared_ptr<rockseis::Sort<float>> Sort_in (new rockseis::Sort<float>());
    Sort_in->setDatafile(Indatafile);
    std::shared_ptr<rockseis::Sort<float>> Sort_out (new rockseis::Sort<float>());
    Sort_out->setDatafile(Surveyfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> gmodel (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

    // Compute record length in samples
    size_t ntrec; 
    ntrec = (size_t) rintf((source->getNt()-1)*source->getDt()/dtrec + 1);

    if(mpi.getRank() == 0) {
	// Master
	Sort_in->createShotmap(Indatafile); 
	Sort_in->setKmapfile("kmap_in.rss");
	Sort_in->setSmapfile("smap_in.rss");
	Sort_in->writeKeymap();
	Sort_in->writeSortmap();

	Sort_out->createShotmap(Surveyfile); 
	Sort_out->setKmapfile("kmap_out.rss");
	Sort_out->setSmapfile("smap_out.rss");
	Sort_out->writeKeymap();
	Sort_out->writeSortmap();

        // Get number of shots
        size_t ngathers =  Sort_out->getNensemb();
        
        if(Precord){
            // Create an empty data file
            Sort_out->createEmptydataset(Precordfile, ntrec, dtrec, 0.0);
        }

        if(Vxrecord){
            // Create an empty data file
            Sort_out->createEmptydataset(Vxrecordfile, ntrec, dtrec, 0.0);
        }
        if(Vyrecord){
            // Create an empty data file
            Sort_out->createEmptydataset(Vyrecordfile, ntrec, dtrec, 0.0);
        }
        if(Vzrecord){
            // Create an empty data file
            Sort_out->createEmptydataset(Vzrecordfile, ntrec, dtrec, 0.0);
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
        std::shared_ptr<rockseis::ModellingAcoustic3D<float>> modelling;
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
	    Sort_in->setKmapfile("kmap_in.rss");
	    Sort_in->setSmapfile("smap_in.rss");
	    Sort_in->readKeymap();
	    Sort_in->readSortmap();

	    Sort_out->setKmapfile("kmap_out.rss");
	    Sort_out->setSmapfile("smap_out.rss");
	    Sort_out->readKeymap();
	    Sort_out->readSortmap();

                // Get the input data
                Sort_in->setDatafile(Indatafile);
                Indata3D = Sort_in->get3DGather(work.id);
                size_t ntr_in = Indata3D->getNtrace();
                lmodel = gmodel->getLocal(Indata3D, apertx, aperty, SMAP);

                // Interpolate input data to avoid instability       
                Indata3Di = std::make_shared<rockseis::Data3D<float>>(ntr_in, source->getNt(), source->getDt(), 0.0);
                interp->interp(Indata3D, Indata3Di);
                Indata3Di->makeMap(lmodel->getGeom(), GMAP);
                Indata3Di->copyGmap2Smap();

                // Get the Datum geometry
                Shotgeom = Sort_out->get3DGather(work.id);
                size_t ntr_out = Shotgeom->getNtrace();

                modelling = std::make_shared<rockseis::ModellingAcoustic3D<float>>(lmodel, Indata3Di, order, snapinc);

                // Set logfile
                modelling->setLogfile("log.txt-" + std::to_string(work.id));

                // Setting Snapshot file 
                if(Psnap){
                    modelling->setSnapP(Psnapfile + "-" + std::to_string(work.id));
                }
                if(Vxsnap){
                    modelling->setSnapVx(Vxsnapfile + "-" + std::to_string(work.id));
                }
                if(Vysnap){
                    modelling->setSnapVy(Vysnapfile + "-" + std::to_string(work.id));
                }
                if(Vzsnap){
                    modelling->setSnapVz(Vzsnapfile + "-" + std::to_string(work.id));
                }

                // Setting Record
                if(Precord){
                    Pdata3D = std::make_shared<rockseis::Data3D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
                    Pdata3D->setField(rockseis::PRESSURE);
                    // Copy geometry to Data
                    Pdata3D->copyCoords(Shotgeom);
                    Pdata3D->makeMap(lmodel->getGeom());
                    modelling->setRecP(Pdata3D);
                }
                if(Vxrecord){
                    Vxdata3D = std::make_shared<rockseis::Data3D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
                    Vxdata3D->setField(rockseis::VX);
                    // Copy geometry to Data
                    Vxdata3D->copyCoords(Shotgeom);
                    Vxdata3D->makeMap(lmodel->getGeom());
                    modelling->setRecVx(Vxdata3D);
                }
                if(Vyrecord){
                    Vydata3D = std::make_shared<rockseis::Data3D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
                    Vydata3D->setField(rockseis::VY);
                    // Copy geometry to Data
                    Vydata3D->copyCoords(Shotgeom);
                    Vydata3D->makeMap(lmodel->getGeom());
                    modelling->setRecVy(Vydata3D);
                }
                if(Vzrecord){
                    Vzdata3D = std::make_shared<rockseis::Data3D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
                    Vzdata3D->setField(rockseis::VZ);
                    // Copy geometry to Data
                    Vzdata3D->copyCoords(Shotgeom);
                    Vzdata3D->makeMap(lmodel->getGeom());
                    modelling->setRecVz(Vzdata3D);
                }

                // Stagger model
                lmodel->staggerModels();

                // Run modelling 
                modelling->run();

                // Output record
                if(Precord){
                    Pdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr_out, ntrec, dtrec, 0.0);
                    Pdata3Di->setFile(Precordfile);
                    interp->interp(Pdata3D, Pdata3Di);
                    Sort_out->put3DGather(Pdata3Di, work.id);
                }
                if(Vxrecord){
                    Vxdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr_out, ntrec, dtrec, 0.0);
                    Vxdata3Di->setFile(Vxrecordfile);
                    interp->interp(Vxdata3D, Vxdata3Di);
                    Sort_out->put3DGather(Vxdata3Di, work.id);
                }
                if(Vyrecord){
                    Vydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr_out, ntrec, dtrec, 0.0);
                    Vydata3Di->setFile(Vyrecordfile);
                    interp->interp(Vydata3D, Vydata3Di);
                    Sort_out->put3DGather(Vydata3Di, work.id);
                }
                if(Vzrecord){
                    Vzdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr_out, ntrec, dtrec, 0.0);
                    Vzdata3Di->setFile(Vzrecordfile);
                    interp->interp(Vzdata3D, Vzdata3Di);
                    Sort_out->put3DGather(Vzdata3Di, work.id);
                }

                // Reset all classes
                Pdata3D.reset();
                Pdata3Di.reset();
                Shotgeom.reset();
                lmodel.reset();
                modelling.reset();
                if(Precord){
                    Pdata3D.reset();
                    Pdata3Di.reset();
                }
                if(Vxrecord){
                    Vxdata3D.reset();
                    Vxdata3Di.reset();
                }
                if(Vyrecord){
                    Vydata3D.reset();
                    Vydata3Di.reset();
                }
                if(Vzrecord){
                    Vzdata3D.reset();
                    Vzdata3Di.reset();
                }
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

