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
			PRINT_DOC(# MPI 3d elastic reverse-time migration configuration file);
			PRINT_DOC();
			PRINT_DOC(# Modelling parameters);
			PRINT_DOC(freesurface = "false";  # True if free surface should be on);
			PRINT_DOC(order = "8";  # Order of finite difference stencil);
			PRINT_DOC(lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(source_type = "0"; # Source type 0 - pressure. 1 for Vx. 2 for Vy. 3 for Vz.);
			PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
			PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
			PRINT_DOC(aperty = "1800"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Checkpointing parameters);
			PRINT_DOC(snapmethod = "0";  );
			PRINT_DOC(nsnaps = "11";);
			PRINT_DOC(incore = "true";  );
			PRINT_DOC();
			PRINT_DOC(# Migration parameters);
			PRINT_DOC(nhx = "1";);
			PRINT_DOC(nhy = "1";);
			PRINT_DOC(nhz = "1";);
			PRINT_DOC();
            PRINT_DOC(# Booleans);
            PRINT_DOC(Pimaging = "true";  # Set these to true if imaging of these events is to be made.);
            PRINT_DOC(Simaging = "true";); 
            PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp3d.rss";);
			PRINT_DOC(Vp = "Vs3d.rss";);
			PRINT_DOC(Rho = "Rho3d.rss";);
			PRINT_DOC(Wavelet = "Wav3d.rss";);
			PRINT_DOC(Pimagefile = "Pimage3d.rss";);
			PRINT_DOC(Simagefile = "Simage3d.rss";);
            PRINT_DOC(Vxrecordfile = "Vxshot.rss";);
            PRINT_DOC(Vyrecordfile = "Vyshot.rss";);
            PRINT_DOC(Vzrecordfile = "Vzshot.rss";);
            PRINT_DOC(Snapfile = "snaps.rss";);
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
    float aperty;
    int stype;
    int nhx=1, nhy=1, nhz=1;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Snapfile;

    bool Pimaging;
    std::string Pimagefile;
    std::shared_ptr<rockseis::Image3D<float>> pimage;

    bool Simaging;
    std::string Simagefile;
    std::shared_ptr<rockseis::Image3D<float>> simage;

    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3Di;

    std::string Vyrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3Di;

    std::string Vzrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3Di;

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
    if(Inpar->getPar("source_type", &stype) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Pimaging", &Pimaging) == INPARSE_ERR) status = true;
    if(Pimaging){
        if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Simaging", &Simaging) == INPARSE_ERR) status = true;
    if(Simaging){
        if(Inpar->getPar("Simagefile", &Simagefile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Snapfile", &Snapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vyrecordfile", &Vyrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhy", &nhy) == INPARSE_ERR) status = true;
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
        default:
            rockseis::rs_error("Invalid option of snapshot saving (snapmethod)."); 
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Vxrecordfile);
	
    // Create a global model class
	std::shared_ptr<rockseis::ModelElastic3D<float>> gmodel (new rockseis::ModelElastic3D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));
    // Create a local model class
	std::shared_ptr<rockseis::ModelElastic3D<float>> lmodel (new rockseis::ModelElastic3D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Vxrecordfile); 
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

        // Images
        if(Pimaging){
            pimage = std::make_shared<rockseis::Image3D<float>>(Pimagefile, gmodel, nhx, nhy, nhz);
            pimage->createEmpty();
            for(unsigned long int i=0; i<ngathers; i++) {
                pimage->stackImage(Pimagefile + "-" + std::to_string(i));
            }
        }

        if(Simaging){
            simage = std::make_shared<rockseis::Image3D<float>>(Simagefile, gmodel, nhx, nhy, nhz);
            simage->createEmpty();
            for(unsigned long int i=0; i<ngathers; i++) {
                simage->stackImage(Simagefile + "-" + std::to_string(i));
            }
        }


		
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::RtmElastic3D<float>> rtm;
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
                Sort->setDatafile(Vxrecordfile);
                Vxdata3D = Sort->get3DGather(work.id);
                size_t ntr = Vxdata3D->getNtrace();

                Sort->setDatafile(Vyrecordfile);
                Vydata3D = Sort->get3DGather(work.id);

                Sort->setDatafile(Vzrecordfile);
                Vzdata3D = Sort->get3DGather(work.id);

                lmodel = gmodel->getLocal(Vxdata3D, apertx, aperty, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Vxdata3D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                switch(stype){
                    case 0:
                        source->setField(PRESSURE);
                        break;
                    case 1:
                        source->setField(VX);
                        break;
                    case 2:
                        source->setField(VY);
                        break;
                    case 3:
                        source->setField(VZ);
                        break;
                    default:
                        rs_error("Unknown source type: ", std::to_string(stype));
                        break;
                }

                // Interpolate shot
                Vxdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vxdata3D, Vxdata3Di);
                Vxdata3Di->makeMap(lmodel->getGeom(), GMAP);
                Vxdata3Di->setField(rockseis::VX);

                Vydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vydata3D, Vydata3Di);
                Vydata3Di->makeMap(lmodel->getGeom(), GMAP);
                Vydata3Di->setField(rockseis::VY);

                Vzdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vzdata3D, Vzdata3Di);
                Vzdata3Di->makeMap(lmodel->getGeom(), GMAP);
                Vzdata3Di->setField(rockseis::VZ);

                rtm = std::make_shared<rockseis::RtmElastic3D<float>>(lmodel, source, Vxdata3Di, Vydata3Di, Vzdata3Di, order, snapinc);
   
                // Setting Image objects
                if(Pimaging){ 
                    pimage = std::make_shared<rockseis::Image3D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                    rtm->setPimage(pimage);
                }
                if(Simaging){ 
                    simage = std::make_shared<rockseis::Image3D<float>>(Simagefile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                    rtm->setSimage(simage);
                }

                // Setting Snapshot file 
                rtm->setSnapfile(Snapfile + "-" + std::to_string(work.id));

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
                    default:
                        rockseis::rs_error("Invalid option of snapshot saving."); 
                }

                // Output image
                if(Pimaging){
                    pimage->write();
                }

                if(Simaging){
                    simage->write();
                }

                // Reset all classes
                Vxdata3D.reset();
                Vxdata3Di.reset();
                Vydata3D.reset();
                Vydata3Di.reset();
                Vzdata3D.reset();
                Vzdata3Di.reset();
                lmodel.reset();
                pimage.reset();
                simage.reset();
                rtm.reset();
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

