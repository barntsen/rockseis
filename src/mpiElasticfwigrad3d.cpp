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
			PRINT_DOC(# MPI 3d elastic full-waveform inversion gradient configuration file);
			PRINT_DOC();
			PRINT_DOC(# Modelling parameters);
			PRINT_DOC(freesurface = "false";  # True if free surface should be on);
			PRINT_DOC(order = "8";  # Order of finite difference stencil);
			PRINT_DOC(lpml = "10"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(source_type = "0"; # Source type: 0 - pressure; 1 for Vx; 2 for Vy; 3 for Vz.);
			PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
			PRINT_DOC(apertx = "900"; # Aperture for local model (source is in the middle));
			PRINT_DOC(aperty = "900"; # Aperture for local model (source is in the middle));
			PRINT_DOC();
			PRINT_DOC(# Checkpointing parameters);
			PRINT_DOC(snapmethod = "0"; # 0 - fullcheckpointing; 1 - optimal checkpointing );
			PRINT_DOC(nsnaps = "11"; # Number of checkpoints to store);
			PRINT_DOC(incore = "true"; # Do checkpointing in memory);
			PRINT_DOC();
            PRINT_DOC(# Booleans);
            PRINT_DOC(Vpgrad = "true";  # Set these to true if imaging of these events is to be made.);
            PRINT_DOC(Vsgrad = "true"; );
            PRINT_DOC(Rhograd = "true"; ); 
            PRINT_DOC(Wavgrad = "true"; );
			PRINT_DOC();
			PRINT_DOC(# Fwi parameters);
			PRINT_DOC(misfit_type= "0"; # 0 - Difference; 1 - Correlation );
            PRINT_DOC();
			PRINT_DOC(# Files);
			PRINT_DOC(Vp = "Vp3d.rss"; # Input Vp model);
			PRINT_DOC(Vs = "Vs3d.rss"; # Input Vs model);
			PRINT_DOC(Rho = "Rho3d.rss"; # Input Rho model);
			PRINT_DOC(Wavelet = "Wav3d.rss"; # Input Wav model);
			PRINT_DOC(Vpgradfile = "Vpgrad3d.rss"; # File to output gradient with respect to Vp);
			PRINT_DOC(Vsgradfile = "Vsgrad3d.rss"; # File to output gradient with respect to Vs);
			PRINT_DOC(Rhogradfile = "Rhograd3d.rss"; # File to output gradient with respect to Rho);
			PRINT_DOC(Wavgradfile = "Wavgrad3d.rss"; # File to output gradient with respect to Wav);
            PRINT_DOC(Vxrecordfile = "Vxshot.rss"; # Input observed data);
			PRINT_DOC(Vxmodelledfile = "Vxmod3d.rss"; # File to output modelled data);
			PRINT_DOC(Vxresidualfile = "Vxres3d.rss"; # File to output residuals);
            PRINT_DOC(Vyrecordfile = "Vyshot.rss"; # Input observed data);
			PRINT_DOC(Vymodelledfile = "Vymod3d.rss"; # File to output modelled data);
			PRINT_DOC(Vyresidualfile = "Vyres3d.rss"; # File to output residuals);
            PRINT_DOC(Vzrecordfile = "Vzshot.rss"; # Input observed data);
			PRINT_DOC(Vzmodelledfile = "Vzmod3d.rss"; # File to output modelled data);
			PRINT_DOC(Vzresidualfile = "Vzres3d.rss"; # File to output residuals);
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
    bool Vpgrad, Vsgrad, Rhograd, Wavgrad;
	int order;
	int snapinc;
	int nsnaps = 0;
	int snapmethod;
	int misfit_type;
    float apertx;
    float aperty;
    int stype;
    int nhx=1, nhy=1, nhz=1;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Snapfile;

    std::string Vpgradfile;
    std::shared_ptr<rockseis::Image3D<float>> vpgrad;

    std::string Vsgradfile;
    std::shared_ptr<rockseis::Image3D<float>> vsgrad;

    std::string Rhogradfile;
    std::shared_ptr<rockseis::Image3D<float>> rhograd;

    std::string Wavgradfile;
    std::shared_ptr<rockseis::Data3D<float>> wavgrad;

    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata3Di;

    std::string Vxmodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Vxdatamod3Di;

    std::string Vxresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Vxdatares3Di;

    std::string Vyrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vydata3Di;

    std::string Vymodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Vydatamod3Di;

    std::string Vyresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Vydatares3Di;

    std::string Vzrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata3Di;

    std::string Vzmodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Vzdatamod3Di;

    std::string Vzresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Vzdatares3Di;

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
    if(Inpar->getPar("Vpgrad", &Vpgrad) == INPARSE_ERR) status = true;
    if(Vpgrad){
        if(Inpar->getPar("Vpgradfile", &Vpgradfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vsgrad", &Vsgrad) == INPARSE_ERR) status = true;
    if(Vsgrad){
        if(Inpar->getPar("Vsgradfile", &Vsgradfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Rhograd", &Rhograd) == INPARSE_ERR) status = true;
    if(Rhograd){
        if(Inpar->getPar("Rhogradfile", &Rhogradfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Wavgrad", &Wavgrad) == INPARSE_ERR) status = true;
    if(Wavgrad){
        if(Inpar->getPar("Wavgradfile", &Wavgradfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Snapfile", &Snapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vxresidualfile", &Vxresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vxmodelledfile", &Vxmodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vyrecordfile", &Vyrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vyresidualfile", &Vyresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vymodelledfile", &Vymodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vzresidualfile", &Vzresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vzmodelledfile", &Vzmodelledfile) == INPARSE_ERR) status = true;
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

        // Create empty file for wavelet gradients
        if(Wavgrad){
            wavgrad = std::make_shared<rockseis::Data3D<float>>(1, source->getNt(), source->getDt(), 0.0);
            wavgrad->setFile(Wavgradfile);
            wavgrad->createEmpty(ngathers);
        }

        // Create a data class for the recorded data in order to get parameters from file
        std::shared_ptr<rockseis::Data3D<float>> Vxdata3D (new rockseis::Data3D<float>(Vxrecordfile));

        // Create modelling and residual data files
        Vxdatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vxdatamod3D->setFile(Vxmodelledfile);
        Vxdatamod3D->createEmpty(Vxdata3D->getNtrace());
        Vxdatares3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vxdatares3D->setFile(Vxresidualfile);
        Vxdatares3D->createEmpty(Vxdata3D->getNtrace());

        Vydatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vydatamod3D->setFile(Vymodelledfile);
        Vydatamod3D->createEmpty(Vxdata3D->getNtrace());
        Vydatares3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vydatares3D->setFile(Vyresidualfile);
        Vydatares3D->createEmpty(Vxdata3D->getNtrace());

        Vzdatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vzdatamod3D->setFile(Vzmodelledfile);
        Vzdatamod3D->createEmpty(Vxdata3D->getNtrace());
        Vzdatares3D = std::make_shared<rockseis::Data3D<float>>(1, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
        Vzdatares3D->setFile(Vzresidualfile);
        Vzdatares3D->createEmpty(Vxdata3D->getNtrace());
       
		// Create work queue
		for(unsigned long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0,{'\0'}});
			mpi.addWork(work);
		}

		// Perform work in parallel
		mpi.performWork();

        // Images
        if(Vpgrad){
            vpgrad = std::make_shared<rockseis::Image3D<float>>(Vpgradfile, gmodel, nhx, nhy, nhz);
            vpgrad->createEmpty();
        }
        if(Vsgrad){

            vsgrad = std::make_shared<rockseis::Image3D<float>>(Vsgradfile, gmodel, nhx, nhy, nhz);
            vsgrad->createEmpty();
        }
        if(Rhograd){

            rhograd = std::make_shared<rockseis::Image3D<float>>(Rhogradfile, gmodel, nhx, nhy, nhz);
            rhograd->createEmpty();
        }

        for(unsigned long int i=0; i<ngathers; i++) {
            if(Vpgrad){
                vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
            }
            if(Vsgrad){
                vsgrad->stackImage(Vsgradfile + "-" + std::to_string(i));
            }
            if(Rhograd){
                rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
            }
        }
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::FwiElastic3D<float>> fwi;
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

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic3D<float>>(lmodel, source, Vxdata3Di, Vydata3Di, Vzdata3Di, order, snapinc);

                // Create modelled and residual data objects 
                Vxdatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatamod3D->copyCoords(Vxdata3D);
                Vxdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatamod3D->setField(rockseis::VX);
                fwi->setDatamodVx(Vxdatamod3D);
                Vxdatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatares3D->copyCoords(Vxdata3D);
                Vxdatares3D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatares3D->setField(rockseis::VX);
                fwi->setDataresVx(Vxdatares3D);

                Vydatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vydatamod3D->copyCoords(Vydata3D);
                Vydatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Vydatamod3D->setField(rockseis::VY);
                fwi->setDatamodVy(Vydatamod3D);
                Vydatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vydatares3D->copyCoords(Vydata3D);
                Vydatares3D->makeMap(lmodel->getGeom(), GMAP);
                Vydatares3D->setField(rockseis::VY);
                fwi->setDataresVy(Vydatares3D);

                Vzdatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatamod3D->copyCoords(Vzdata3D);
                Vzdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatamod3D->setField(rockseis::VZ);
                fwi->setDatamodVz(Vzdatamod3D);
                Vzdatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatares3D->copyCoords(Vzdata3D);
                Vzdatares3D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatares3D->setField(rockseis::VZ);
                fwi->setDataresVz(Vzdatares3D);

                // Setting misfit type
                fwi->setMisfit_type(fwimisfit);

                // Creating gradient objects and setting them up in fwi class
                if(Vpgrad){
                    vpgrad = std::make_shared<rockseis::Image3D<float>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                fwi->setVpgrad(vpgrad);
                }
                if(Vsgrad){
                    vsgrad = std::make_shared<rockseis::Image3D<float>>(Vsgradfile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                fwi->setVsgrad(vsgrad);
                }
                if(Rhograd){
                    rhograd = std::make_shared<rockseis::Image3D<float>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
                fwi->setRhograd(rhograd);
                }
   

                if(Wavgrad){
                    wavgrad = std::make_shared<rockseis::Data3D<float>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
                    // Copy geometry
                    wavgrad->copyCoords(source);
                    wavgrad->makeMap(lmodel->getGeom(), SMAP);
                    switch(stype){
                        case 0:
                            wavgrad->setField(PRESSURE);
                            break;
                        case 1:
                            wavgrad->setField(VX);
                            break;
                        case 2:
                            wavgrad->setField(VY);
                            break;
                        case 3:
                            wavgrad->setField(VZ);
                            break;
                        default:
                            rs_error("Unknown wavgrad type: ", std::to_string(stype));
                            break;
                    }
                    fwi->setWavgrad(wavgrad);
                }

                // Setting Snapshot file 
                fwi->setSnapfile(Snapfile + "-" + std::to_string(work.id));

                // Setting Snapshot parameters
                fwi->setNcheck(nsnaps);
                fwi->setIncore(incore);

                // Set logfile
                fwi->setLogfile("log.txt-" + std::to_string(work.id));

                // Stagger model
                lmodel->staggerModels();

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
                if(Vpgrad){
                    vpgrad->write();
                }
                if(Vsgrad){
                    vsgrad->write();
                }
                if(Rhograd){
                    rhograd->write();
                }
                if(Wavgrad){
                    wavgrad->putTrace(Wavgradfile, work.id);
                }

                // Output modelled and residual data
                Vxdatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vxdatamod3Di->setFile(Vxmodelledfile);
                interp->interp(Vxdatamod3D, Vxdatamod3Di);
                Sort->put3DGather(Vxdatamod3Di, work.id);

                Vxdatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vxdatares3Di->setFile(Vxresidualfile);
                interp->interp(Vxdatares3D, Vxdatares3Di);
                Sort->put3DGather(Vxdatares3Di, work.id);

                Vydatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vydatamod3Di->setFile(Vymodelledfile);
                interp->interp(Vydatamod3D, Vydatamod3Di);
                Sort->put3DGather(Vydatamod3Di, work.id);

                Vydatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vydatares3Di->setFile(Vyresidualfile);
                interp->interp(Vydatares3D, Vydatares3Di);
                Sort->put3DGather(Vydatares3Di, work.id);

                Vzdatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vzdatamod3Di->setFile(Vzmodelledfile);
                interp->interp(Vzdatamod3D, Vzdatamod3Di);
                Sort->put3DGather(Vzdatamod3Di, work.id);

                Vzdatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Vxdata3D->getNt(), Vxdata3D->getDt(), Vxdata3D->getOt());
                Vzdatares3Di->setFile(Vzresidualfile);
                interp->interp(Vzdatares3D, Vzdatares3Di);
                Sort->put3DGather(Vzdatares3Di, work.id);

                // Reset all classes
                Vxdata3D.reset();
                Vxdata3Di.reset();
                Vxdatamod3D.reset();
                Vxdatamod3Di.reset();
                Vxdatares3D.reset();
                Vxdatares3Di.reset();
                Vydata3D.reset();
                Vydata3Di.reset();
                Vydatamod3D.reset();
                Vydatamod3Di.reset();
                Vydatares3D.reset();
                Vydatares3Di.reset();
                Vzdata3D.reset();
                Vzdata3Di.reset();
                Vzdatamod3D.reset();
                Vzdatamod3Di.reset();
                Vzdatares3D.reset();
                Vzdatares3Di.reset();
                lmodel.reset();
                vpgrad.reset();
                vsgrad.reset();
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

