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
            PRINT_DOC(source_type = "0"; # Source type: 0 - pressure; 1 for Ux; 2 for Uy; 3 for Uz.);
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
            PRINT_DOC(Uxrecordfile = "Uxshot.rss"; # Input observed data);
			PRINT_DOC(Uxmodelledfile = "Uxmod3d.rss"; # File to output modelled data);
			PRINT_DOC(Uxresidualfile = "Uxres3d.rss"; # File to output residuals);
            PRINT_DOC(Uyrecordfile = "Uyshot.rss"; # Input observed data);
			PRINT_DOC(Uymodelledfile = "Uymod3d.rss"; # File to output modelled data);
			PRINT_DOC(Uyresidualfile = "Uyres3d.rss"; # File to output residuals);
            PRINT_DOC(Uzrecordfile = "Uzshot.rss"; # Input observed data);
			PRINT_DOC(Uzmodelledfile = "Uzmod3d.rss"; # File to output modelled data);
			PRINT_DOC(Uzresidualfile = "Uzres3d.rss"; # File to output residuals);
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

    std::string Uxrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Uxdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Uxdata3Di;

    std::string Uxmodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Uxdatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Uxdatamod3Di;

    std::string Uxresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Uxdatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Uxdatares3Di;

    std::string Uyrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Uydata3D;
    std::shared_ptr<rockseis::Data3D<float>> Uydata3Di;

    std::string Uymodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Uydatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Uydatamod3Di;

    std::string Uyresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Uydatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Uydatares3Di;

    std::string Uzrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Uzdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Uzdata3Di;

    std::string Uzmodelledfile;
    std::shared_ptr<rockseis::Data3D<float>> Uzdatamod3D;
    std::shared_ptr<rockseis::Data3D<float>> Uzdatamod3Di;

    std::string Uzresidualfile;
    std::shared_ptr<rockseis::Data3D<float>> Uzdatares3D;
    std::shared_ptr<rockseis::Data3D<float>> Uzdatares3Di;

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
    if(Inpar->getPar("Uxrecordfile", &Uxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uxresidualfile", &Uxresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uxmodelledfile", &Uxmodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uyrecordfile", &Uyrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uyresidualfile", &Uyresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uymodelledfile", &Uymodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uzrecordfile", &Uzrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uzresidualfile", &Uzresidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uzmodelledfile", &Uzmodelledfile) == INPARSE_ERR) status = true;
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
    Sort->setDatafile(Uxrecordfile);
	
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
        Sort->createShotmap(Uxrecordfile); 
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
        std::shared_ptr<rockseis::Data3D<float>> Uxdata3D (new rockseis::Data3D<float>(Uxrecordfile));

        // Create modelling and residual data files
        Uxdatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uxdatamod3D->setFile(Uxmodelledfile);
        Uxdatamod3D->createEmpty(Uxdata3D->getNtrace());
        Uxdatares3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uxdatares3D->setFile(Uxresidualfile);
        Uxdatares3D->createEmpty(Uxdata3D->getNtrace());

        Uydatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uydatamod3D->setFile(Uymodelledfile);
        Uydatamod3D->createEmpty(Uxdata3D->getNtrace());
        Uydatares3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uydatares3D->setFile(Uyresidualfile);
        Uydatares3D->createEmpty(Uxdata3D->getNtrace());

        Uzdatamod3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uzdatamod3D->setFile(Uzmodelledfile);
        Uzdatamod3D->createEmpty(Uxdata3D->getNtrace());
        Uzdatares3D = std::make_shared<rockseis::Data3D<float>>(1, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
        Uzdatares3D->setFile(Uzresidualfile);
        Uzdatares3D->createEmpty(Uxdata3D->getNtrace());
       
		// Create work queue
		for(long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
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

        for(long int i=0; i<ngathers; i++) {
            if(Vpgrad){
                vpgrad->stackImage(Vpgradfile + "-" + std::to_string(i));
                remove_file(Vpgradfile + "-" + std::to_string(i));
            }
            if(Vsgrad){
                vsgrad->stackImage(Vsgradfile + "-" + std::to_string(i));
                remove_file(Vsgradfile + "-" + std::to_string(i));
            }
            if(Rhograd){
                rhograd->stackImage(Rhogradfile + "-" + std::to_string(i));
                remove_file(Rhogradfile + "-" + std::to_string(i));
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
                Sort->setDatafile(Uxrecordfile);
                Uxdata3D = Sort->get3DGather(work.id);
                size_t ntr = Uxdata3D->getNtrace();

                Sort->setDatafile(Uyrecordfile);
                Uydata3D = Sort->get3DGather(work.id);

                Sort->setDatafile(Uzrecordfile);
                Uzdata3D = Sort->get3DGather(work.id);

                lmodel = gmodel->getLocal(Uxdata3D, apertx, aperty, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Uxdata3D);
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
                Uxdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uxdata3D, Uxdata3Di);
                Uxdata3Di->makeMap(lmodel->getGeom(), GMAP);
                Uxdata3Di->setField(rockseis::VX);

                Uydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uydata3D, Uydata3Di);
                Uydata3Di->makeMap(lmodel->getGeom(), GMAP);
                Uydata3Di->setField(rockseis::VY);

                Uzdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Uzdata3D, Uzdata3Di);
                Uzdata3Di->makeMap(lmodel->getGeom(), GMAP);
                Uzdata3Di->setField(rockseis::VZ);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic3D<float>>(lmodel, source, Uxdata3Di, Uydata3Di, Uzdata3Di, order, snapinc);

                // Create modelled and residual data objects 
                Uxdatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatamod3D->copyCoords(Uxdata3D);
                Uxdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uxdatamod3D->setField(rockseis::VX);
                fwi->setDatamodUx(Uxdatamod3D);
                Uxdatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uxdatares3D->copyCoords(Uxdata3D);
                Uxdatares3D->makeMap(lmodel->getGeom(), GMAP);
                Uxdatares3D->setField(rockseis::VX);
                fwi->setDataresUx(Uxdatares3D);

                Uydatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uydatamod3D->copyCoords(Uydata3D);
                Uydatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uydatamod3D->setField(rockseis::VY);
                fwi->setDatamodUy(Uydatamod3D);
                Uydatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uydatares3D->copyCoords(Uydata3D);
                Uydatares3D->makeMap(lmodel->getGeom(), GMAP);
                Uydatares3D->setField(rockseis::VY);
                fwi->setDataresUy(Uydatares3D);

                Uzdatamod3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatamod3D->copyCoords(Uzdata3D);
                Uzdatamod3D->makeMap(lmodel->getGeom(), GMAP);
                Uzdatamod3D->setField(rockseis::VZ);
                fwi->setDatamodUz(Uzdatamod3D);
                Uzdatares3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Uzdatares3D->copyCoords(Uzdata3D);
                Uzdatares3D->makeMap(lmodel->getGeom(), GMAP);
                Uzdatares3D->setField(rockseis::VZ);
                fwi->setDataresUz(Uzdatares3D);

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
                Uxdatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uxdatamod3Di->setFile(Uxmodelledfile);
                interp->interp(Uxdatamod3D, Uxdatamod3Di);
                Sort->put3DGather(Uxdatamod3Di, work.id);

                Uxdatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uxdatares3Di->setFile(Uxresidualfile);
                interp->interp(Uxdatares3D, Uxdatares3Di);
                Sort->put3DGather(Uxdatares3Di, work.id);

                Uydatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uydatamod3Di->setFile(Uymodelledfile);
                interp->interp(Uydatamod3D, Uydatamod3Di);
                Sort->put3DGather(Uydatamod3Di, work.id);

                Uydatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uydatares3Di->setFile(Uyresidualfile);
                interp->interp(Uydatares3D, Uydatares3Di);
                Sort->put3DGather(Uydatares3Di, work.id);

                Uzdatamod3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uzdatamod3Di->setFile(Uzmodelledfile);
                interp->interp(Uzdatamod3D, Uzdatamod3Di);
                Sort->put3DGather(Uzdatamod3Di, work.id);

                Uzdatares3Di = std::make_shared<rockseis::Data3D<float>>(ntr, Uxdata3D->getNt(), Uxdata3D->getDt(), Uxdata3D->getOt());
                Uzdatares3Di->setFile(Uzresidualfile);
                interp->interp(Uzdatares3D, Uzdatares3Di);
                Sort->put3DGather(Uzdatares3Di, work.id);

                // Reset all classes
                Uxdata3D.reset();
                Uxdata3Di.reset();
                Uxdatamod3D.reset();
                Uxdatamod3Di.reset();
                Uxdatares3D.reset();
                Uxdatares3Di.reset();
                Uydata3D.reset();
                Uydata3Di.reset();
                Uydatamod3D.reset();
                Uydatamod3Di.reset();
                Uydatares3D.reset();
                Uydatares3Di.reset();
                Uzdata3D.reset();
                Uzdata3Di.reset();
                Uzdatamod3D.reset();
                Uzdatamod3Di.reset();
                Uzdatares3D.reset();
                Uzdatares3Di.reset();
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

