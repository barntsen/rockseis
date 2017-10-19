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
			PRINT_DOC(# MPI 2d elastic full-waveform inversion gradient configuration file);
			PRINT_DOC();
			PRINT_DOC(# Modelling parameters);
			PRINT_DOC(freesurface = "false";  # True if free surface should be on);
			PRINT_DOC(order = "8";  # Order of finite difference stencil);
			PRINT_DOC(lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(source_type = "0"; # Source type 0 - pressure. 1 for Vx. 3 for Vz.);
			PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
			PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
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
			PRINT_DOC(Vp = "Vp2d.rss"; # Input Vp model);
			PRINT_DOC(Vs = "Vs2d.rss"; # Input Vs model);
			PRINT_DOC(Rho = "Rho2d.rss"; # Input Rho model);
			PRINT_DOC(Wavelet = "Wav2d.rss"; # Input Wav model);
			PRINT_DOC(Vpgradfile = "Vpgrad2d.rss"; # File to output gradient with respect to Vp);
			PRINT_DOC(Vsgradfile = "Vsgrad2d.rss"; # File to output gradient with respect to Vs);
			PRINT_DOC(Rhogradfile = "Rhograd2d.rss"; # File to output gradient with respect to Rho);
			PRINT_DOC(Wavgradfile = "Wavgrad2d.rss"; # File to output gradient with respect to Wav);
            PRINT_DOC(Vxrecordfile = "Vxshot.rss"; # Input observed data);
			PRINT_DOC(Vxmodelledfile = "Vxmod2d.rss"; # File to output modelled data);
			PRINT_DOC(Vxresidualfile = "Vxres2d.rss"; # File to output residuals);
            PRINT_DOC(Vzrecordfile = "Vzshot.rss"; # Input observed data);
			PRINT_DOC(Vzmodelledfile = "Vzmod2d.rss"; # File to output modelled data);
			PRINT_DOC(Vzresidualfile = "Vzres2d.rss"; # File to output residuals);
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
    int stype;
    int nhx=1, nhz=1;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Snapfile;

    std::string Vpgradfile;
    std::shared_ptr<rockseis::Image2D<float>> vpgrad;

    std::string Vsgradfile;
    std::shared_ptr<rockseis::Image2D<float>> vsgrad;

    std::string Rhogradfile;
    std::shared_ptr<rockseis::Image2D<float>> rhograd;

    std::string Wavgradfile;
    std::shared_ptr<rockseis::Data2D<float>> wavgrad;

    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> Vxdata2D;
    std::shared_ptr<rockseis::Data2D<float>> Vxdata2Di;

    std::string Vxmodelledfile;
    std::shared_ptr<rockseis::Data2D<float>> Vxdatamod2D;
    std::shared_ptr<rockseis::Data2D<float>> Vxdatamod2Di;

    std::string Vxresidualfile;
    std::shared_ptr<rockseis::Data2D<float>> Vxdatares2D;
    std::shared_ptr<rockseis::Data2D<float>> Vxdatares2Di;

    std::string Vzrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> Vzdata2D;
    std::shared_ptr<rockseis::Data2D<float>> Vzdata2Di;

    std::string Vzmodelledfile;
    std::shared_ptr<rockseis::Data2D<float>> Vzdatamod2D;
    std::shared_ptr<rockseis::Data2D<float>> Vzdatamod2Di;

    std::string Vzresidualfile;
    std::shared_ptr<rockseis::Data2D<float>> Vzdatares2D;
    std::shared_ptr<rockseis::Data2D<float>> Vzdatares2Di;

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
	std::shared_ptr<rockseis::ModelElastic2D<float>> gmodel (new rockseis::ModelElastic2D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));
    // Create a local model class
	std::shared_ptr<rockseis::ModelElastic2D<float>> lmodel (new rockseis::ModelElastic2D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(Waveletfile));

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
            wavgrad = std::make_shared<rockseis::Data2D<float>>(1, source->getNt(), source->getDt(), 0.0);
            wavgrad->setFile(Wavgradfile);
            wavgrad->createEmpty(ngathers);
        }

        // Create a data class for the recorded data in order to get parameters from file
        std::shared_ptr<rockseis::Data2D<float>> Vxdata2D (new rockseis::Data2D<float>(Vxrecordfile));

        // Create modelling and residual data files
        Vxdatamod2D = std::make_shared<rockseis::Data2D<float>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
        Vxdatamod2D->setFile(Vxmodelledfile);
        Vxdatamod2D->createEmpty(Vxdata2D->getNtrace());
        Vxdatares2D = std::make_shared<rockseis::Data2D<float>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
        Vxdatares2D->setFile(Vxresidualfile);
        Vxdatares2D->createEmpty(Vxdata2D->getNtrace());

        Vzdatamod2D = std::make_shared<rockseis::Data2D<float>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
        Vzdatamod2D->setFile(Vzmodelledfile);
        Vzdatamod2D->createEmpty(Vxdata2D->getNtrace());
        Vzdatares2D = std::make_shared<rockseis::Data2D<float>>(1, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
        Vzdatares2D->setFile(Vzresidualfile);
        Vzdatares2D->createEmpty(Vxdata2D->getNtrace());
       
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
            vpgrad = std::make_shared<rockseis::Image2D<float>>(Vpgradfile, gmodel, nhx, nhz);
            vpgrad->createEmpty();
        }
        if(Vsgrad){

            vsgrad = std::make_shared<rockseis::Image2D<float>>(Vsgradfile, gmodel, nhx, nhz);
            vsgrad->createEmpty();
        }
        if(Rhograd){

            rhograd = std::make_shared<rockseis::Image2D<float>>(Rhogradfile, gmodel, nhx, nhz);
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
        std::shared_ptr<rockseis::FwiElastic2D<float>> fwi;
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
                Vxdata2D = Sort->get2DGather(work.id);
                size_t ntr = Vxdata2D->getNtrace();

                Sort->setDatafile(Vzrecordfile);
                Vzdata2D = Sort->get2DGather(work.id);

                lmodel = gmodel->getLocal(Vxdata2D, apertx, SMAP);

                // Read wavelet data, set shot coordinates and make a map
                source->read();
                source->copyCoords(Vxdata2D);
                source->makeMap(lmodel->getGeom(), SMAP);

                //Setting sourcetype 
                switch(stype){
                    case 0:
                        source->setField(PRESSURE);
                        break;
                    case 1:
                        source->setField(VX);
                        break;
                    case 3:
                        source->setField(VZ);
                        break;
                    default:
                        rs_error("Unknown source type: ", std::to_string(stype));
                        break;
                }

                // Interpolate shot
                Vxdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vxdata2D, Vxdata2Di);
                Vxdata2Di->makeMap(lmodel->getGeom(), GMAP);
                Vxdata2Di->setField(rockseis::VX);

                Vzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                interp->interp(Vzdata2D, Vzdata2Di);
                Vzdata2Di->makeMap(lmodel->getGeom(), GMAP);
                Vzdata2Di->setField(rockseis::VZ);

                // Create fwi object
                fwi = std::make_shared<rockseis::FwiElastic2D<float>>(lmodel, source, Vxdata2Di, Vzdata2Di, order, snapinc);

                // Create modelled and residual data objects 
                Vxdatamod2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatamod2D->copyCoords(Vxdata2D);
                Vxdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatamod2D->setField(rockseis::VX);
                fwi->setDatamodVx(Vxdatamod2D);
                Vxdatares2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vxdatares2D->copyCoords(Vxdata2D);
                Vxdatares2D->makeMap(lmodel->getGeom(), GMAP);
                Vxdatares2D->setField(rockseis::VX);
                fwi->setDataresVx(Vxdatares2D);

                Vzdatamod2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatamod2D->copyCoords(Vzdata2D);
                Vzdatamod2D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatamod2D->setField(rockseis::VZ);
                fwi->setDatamodVz(Vzdatamod2D);
                Vzdatares2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                Vzdatares2D->copyCoords(Vzdata2D);
                Vzdatares2D->makeMap(lmodel->getGeom(), GMAP);
                Vzdatares2D->setField(rockseis::VZ);
                fwi->setDataresVz(Vzdatares2D);

                // Setting misfit type
                fwi->setMisfit_type(fwimisfit);

                // Creating gradient objects and setting them up in fwi class
                if(Vpgrad){
                    vpgrad = std::make_shared<rockseis::Image2D<float>>(Vpgradfile + "-" + std::to_string(work.id), lmodel, nhx, nhz);
                fwi->setVpgrad(vpgrad);
                }
                if(Vsgrad){
                    vsgrad = std::make_shared<rockseis::Image2D<float>>(Vsgradfile + "-" + std::to_string(work.id), lmodel, nhx, nhz);
                fwi->setVsgrad(vsgrad);
                }
                if(Rhograd){
                    rhograd = std::make_shared<rockseis::Image2D<float>>(Rhogradfile + "-" + std::to_string(work.id), lmodel, nhx, nhz);
                fwi->setRhograd(rhograd);
                }
   

                if(Wavgrad){
                    wavgrad = std::make_shared<rockseis::Data2D<float>>(source->getNtrace(), source->getNt(), source->getDt(), 0.0);
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
                Vxdatamod2Di = std::make_shared<rockseis::Data2D<float>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vxdatamod2Di->setFile(Vxmodelledfile);
                interp->interp(Vxdatamod2D, Vxdatamod2Di);
                Sort->put2DGather(Vxdatamod2Di, work.id);

                Vxdatares2Di = std::make_shared<rockseis::Data2D<float>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vxdatares2Di->setFile(Vxresidualfile);
                interp->interp(Vxdatares2D, Vxdatares2Di);
                Sort->put2DGather(Vxdatares2Di, work.id);

                Vzdatamod2Di = std::make_shared<rockseis::Data2D<float>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vzdatamod2Di->setFile(Vzmodelledfile);
                interp->interp(Vzdatamod2D, Vzdatamod2Di);
                Sort->put2DGather(Vzdatamod2Di, work.id);

                Vzdatares2Di = std::make_shared<rockseis::Data2D<float>>(ntr, Vxdata2D->getNt(), Vxdata2D->getDt(), Vxdata2D->getOt());
                Vzdatares2Di->setFile(Vzresidualfile);
                interp->interp(Vzdatares2D, Vzdatares2Di);
                Sort->put2DGather(Vzdatares2Di, work.id);

                // Reset all classes
                Vxdata2D.reset();
                Vxdata2Di.reset();
                Vxdatamod2D.reset();
                Vxdatamod2Di.reset();
                Vxdatares2D.reset();
                Vxdatares2Di.reset();
                Vzdata2D.reset();
                Vzdata2Di.reset();
                Vzdatamod2D.reset();
                Vzdatamod2Di.reset();
                Vzdatares2D.reset();
                Vzdatares2Di.reset();
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

