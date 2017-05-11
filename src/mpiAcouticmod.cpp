#include <iostream>
#include "parallel.h"
#include <mpi.h>
#include "model.h"
#include "modelling.h"
#include "pml.h"
#include "waves.h"
#include "utils.h"
#include "der.h"
#include "sort.h"
#include "data.h"
#include "file.h"
#include <memory>
#include <fstream>
#include <math.h>
#include <config4cpp/Configuration.h>

using namespace rockseis;

int main(int argc, char** argv) {
	// Initializing MPI
	MPImodeling mpi = MPImodeling(&argc,&argv);
    if(mpi.getNrank() < 2){
        rs_error("This is a parallel program, it must run with at least 2 processors, use mpirun.");
    }

	bool status;
	/* General input parameters */
	int lpml=0;
	bool fs=0;
	int order=0;
	int snapinc=0;
    int dim;
    float apertx=900;
    float aperty=0;
    std::string Surveyfile;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    bool Psnap=0, Precord=0;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data2D<float>> Pdata2D;
    std::shared_ptr<rockseis::Data3D<float>> Pdata3D;

    bool Axsnap=0, Axrecord=0;
    std::string Axsnapfile;
    std::string Axrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> Axdata2D;
    std::shared_ptr<rockseis::Data3D<float>> Axdata3D;

    bool Aysnap=0, Ayrecord=0;
    std::string Aysnapfile;
    std::string Ayrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Aydata3D;

    bool Azsnap=0, Azrecord=0;
    std::string Azsnapfile;
    std::string Azrecordfile;
    std::shared_ptr<rockseis::Data2D<float>> Azdata2D;
    std::shared_ptr<rockseis::Data3D<float>> Azdata3D;

    /* Get parameters from configuration file */
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "mpiacumod.cfg";

    status = 0;
    try {
        cfg->parse(configFile);
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        return 1;
    }

    try {
        lpml = cfg->lookupInt(scope, "lpml");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        order = cfg->lookupInt(scope, "order");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        snapinc = cfg->lookupInt(scope, "snapinc");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        dim = cfg->lookupInt(scope, "dim");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        fs = cfg->lookupBoolean(scope, "freesurface");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Vpfile = cfg->lookupString(scope, "Vp");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        Rhofile = cfg->lookupString(scope, "Rho");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Waveletfile = cfg->lookupString(scope, "wavelet");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Surveyfile = cfg->lookupString(scope, "survey");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Psnap = cfg->lookupBoolean(scope, "Psnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        apertx = cfg->lookupFloat(scope, "apertx");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(dim == 3){
        try {
            aperty = cfg->lookupFloat(scope, "aperty");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    if(Psnap){
        try {
            Psnapfile = cfg->lookupString(scope, "Psnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Precord = cfg->lookupBoolean(scope, "Precord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Precord){
        try {
            Precordfile = cfg->lookupString(scope, "Precordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
   try {
        Axsnap = cfg->lookupBoolean(scope, "Axsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Axsnap){
        try {
            Axsnapfile = cfg->lookupString(scope, "Axsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Axrecord = cfg->lookupBoolean(scope, "Axrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Axrecord){
        try {
            Axrecordfile = cfg->lookupString(scope, "Axrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
   }
   try {
        Aysnap = cfg->lookupBoolean(scope, "Aysnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Aysnap){
        try {
            Aysnapfile = cfg->lookupString(scope, "Aysnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Ayrecord = cfg->lookupBoolean(scope, "Ayrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Ayrecord){
        try {
            Ayrecordfile = cfg->lookupString(scope, "Ayrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }

   try {
        Azsnap = cfg->lookupBoolean(scope, "Azsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Azsnap){
        try {
            Azsnapfile = cfg->lookupString(scope, "Azsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Azrecord = cfg->lookupBoolean(scope, "Azrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Azrecord){
        try {
            Azrecordfile = cfg->lookupString(scope, "Azrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }

	// Destroy cfg
	cfg->destroy();

	if(status == 1){
		rs_error("Program terminated due to input errors.");
	}

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Surveyfile);
	

    // Create a global model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> gmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Rhofile, lpml ,fs));
    // Create a local model class
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> lmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Rhofile, lpml ,fs));

	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(Waveletfile));
	if(mpi.getRank() == 0) {
		// Master
        Sort->createShotmap(Surveyfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // Get number of shots
        int ngathers =  Sort->getNensemb();
        size_t ntraces = Sort->getNtraces();

		// Create work queue
		for(unsigned long int i=0; i<ngathers; i++) {
			// Work struct
			std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
			mpi.addWork(work);
		}

		// Print work queue
		std::cerr << "Work queue before parallelization" << std::endl;
		mpi.printWork();

		// Perform work in parallel
		mpi.performWork();
	
		// Print work queue
		std::cerr << "Work queue after parallelization" << std::endl;
		mpi.printWork();
    }
    else {
        /* Slave */
        std::shared_ptr<rockseis::Data2D<float>> OneShot;
        std::shared_ptr<rockseis::ModellingAcoustic2D<float>> modelling;
        rockseis::Point2D<float> *scoords;
        float xs,zs; // Source coordinates
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

                OneShot = Sort->get2DGather(work.id);
                size_t ntr = OneShot->getNtrace();
                lmodel = gmodel->getLocal(OneShot, apertx, SMAP);
                scoords = (OneShot->getGeom())->getScoords();
                xs = scoords[0].x;
                zs = scoords[0].y;

                scoords = (source->getGeom())->getScoords();
                scoords[0].x = xs;
                scoords[0].y = zs;

                // Read wavelet data and coordinates and make a map
                source->read();
                source->makeMap(lmodel->getGeom());

                modelling = std::make_shared<rockseis::ModellingAcoustic2D<float>>(lmodel, source, order, snapinc);

                // Setting Snapshot file 
                if(Psnap){
                    modelling->setSnapP(Psnapfile);
                }
                if(Axsnap){
                    modelling->setSnapAx(Axsnapfile);
                }
                if(Azsnap){
                    modelling->setSnapAz(Azsnapfile);
                }

                // Setting Record
                if(Precord){
                    Pdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
                    Pdata2D->setField(rockseis::PRESSURE);
                    std:: string datafile = std::to_string(work.id) + "_" + Precordfile;
                    Pdata2D->setFile(datafile);
                    // Copy geometry to Data
                    Pdata2D->copyCoords(OneShot);
                    Pdata2D->makeMap(lmodel->getGeom());
                    modelling->setRecP(Pdata2D);
                }

                // Stagger model
                lmodel->staggerModels();

                // Run modelling 
                modelling->run();

                // Output record
                if(Precord){
                    Pdata2D->write();
                }


                // Reset all classes
                OneShot.reset();
                lmodel.reset();
                modelling.reset();
                if(Precord){
                    Pdata2D.reset();
                }
                work.status = WORK_FINISHED;

                // Send result back
                mpi.sendResult(work);		
            }
        }
    }
}

