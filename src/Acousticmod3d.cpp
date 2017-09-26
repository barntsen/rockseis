#include <iostream>
#include <fstream>
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

using namespace rockseis;

int main(int argc, char** argv) {
    if(argc < 2){
        PRINT_DOC(# 3d acoustic modelling default configuration file);
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
        PRINT_DOC(            Axrecord = "false";);
        PRINT_DOC(            Ayrecord = "false";);
        PRINT_DOC(        Azrecord = "false";);
        PRINT_DOC(        Psnap = "false";);
        PRINT_DOC(        Axsnap = "false";);
        PRINT_DOC(        Aysnap = "false";);
        PRINT_DOC(        Azsnap = "false";);
        PRINT_DOC();
        PRINT_DOC(# Files);
        PRINT_DOC(        Vp = "Vp3d.rss";);
        PRINT_DOC(        Rho = "Rho3d.rss";);
        PRINT_DOC(        Wavelet = "Wav3d.rss";);
        PRINT_DOC(        Survey = "3DSurvey.rss";);
        PRINT_DOC(        Precordfile = "Pshot.rss";);
        PRINT_DOC(        Axrecordfile = "Axshot.rss";);
        PRINT_DOC(        Ayrecordfile = "Ayshot.rss";);
        PRINT_DOC(        Azrecordfile = "Azshot.rss";);
        PRINT_DOC(        Psnapfile = "Psnaps.rss";);
        PRINT_DOC(        Axsnapfile = "Axsnaps.rss";);
        PRINT_DOC(        Aysnapfile = "Aysnaps.rss";);
        PRINT_DOC(        Azsnapfile = "Azsnaps.rss";);
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
    std::string Surveyfile;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    bool Psnap=0, Precord=0;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data3D<float>> Pdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Pdata3Di;

    bool Axsnap=0, Axrecord=0;
    std::string Axsnapfile;
    std::string Axrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Axdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Axdata3Di;

    bool Aysnap=0, Ayrecord=0;
    std::string Aysnapfile;
    std::string Ayrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Aydata3D;
    std::shared_ptr<rockseis::Data3D<float>> Aydata3Di;

    bool Azsnap=0, Azrecord=0;
    std::string Azsnapfile;
    std::string Azrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Azdata3D;
    std::shared_ptr<rockseis::Data3D<float>> Azdata3Di;

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
    if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Psnap", &Psnap) == INPARSE_ERR) status = true;
    if(Psnap){
        if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Axsnap", &Axsnap) == INPARSE_ERR) status = true;
    if(Axsnap){
        if(Inpar->getPar("Axsnapfile", &Axsnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Aysnap", &Aysnap) == INPARSE_ERR) status = true;
    if(Aysnap){
        if(Inpar->getPar("Aysnapfile", &Aysnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Azsnap", &Azsnap) == INPARSE_ERR) status = true;
    if(Azsnap){
        if(Inpar->getPar("Azsnapfile", &Azsnapfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Precord", &Precord) == INPARSE_ERR) status = true;
    if(Precord){
        if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Axrecord", &Axrecord) == INPARSE_ERR) status = true;
    if(Axrecord){
        if(Inpar->getPar("Axrecordfile", &Axrecordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Azrecord", &Azrecord) == INPARSE_ERR) status = true;
    if(Azrecord){
        if(Inpar->getPar("Azrecordfile", &Azrecordfile) == INPARSE_ERR) status = true;
    }

    if(status == true){
        rs_error("Program terminated due to input errors.");
    }

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Surveyfile);

    // Create a global model class
    std::shared_ptr<rockseis::ModelAcoustic3D<float>> gmodel (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));
    // Create a local model class
    std::shared_ptr<rockseis::ModelAcoustic3D<float>> lmodel (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));

    // Create a data class for the source wavelet
    std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

    // Create an interpolation class
    std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

    // Compute record length in samples
    size_t ntrec; 
    ntrec = (size_t) rintf((source->getNt()-1)*source->getDt()/dtrec + 1);

    Sort->createShotmap(Surveyfile); 
    Sort->writeKeymap();
    Sort->writeSortmap();

    // Get number of shots
    size_t ngathers =  Sort->getNensemb();

    if(Precord){
        // Create an empty data file
        Sort->createEmptydataset(Precordfile, ntrec, dtrec, 0.0);
    }

    if(Axrecord){
        // Create an empty data file
        Sort->createEmptydataset(Axrecordfile, ntrec, dtrec, 0.0);
    }
    if(Ayrecord){
        // Create an empty data file
        Sort->createEmptydataset(Ayrecordfile, ntrec, dtrec, 0.0);
    }
    if(Azrecord){
        // Create an empty data file
        Sort->createEmptydataset(Azrecordfile, ntrec, dtrec, 0.0);
    }

    std::shared_ptr<rockseis::Data3D<float>> Shotgeom;
    std::shared_ptr<rockseis::ModellingAcoustic3D<float>> modelling;
    // Loop over shot gathers
    for(unsigned long int i=0; i<ngathers; i++) {
        Sort->readKeymap();
        Sort->readSortmap();

        Shotgeom = Sort->get3DGather(i);
        size_t ntr = Shotgeom->getNtrace();
        lmodel = gmodel->getLocal(Shotgeom, apertx, aperty, SMAP);

        // Read wavelet data, set shot and receiver coordinates and make a map
        source->read();
        source->copyCoords(Shotgeom);
        source->makeMap(lmodel->getGeom(), SMAP);

        modelling = std::make_shared<rockseis::ModellingAcoustic3D<float>>(lmodel, source, order, snapinc);

        // Set logfile
        modelling->setLogfile("log.txt-" + std::to_string(i));

        // Setting Snapshot file 
        if(Psnap){
            modelling->setSnapP(Psnapfile + "-" + std::to_string(i));
        }
        if(Axsnap){
            modelling->setSnapAx(Axsnapfile + "-" + std::to_string(i));
        }
        if(Aysnap){
            modelling->setSnapAy(Aysnapfile + "-" + std::to_string(i));
        }
        if(Azsnap){
            modelling->setSnapAz(Azsnapfile + "-" + std::to_string(i));
        }

        // Setting Record
        if(Precord){
            Pdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Pdata3D->setField(rockseis::PRESSURE);
            // Copy geometry to Data
            Pdata3D->copyCoords(Shotgeom);
            Pdata3D->makeMap(lmodel->getGeom());
            modelling->setRecP(Pdata3D);
        }
        if(Axrecord){
            Axdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Axdata3D->setField(rockseis::VX);
            // Copy geometry to Data
            Axdata3D->copyCoords(Shotgeom);
            Axdata3D->makeMap(lmodel->getGeom());
            modelling->setRecAx(Axdata3D);
        }
        if(Ayrecord){
            Aydata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Aydata3D->setField(rockseis::VY);
            // Copy geometry to Data
            Aydata3D->copyCoords(Shotgeom);
            Aydata3D->makeMap(lmodel->getGeom());
            modelling->setRecAy(Aydata3D);
        }
        if(Azrecord){
            Azdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Azdata3D->setField(rockseis::VZ);
            // Copy geometry to Data
            Azdata3D->copyCoords(Shotgeom);
            Azdata3D->makeMap(lmodel->getGeom());
            modelling->setRecAz(Azdata3D);
        }

        // Stagger model
        lmodel->staggerModels();

        // Run modelling 
        modelling->run();

        // Output record
        if(Precord){
            Pdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Pdata3Di->setFile(Precordfile);
            interp->interp(Pdata3D, Pdata3Di);
            Sort->put3DGather(Pdata3Di, i);
        }
        if(Axrecord){
            Axdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Axdata3Di->setFile(Axrecordfile);
            interp->interp(Axdata3D, Axdata3Di);
            Sort->put3DGather(Axdata3Di, i);
        }
        if(Ayrecord){
            Aydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Aydata3Di->setFile(Ayrecordfile);
            interp->interp(Aydata3D, Aydata3Di);
            Sort->put3DGather(Aydata3Di, i);
        }
        if(Azrecord){
            Azdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Azdata3Di->setFile(Azrecordfile);
            interp->interp(Azdata3D, Azdata3Di);
            Sort->put3DGather(Azdata3Di, i);
        }

        // Reset all classes
        Shotgeom.reset();
        lmodel.reset();
        modelling.reset();
        if(Precord){
            Pdata3D.reset();
            Pdata3Di.reset();
        }
        if(Axrecord){
            Axdata3D.reset();
            Axdata3Di.reset();
        }
        if(Ayrecord){
            Aydata3D.reset();
            Aydata3Di.reset();
        }
        if(Azrecord){
            Azdata3D.reset();
            Azdata3Di.reset();
        }

    }
}
