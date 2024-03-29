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
        PRINT_DOC(# 3d elastic modelling default configuration file);
        PRINT_DOC();
        PRINT_DOC(# Modelling parameters);
        PRINT_DOC(        freesurface = "true";  # True if free surface should be on);
        PRINT_DOC(            order = "8";  # Order of finite difference stencil 2-8);
        PRINT_DOC(            lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
        PRINT_DOC(            source_type = "PRESSURE"; # Source type );
        PRINT_DOC(            snapinc = "10"; # Snap interval in multiples of modelling interval);
        PRINT_DOC(            dtrec = "4e-3"; # Recording interval in seconds);
        PRINT_DOC(            apertx = "900"; # Aperture for local model (source is in the middle));
        PRINT_DOC(            aperty = "900"; # Aperture for local model (source is in the middle));
        PRINT_DOC();
        PRINT_DOC(# Booleans);
        PRINT_DOC(            Precord = "true";  # Set these to true if recording or snapshoting is to be made.);
        PRINT_DOC(            Vxrecord = "true";);
        PRINT_DOC(            Vyrecord = "true";);
        PRINT_DOC(        Vzrecord = "true";);
        PRINT_DOC(        Psnap = "false";);
        PRINT_DOC(        Vxsnap = "false";);
        PRINT_DOC(        Vysnap = "false";);
        PRINT_DOC(        Vzsnap = "false";);
        PRINT_DOC();
        PRINT_DOC(# Files);
        PRINT_DOC(        Vp = "Vp3d.rss";);
        PRINT_DOC(        Vs = "Vs3d.rss";);
        PRINT_DOC(        Rho = "Rho3d.rss";);
        PRINT_DOC(        Wavelet = "Wav3d.rss";);
        PRINT_DOC(        Survey = "3DSurvey.rss";);
        PRINT_DOC(        Precordfile = "Pshot.rss";);
        PRINT_DOC(        Vxrecordfile = "Vxshot.rss";);
        PRINT_DOC(        Vyrecordfile = "Vyshot.rss";);
        PRINT_DOC(        Vzrecordfile = "Vzshot.rss";);
        PRINT_DOC(        Psnapfile = "Psnaps.rss";);
        PRINT_DOC(        Vxsnapfile = "Vxsnaps.rss";);
        PRINT_DOC(        Vysnapfile = "Vysnaps.rss";);
        PRINT_DOC(        Vzsnapfile = "Vzsnaps.rss";);
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
    int stype;
    std::string Surveyfile;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
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

    // Create a local model class pointer
    std::shared_ptr<rockseis::ModelElastic3D<float>> lmodel;

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
    if(Inpar->getPar("source_type", &stype) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
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
    if(Inpar->getPar("Vyrecord", &Vyrecord) == INPARSE_ERR) status = true;
    if(Vyrecord){
        if(Inpar->getPar("Vyrecordfile", &Vyrecordfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Vzrecord", &Vzrecord) == INPARSE_ERR) status = true;
    if(Vzrecord){
        if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
    }

    if(status == true){
        rs_error("Program terminated due to input errors.");
    }

    // Create a sort class
    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(Surveyfile);

    // Create a global model class
    std::shared_ptr<rockseis::ModelElastic3D<float>> gmodel (new rockseis::ModelElastic3D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));

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

    if(Vxrecord){
        // Create an empty data file
        Sort->createEmptydataset(Vxrecordfile, ntrec, dtrec, 0.0);
    }
    if(Vyrecord){
        // Create an empty data file
        Sort->createEmptydataset(Vyrecordfile, ntrec, dtrec, 0.0);
    }
    if(Vzrecord){
        // Create an empty data file
        Sort->createEmptydataset(Vzrecordfile, ntrec, dtrec, 0.0);
    }

    std::shared_ptr<rockseis::Data3D<float>> Shotgeom;
    std::shared_ptr<rockseis::ModellingElastic3D<float>> modelling;

    for(unsigned long int i=0; i<ngathers; i++) {
        Sort->readKeymap();
        Sort->readSortmap();

        Shotgeom = Sort->get3DGather(i);
        size_t ntr = Shotgeom->getNtrace();
        lmodel = gmodel->getLocal(Shotgeom, apertx, aperty, SMAP);

        // Read wavelet data, set shot coordinates and make a map
        source->read();
        source->copyCoords(Shotgeom);
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

        modelling = std::make_shared<rockseis::ModellingElastic3D<float>>(lmodel, source, order, snapinc);

        // Set logfile
        modelling->setLogfile("log.txt-" + std::to_string(i));

        // Setting Snapshot file 
        if(Psnap){
            modelling->setSnapP(Psnapfile + "-" + std::to_string(i));
        }
        if(Vxsnap){
            modelling->setSnapVx(Vxsnapfile + "-" + std::to_string(i));
        }
        if(Vysnap){
            modelling->setSnapVy(Vysnapfile + "-" + std::to_string(i));
        }
        if(Vzsnap){
            modelling->setSnapVz(Vzsnapfile + "-" + std::to_string(i));
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
        if(Vxrecord){
            Vxdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Vxdata3D->setField(rockseis::VX);
            // Copy geometry to Data
            Vxdata3D->copyCoords(Shotgeom);
            Vxdata3D->makeMap(lmodel->getGeom());
            modelling->setRecVx(Vxdata3D);
        }
        if(Vyrecord){
            Vydata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            Vydata3D->setField(rockseis::VY);
            // Copy geometry to Data
            Vydata3D->copyCoords(Shotgeom);
            Vydata3D->makeMap(lmodel->getGeom());
            modelling->setRecVy(Vydata3D);
        }
        if(Vzrecord){
            Vzdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
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
            Pdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Pdata3Di->setFile(Precordfile);
            interp->interp(Pdata3D, Pdata3Di);
            Sort->put3DGather(Pdata3Di, i);
        }
        if(Vxrecord){
            Vxdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Vxdata3Di->setFile(Vxrecordfile);
            interp->interp(Vxdata3D, Vxdata3Di);
            Sort->put3DGather(Vxdata3Di, i);
        }
        if(Vyrecord){
            Vydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Vydata3Di->setFile(Vyrecordfile);
            interp->interp(Vydata3D, Vydata3Di);
            Sort->put3DGather(Vydata3Di, i);
        }
        if(Vzrecord){
            Vzdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
            Vzdata3Di->setFile(Vzrecordfile);
            interp->interp(Vzdata3D, Vzdata3Di);
            Sort->put3DGather(Vzdata3Di, i);
        }

        // Reset all classes
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

    }
}

