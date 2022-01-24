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
   MPIdomaindecomp mpi = MPIdomaindecomp(&argc,&argv);
   if(mpi.getNrank() < 2){
      rs_error("This is a parallel program, it must run with at least 2 processors, use mpirun.");
   }

   if(argc < 2){
      if(mpi.getRank() == 0){
         PRINT_DOC(# MPI 3d elastic modelling default configuration file);
         PRINT_DOC();
         PRINT_DOC(# Domain decomposition parameter);
         PRINT_DOC(        ndomain0 = "1";  # Number of domains along x direction to split the model into);
         PRINT_DOC(        ndomain1 = "1";  # Number of domains along y direction to split the model into);
         PRINT_DOC(        ndomain2 = "1";  # Number of domains along z direction to split the model into);
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(        freesurface = "true";  # True if free surface should be on);
         PRINT_DOC(            order = "8";  # Order of finite difference stencil 2-8);
         PRINT_DOC(            lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
         PRINT_DOC(            source_type = "0"; # Source type 0 - pressure. 1 for Vx. 2 for Vy. 3 for Vz.);
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
         PRINT_DOC(        C11 = "C11.rss";);
         PRINT_DOC(        C12 = "C12.rss";);
         PRINT_DOC(        C13 = "C13.rss";);
         PRINT_DOC(        C22 = "C22.rss";);
         PRINT_DOC(        C23 = "C23.rss";);
         PRINT_DOC(        C33 = "C33.rss";);
         PRINT_DOC(        C44 = "C44.rss";);
         PRINT_DOC(        C55 = "C55.rss";);
         PRINT_DOC(        C66 = "C66.rss";);
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
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   int lpml;
   bool fs;
   int order;
   int ndomain0;
   int ndomain1;
   int ndomain2;
   int snapinc;
   float apertx;
   float aperty;
   float dtrec;
   int stype;
   std::string Surveyfile;
   std::string Waveletfile;
   std::string C11file;
   std::string C12file;
   std::string C13file;
   std::string C22file;
   std::string C23file;
   std::string C33file;
   std::string C44file;
   std::string C55file;
   std::string C66file;
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
   std::shared_ptr<rockseis::ModelOrtho3D<float>> lmodel;

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
   if(Inpar->getPar("ndomain0", &ndomain0) == INPARSE_ERR) status = true;
   if(Inpar->getPar("ndomain1", &ndomain1) == INPARSE_ERR) status = true;
   if(Inpar->getPar("ndomain2", &ndomain2) == INPARSE_ERR) status = true;
   if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("source_type", &stype) == INPARSE_ERR) status = true;
   if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C11", &C11file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C12", &C12file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C13", &C13file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C22", &C22file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C23", &C23file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C33", &C33file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C44", &C44file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C55", &C55file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C66", &C66file) == INPARSE_ERR) status = true;
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

   // Setup Domain decomposition
   mpi.setNdomain(ndomain0*ndomain1*ndomain2);
   mpi.splitDomains();

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Surveyfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelOrtho3D<float>> gmodel (new rockseis::ModelOrtho3D<float>(C11file, C12file, C13file, C22file, C23file, C33file, C44file, C55file, C66file, Rhofile, lpml ,fs));

   // Create a data class for the source wavelet
   std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

   // Create an interpolation class
   std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

   // Compute record length in samples
   size_t ntrec; 
   ntrec = (size_t) rintf((source->getNt()-1)*source->getDt()/dtrec + 1);

   if(mpi.getRank() == 0) {
      // Master
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
      std::shared_ptr<rockseis::ModellingOrtho3D<float>> modelling;
      while(1) {
         if(!mpi.ifActive()){
            break;
         }
         workModeling_t work = mpi.receiveWork();

         if(work.MPItag == MPI_TAG_DIE) {
            break;
         }

         if(work.MPItag == MPI_TAG_NO_WORK) {
            mpi.sendNoWork(mpi.getMasterComm(), 0);
         }
         else {
            // Do some work
            Sort->readKeymap();
            Sort->readSortmap();

            Shotgeom = Sort->get3DGather(work.id);
            size_t ntr = Shotgeom->getNtrace();
            lmodel = gmodel->getDomainmodel(Shotgeom, apertx, aperty, SMAP, mpi.getDomainrank(), ndomain0,ndomain1,ndomain2, order);
            (lmodel->getDomain())->setMpi(&mpi);

            // Read wavelet data, set shot coordinates and make a map
            source->read();
            source->copyCoords(Shotgeom);

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
            source->makeMap(lmodel->getGeom(), SMAP);

            modelling = std::make_shared<rockseis::ModellingOrtho3D<float>>(lmodel, source, order, snapinc);

            // Set logfile
            modelling->setLogfile("log.txt-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));

            // Setting Snapshot file 
            if(Psnap){
               modelling->setSnapP(Psnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Vxsnap){
               modelling->setSnapVx(Vxsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Vysnap){
               modelling->setSnapVy(Vysnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Vzsnap){
               modelling->setSnapVz(Vzsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }

            // Setting Record
            if(Precord){
               Pdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Pdata3D->setField(rockseis::PRESSURE);
               // Copy geometry to Data
               Pdata3D->copyCoords(Shotgeom);
               Pdata3D->makeMap(lmodel->getGeom(),SMAP,0,0,0);
               Pdata3D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(1),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(1),(lmodel->getDomain())->getPadh(2));
               modelling->setRecP(Pdata3D);
            }
            if(Vxrecord){
               Vxdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Vxdata3D->setField(rockseis::VX);
               // Copy geometry to Data
               Vxdata3D->copyCoords(Shotgeom);
               Vxdata3D->makeMap(lmodel->getGeom(),SMAP,0,0,0);
               Vxdata3D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(1),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(1),(lmodel->getDomain())->getPadh(2));
               modelling->setRecVx(Vxdata3D);
            }
            if(Vyrecord){
               Vydata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Vydata3D->setField(rockseis::VY);
               // Copy geometry to Data
               Vydata3D->copyCoords(Shotgeom);
               Vydata3D->makeMap(lmodel->getGeom(),SMAP,0,0,0);
               Vydata3D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(1),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(1),(lmodel->getDomain())->getPadh(2));
               modelling->setRecVy(Vydata3D);
            }
            if(Vzrecord){
               Vzdata3D = std::make_shared<rockseis::Data3D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Vzdata3D->setField(rockseis::VZ);
               // Copy geometry to Data
               Vzdata3D->copyCoords(Shotgeom);
               Vzdata3D->makeMap(lmodel->getGeom(),SMAP,0,0,0);
               Vzdata3D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(1),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(1),(lmodel->getDomain())->getPadh(2));
               modelling->setRecVz(Vzdata3D);
            }

            // Run modelling 
            modelling->run();

            // Output record
            if(Precord){
               Pdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
               Pdata3Di->setFile(Precordfile);
               interp->interp(Pdata3D, Pdata3Di);
               Sort->put3DGather(Pdata3Di, work.id, (Pdata3D->getGeom())->getGmap());
            }
            if(Vxrecord){
               Vxdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
               Vxdata3Di->setFile(Vxrecordfile);
               interp->interp(Vxdata3D, Vxdata3Di);
               Sort->put3DGather(Vxdata3Di, work.id, (Vxdata3D->getGeom())->getGmap());
            }
            if(Vyrecord){
               Vydata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
               Vydata3Di->setFile(Vyrecordfile);
               interp->interp(Vydata3D, Vydata3Di);
               Sort->put3DGather(Vydata3Di, work.id, (Vydata3D->getGeom())->getGmap());
            }
            if(Vzrecord){
               Vzdata3Di = std::make_shared<rockseis::Data3D<float>>(ntr, ntrec, dtrec, 0.0);
               Vzdata3Di->setFile(Vzrecordfile);
               interp->interp(Vzdata3D, Vzdata3Di);
               Sort->put3DGather(Vzdata3Di, work.id, (Vzdata3D->getGeom())->getGmap());
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
            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }
   }
}

