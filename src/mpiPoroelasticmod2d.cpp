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
         PRINT_DOC(# MPI 2d poroelastic modelling default configuration file);
         PRINT_DOC();
         PRINT_DOC(# Domain decomposition parameter);
         PRINT_DOC(        ndomain0 = "1";  # Number of domains along x direction to split the model into);
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(        freesurface = "true";  # True if free surface should be on);
         PRINT_DOC(            order = "8";  # Order of finite difference stencil 2-8);
         PRINT_DOC(            lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
         PRINT_DOC(            source_type = "0"; # Source type 0 - pressure. 1 for Vx. 3 for Vz.);
         PRINT_DOC(            snapinc = "10"; # Snap interval in multiples of modelling interval);
         PRINT_DOC(            dtrec = "4e-3"; # Recording interval in seconds);
         PRINT_DOC(            apertx = "900"; # Aperture for local model (source is in the middle));
         PRINT_DOC(            f0 = "20"; # Dominant frequency of the source (For PML parameters));
         PRINT_DOC();
         PRINT_DOC(# Booleans);
         PRINT_DOC(            Precord = "true";  # Set these to true if recording or snapshoting is to be made.);
         PRINT_DOC(            Qxrecord = "true";);
         PRINT_DOC(        Qzrecord = "true";);
         PRINT_DOC(        Psnap = "false";);
         PRINT_DOC(        Qxsnap = "false";);
         PRINT_DOC(        Qzsnap = "false";);
         PRINT_DOC(            Szzrecord = "true";  # Set these to true if recording or snapshoting is to be made.);
         PRINT_DOC(            Vxrecord = "true";);
         PRINT_DOC(        Vzrecord = "true";);
         PRINT_DOC(        Szzsnap = "false";);
         PRINT_DOC(        Vxsnap = "false";);
         PRINT_DOC(        Vzsnap = "false";);
//         PRINT_DOC();
//         PRINT_DOC(# Input acoustic model files);
//         PRINT_DOC(        Rho_acu = "Rho_acu2d.rss";);
//         PRINT_DOC(        Vp_acu = "Vp_acu2d.rss";);
         PRINT_DOC();
         PRINT_DOC(# Input poroelastic model files);
         PRINT_DOC(        Rho = "Rho2d.rss";);
         PRINT_DOC(        Rhof = "Rhof2d.rss";);
         PRINT_DOC(        Por = "Por2d.rss";);
         PRINT_DOC(        Kd = "Kd2d.rss";);
         PRINT_DOC(        Ks = "Ks2d.rss";);
         PRINT_DOC(        Kf = "Kf2d.rss";);
         PRINT_DOC(        Mu = "Mu2d.rss";);
         PRINT_DOC(        Mob = "Mob2d.rss";);
         PRINT_DOC(        Psi = "Psi2d.rss";);
         PRINT_DOC(        Wavelet = "Wav2d.rss";);
         PRINT_DOC(        Survey = "2DSurvey.rss";);
         PRINT_DOC();
         PRINT_DOC(# Output files);
         PRINT_DOC(        Precordfile = "Pshot.rss";);
         PRINT_DOC(        Qxrecordfile = "Qxshot.rss";);
         PRINT_DOC(        Qzrecordfile = "Qzshot.rss";);
         PRINT_DOC(        Psnapfile = "Psnaps.rss";);
         PRINT_DOC(        Qxsnapfile = "Qxsnaps.rss";);
         PRINT_DOC(        Qzsnapfile = "Qzsnaps.rss";);
         PRINT_DOC(        Szzrecordfile = "Szzshot.rss";);
         PRINT_DOC(        Vxrecordfile = "Vxshot.rss";);
         PRINT_DOC(        Vzrecordfile = "Vzshot.rss";);
         PRINT_DOC(        Psnapfile = "Szzsnaps.rss";);
         PRINT_DOC(        Vxsnapfile = "Vxsnaps.rss";);
         PRINT_DOC(        Vzsnapfile = "Vzsnaps.rss";);
      }
      exit(1);
   }

   bool status;
   /* General input parameters */
   int lpml;
   bool fs;
   int ndomain0;
   int ndomain1;
   int order;
   int snapinc;
   float apertx;
   float f0;
   float dtrec;
   int stype;
//   std::string Acu_rhofile;
//   std::string Acu_vpfile;
   std::string Surveyfile;
   std::string Waveletfile;
   std::string Rhofile;
   std::string Rhoffile;
   std::string Porfile;
   std::string Kdfile;
   std::string Ksfile;
   std::string Kffile;
   std::string Mufile;
   std::string Mobfile;
   std::string Psifile;
   bool Psnap=0, Precord=0;
   std::string Psnapfile;
   std::string Precordfile;
   std::shared_ptr<rockseis::Data2D<float>> Pdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Pdata2Di;

   bool Vxsnap=0, Vxrecord=0;
   std::string Vxsnapfile;
   std::string Vxrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Vxdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Vxdata2Di;

   bool Vzsnap=0, Vzrecord=0;
   std::string Vzsnapfile;
   std::string Vzrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Vzdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Vzdata2Di;

   bool Szzsnap=0, Szzrecord=0;
   std::string Szzsnapfile;
   std::string Szzrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Szzdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Szzdata2Di;

   bool Qxsnap=0, Qxrecord=0;
   std::string Qxsnapfile;
   std::string Qxrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Qxdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Qxdata2Di;

   bool Qzsnap=0, Qzrecord=0;
   std::string Qzsnapfile;
   std::string Qzrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Qzdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Qzdata2Di;



   // Create a local model class
   std::shared_ptr<rockseis::ModelPoroelastic2D<float>> poro_lmodel;
   std::shared_ptr<rockseis::ModelAcoustic2D<float>> acu_lmodel;

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
   if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("source_type", &stype) == INPARSE_ERR) status = true;
   if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
//   if(Inpar->getPar("Vp_acu", &Acu_vpfile) == INPARSE_ERR) status = true;
//   if(Inpar->getPar("Rho_acu", &Acu_rhofile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Rhof", &Rhoffile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Por", &Porfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Kd", &Kdfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Ks", &Ksfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Kf", &Kffile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Mu", &Mufile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Mob", &Mobfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Psi", &Psifile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("f0", &f0) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Szzsnap", &Szzsnap) == INPARSE_ERR) status = true;
   if(Szzsnap){
      if(Inpar->getPar("Szzsnapfile", &Szzsnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vxsnap", &Vxsnap) == INPARSE_ERR) status = true;
   if(Vxsnap){
      if(Inpar->getPar("Vxsnapfile", &Vxsnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vzsnap", &Vzsnap) == INPARSE_ERR) status = true;
   if(Vzsnap){
      if(Inpar->getPar("Vzsnapfile", &Vzsnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Szzrecord", &Szzrecord) == INPARSE_ERR) status = true;
   if(Szzrecord){
      if(Inpar->getPar("Szzrecordfile", &Szzrecordfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vxrecord", &Vxrecord) == INPARSE_ERR) status = true;
   if(Vxrecord){
      if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vzrecord", &Vzrecord) == INPARSE_ERR) status = true;
   if(Vzrecord){
      if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Psnap", &Psnap) == INPARSE_ERR) status = true;
   if(Psnap){
      if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Qxsnap", &Qxsnap) == INPARSE_ERR) status = true;
   if(Qxsnap){
      if(Inpar->getPar("Qxsnapfile", &Qxsnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Qzsnap", &Qzsnap) == INPARSE_ERR) status = true;
   if(Qzsnap){
      if(Inpar->getPar("Qzsnapfile", &Qzsnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Precord", &Precord) == INPARSE_ERR) status = true;
   if(Precord){
      if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Qxrecord", &Qxrecord) == INPARSE_ERR) status = true;
   if(Qxrecord){
      if(Inpar->getPar("Qxrecordfile", &Qxrecordfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Qzrecord", &Qzrecord) == INPARSE_ERR) status = true;
   if(Qzrecord){
      if(Inpar->getPar("Qzrecordfile", &Qzrecordfile) == INPARSE_ERR) status = true;
   }


   if(status == true){
      rs_error("Program terminated due to input errors.");
   }

   ndomain1=1;
   // Setup Domain decomposition
   mpi.setNdomain(ndomain0*ndomain1);
   mpi.splitDomains();

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Surveyfile);

   // Create global model classes
   std::shared_ptr<rockseis::ModelPoroelastic2D<float>> poro_gmodel (new rockseis::ModelPoroelastic2D<float>(Rhofile,Rhoffile,Porfile,Kdfile,Ksfile,Kffile,Mufile,Mobfile,Psifile,lpml,f0,fs));
//   std::shared_ptr<rockseis::ModelAcoustic2D<float>> acu_gmodel (new rockseis::ModelAcoustic2D<float>(Acu_rhofile,Acu_vpfile,lpml,fs));

   // Create a data class for the source wavelet
   std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(Waveletfile));

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

      if(Szzrecord){
         // Create an empty data file
         Sort->createEmptydataset(Szzrecordfile, ntrec, dtrec, 0.0);
      }

      if(Qxrecord){
         // Create an empty data file
         Sort->createEmptydataset(Qxrecordfile, ntrec, dtrec, 0.0);
      }
      if(Qzrecord){
         // Create an empty data file
         Sort->createEmptydataset(Qzrecordfile, ntrec, dtrec, 0.0);
      }

      if(Precord){
         // Create an empty data file
         Sort->createEmptydataset(Precordfile, ntrec, dtrec, 0.0);
      }

      if(Vxrecord){
         // Create an empty data file
         Sort->createEmptydataset(Vxrecordfile, ntrec, dtrec, 0.0);
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
      std::shared_ptr<rockseis::Data2D<float>> Shotgeom;
      std::shared_ptr<rockseis::ModellingPoroelastic2D<float>> modelling;
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

            Shotgeom = Sort->get2DGather(work.id);
            size_t ntr = Shotgeom->getNtrace();
            poro_lmodel = poro_gmodel->getDomainmodel(Shotgeom, apertx, SMAP, mpi.getDomainrank(), ndomain0, ndomain1, order);
            (poro_lmodel->getDomain())->setMpi(&mpi);
//            acu_lmodel = acu_gmodel->getDomainmodel(Shotgeom, apertx, SMAP, mpi.getDomainrank(), ndomain0, ndomain1, order);
//            (acu_lmodel->getDomain())->setMpi(&mpi);

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
               case 3:
                  source->setField(VZ);
                  break;
               default:
                  rs_error("Unknown source type: ", std::to_string(stype));
                  break;
            }
            source->makeMap(poro_lmodel->getGeom(), SMAP);

            modelling = std::make_shared<rockseis::ModellingPoroelastic2D<float>>(poro_lmodel, source, order, snapinc);

            // Set logfile
            modelling->setLogfile("log.txt-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));

            // Setting Snapshot file 
            if(Psnap){
               modelling->setSnapP(Psnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Vxsnap){
               modelling->setSnapVx(Vxsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Vzsnap){
               modelling->setSnapVz(Vzsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Szzsnap){
               modelling->setSnapSzz(Szzsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Qxsnap){
               modelling->setSnapQx(Qxsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }
            if(Qzsnap){
               modelling->setSnapQz(Qzsnapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));
            }

            // Setting Record
            if(Precord){
               Pdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Pdata2D->setField(rockseis::PRESSURE);
               // Copy geometry to Data
               Pdata2D->copyCoords(Shotgeom);
               Pdata2D->makeMap(poro_lmodel->getGeom(),SMAP);
               Pdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecP(Pdata2D);
            }
            if(Vxrecord){
               Vxdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Vxdata2D->setField(rockseis::VX);
               // Copy geometry to Data
               Vxdata2D->copyCoords(Shotgeom);
               Vxdata2D->makeMap(poro_lmodel->getGeom(),SMAP,0,0);
               Vxdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecVx(Vxdata2D);
            }
            if(Vzrecord){
               Vzdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Vzdata2D->setField(rockseis::VZ);
               // Copy geometry to Data
               Vzdata2D->copyCoords(Shotgeom);
               Vzdata2D->makeMap(poro_lmodel->getGeom(),SMAP,0,0);
               Vzdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecVz(Vzdata2D);
            }
            if(Szzrecord){
               Szzdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Szzdata2D->setField(rockseis::SZZ);
               // Copy geometry to Data
               Szzdata2D->copyCoords(Shotgeom);
               Szzdata2D->makeMap(poro_lmodel->getGeom(),SMAP);
               Szzdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecSzz(Szzdata2D);
            }
            if(Qxrecord){
               Qxdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Qxdata2D->setField(rockseis::QX);
               // Copy geometry to Data
               Qxdata2D->copyCoords(Shotgeom);
               Qxdata2D->makeMap(poro_lmodel->getGeom(),SMAP,0,0);
               Qxdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecQx(Qxdata2D);
            }
            if(Qzrecord){
               Qzdata2D = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
               Qzdata2D->setField(rockseis::QZ);
               // Copy geometry to Data
               Qzdata2D->copyCoords(Shotgeom);
               Qzdata2D->makeMap(poro_lmodel->getGeom(),SMAP,0,0);
               Qzdata2D->makeMap(poro_lmodel->getGeom(),GMAP,(poro_lmodel->getDomain())->getPadl(0),(poro_lmodel->getDomain())->getPadl(2),(poro_lmodel->getDomain())->getPadh(0),(poro_lmodel->getDomain())->getPadh(2));
               modelling->setRecQz(Qzdata2D);
            }


            // Run modelling 
            modelling->run();

            // Output record
            if(Szzrecord){
               Szzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Szzdata2Di->setFile(Szzrecordfile);
               interp->interp(Szzdata2D, Szzdata2Di);
               Sort->put2DGather(Szzdata2Di, work.id, (Szzdata2D->getGeom())->getGmap());
            }
            if(Vxrecord){
               Vxdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Vxdata2Di->setFile(Vxrecordfile);
               interp->interp(Vxdata2D, Vxdata2Di);
               Sort->put2DGather(Vxdata2Di, work.id, (Vxdata2D->getGeom())->getGmap());
            }
            if(Vzrecord){
               Vzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Vzdata2Di->setFile(Vzrecordfile);
               interp->interp(Vzdata2D, Vzdata2Di);
               Sort->put2DGather(Vzdata2Di, work.id, (Vzdata2D->getGeom())->getGmap());
            }

            if(Precord){
               Pdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Pdata2Di->setFile(Precordfile);
               interp->interp(Pdata2D, Pdata2Di);
               Sort->put2DGather(Pdata2Di, work.id, (Pdata2D->getGeom())->getGmap());
            }
            if(Qxrecord){
               Qxdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Qxdata2Di->setFile(Qxrecordfile);
               interp->interp(Qxdata2D, Qxdata2Di);
               Sort->put2DGather(Qxdata2Di, work.id, (Qxdata2D->getGeom())->getGmap());
            }
            if(Qzrecord){
               Qzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, ntrec, dtrec, 0.0);
               Qzdata2Di->setFile(Qzrecordfile);
               interp->interp(Qzdata2D, Qzdata2Di);
               Sort->put2DGather(Qzdata2Di, work.id, (Qzdata2D->getGeom())->getGmap());
            }

            // Reset all classes
            Shotgeom.reset();
            acu_lmodel.reset();
            poro_lmodel.reset();
            modelling.reset();
            if(Precord){
               Pdata2D.reset();
               Pdata2Di.reset();
            }
            if(Qxrecord){
               Qxdata2D.reset();
               Qxdata2Di.reset();
            }
            if(Qzrecord){
               Qzdata2D.reset();
               Qzdata2Di.reset();
            }
            if(Szzrecord){
               Szzdata2D.reset();
               Szzdata2Di.reset();
            }
            if(Vxrecord){
               Vxdata2D.reset();
               Vxdata2Di.reset();
            }
            if(Vzrecord){
               Vzdata2D.reset();
               Vzdata2Di.reset();
            }

            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }
   }
}

