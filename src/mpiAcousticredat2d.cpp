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
         PRINT_DOC(# MPI 2d acoustic modelling default configuration file);
         PRINT_DOC();
         PRINT_DOC(# Domain decomposition parameter);
         PRINT_DOC(        ndomain0 = "1";  # Number of domains along x direction to split the model into);
         PRINT_DOC(        ndomain1 = "1";  # Number of domains along z direction to split the model into);
         PRINT_DOC();
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(        freesurface = "true";  # True if free surface should be on);
         PRINT_DOC(            order = "8";  # Order of finite difference stencil 2-8);
         PRINT_DOC(            lpml = "18"; # Size of pml absorbing boundary (should be larger than order + 5 ));
         PRINT_DOC(            snapinc = "10"; # Snap interval in multiples of modelling interval);
         PRINT_DOC(            dtrec = "4e-3"; # Recording interval in seconds);
         PRINT_DOC(            apertx = "900"; # Aperture for local model (source is in the middle));
         PRINT_DOC();
         PRINT_DOC(# Booleans);
         PRINT_DOC(            Precord = "true";  # Set these to true if recording or snapshoting is to be made.);
         PRINT_DOC(            Vxrecord = "true";);
         PRINT_DOC(        Vzrecord = "true";);
         PRINT_DOC(        Psnap = "false";);
         PRINT_DOC(        Vxsnap = "false";);
         PRINT_DOC(        Vzsnap = "false";);
         PRINT_DOC();
         PRINT_DOC(# Files);
         PRINT_DOC(        Vp = "Vp2d.rss";);
         PRINT_DOC(        Rho = "Rho2d.rss";);
         PRINT_DOC(        Wavelet = "Wav2d.rss";);
         PRINT_DOC(        Input_Pdata = "2Ddata.rss";);
         PRINT_DOC(        Datum_survey = "2Dsurvey.rss";);
         PRINT_DOC(        Precordfile = "Pshot.rss";);
         PRINT_DOC(        Vxrecordfile = "Vxshot.rss";);
         PRINT_DOC(        Vzrecordfile = "Vzshot.rss";);
         PRINT_DOC(        Psnapfile = "Psnaps.rss";);
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
   float dtrec;
   std::string Indatafile;
   std::string Surveyfile;
   std::string Waveletfile;

   std::shared_ptr<rockseis::Data2D<float>> Indata2D;
   std::shared_ptr<rockseis::Data2D<float>> Indata2Di;
   std::string Vpfile;
   std::string Rhofile;
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

   // Create a local model class
   std::shared_ptr<rockseis::ModelAcoustic2D<float>> lmodel;

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
   if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Input_Pdata", &Indatafile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Datum_survey", &Surveyfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Psnap", &Psnap) == INPARSE_ERR) status = true;
   if(Psnap){
      if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vxsnap", &Vxsnap) == INPARSE_ERR) status = true;
   if(Vxsnap){
      if(Inpar->getPar("Vxsnapfile", &Vxsnapfile) == INPARSE_ERR) status = true;
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

   // Setup Domain decomposition
   mpi.setNdomain(ndomain0*ndomain1);
   mpi.splitDomains();

   // Create a sort classes
   std::shared_ptr<rockseis::Sort<float>> Sort_in (new rockseis::Sort<float>());
   Sort_in->setDatafile(Indatafile);
   std::shared_ptr<rockseis::Sort<float>> Sort_out (new rockseis::Sort<float>());
   Sort_out->setDatafile(Surveyfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelAcoustic2D<float>> gmodel (new rockseis::ModelAcoustic2D<float>(Vpfile, Rhofile, lpml ,fs));

   // Create a data class for the source wavelet
   std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(Waveletfile));

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
      std::shared_ptr<rockseis::Data2D<float>> Shotgeom;
      std::shared_ptr<rockseis::ModellingAcoustic2D<float>> modelling;
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
            Indata2D = Sort_in->get2DGather(work.id);
            size_t ntr_in = Indata2D->getNtrace();
            lmodel = gmodel->getDomainmodel(Indata2D, apertx, SMAP, mpi.getDomainrank(), ndomain0, ndomain1, order);
            (lmodel->getDomain())->setMpi(&mpi);

            // Interpolate input data to avoid instability       
            Indata2Di = std::make_shared<rockseis::Data2D<float>>(ntr_in, source->getNt(), source->getDt(), 0.0);
            interp->interp(Indata2D, Indata2Di);
            Indata2Di->makeMap(lmodel->getGeom(), GMAP);
            Indata2Di->copyGmap2Smap();

            // Get the Datum geometry
            Shotgeom = Sort_out->get2DGather(work.id);
            size_t ntr_out = Shotgeom->getNtrace();

            modelling = std::make_shared<rockseis::ModellingAcoustic2D<float>>(lmodel, Indata2Di, order, snapinc);
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

            // Setting Record
            if(Precord){
               Pdata2D = std::make_shared<rockseis::Data2D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
               Pdata2D->setField(rockseis::PRESSURE);
               // Copy geometry to Data
               Pdata2D->copyCoords(Shotgeom);
               Pdata2D->makeMap(lmodel->getGeom(),SMAP);
               Pdata2D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(2));
               modelling->setRecP(Pdata2D);
            }
            if(Vxrecord){
               Vxdata2D = std::make_shared<rockseis::Data2D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
               Vxdata2D->setField(rockseis::VX);
               // Copy geometry to Data
               Vxdata2D->copyCoords(Shotgeom);
               Vxdata2D->makeMap(lmodel->getGeom(),SMAP,0,0);
               Vxdata2D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(2));
               modelling->setRecVx(Vxdata2D);
            }
            if(Vzrecord){
               Vzdata2D = std::make_shared<rockseis::Data2D<float>>(ntr_out, source->getNt(), source->getDt(), 0.0);
               Vzdata2D->setField(rockseis::VZ);
               // Copy geometry to Data
               Vzdata2D->copyCoords(Shotgeom);
               Vzdata2D->makeMap(lmodel->getGeom(),SMAP,0,0);
               Vzdata2D->makeMap(lmodel->getGeom(),GMAP,(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadh(2));
               modelling->setRecVz(Vzdata2D);
            }

            // Run modelling 
            modelling->run();

            // Output record
            if(Precord){
               Pdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr_out, ntrec, dtrec, 0.0);
               Pdata2Di->setFile(Precordfile);
               interp->interp(Pdata2D, Pdata2Di);
               Sort_out->put2DGather(Pdata2Di, work.id, (Pdata2D->getGeom())->getGmap());
            }
            if(Vxrecord){
               Vxdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr_out, ntrec, dtrec, 0.0);
               Vxdata2Di->setFile(Vxrecordfile);
               interp->interp(Vxdata2D, Vxdata2Di);
               Sort_out->put2DGather(Vxdata2Di, work.id, (Vxdata2D->getGeom())->getGmap());
            }
            if(Vzrecord){
               Vzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr_out, ntrec, dtrec, 0.0);
               Vzdata2Di->setFile(Vzrecordfile);
               interp->interp(Vzdata2D, Vzdata2Di);
               Sort_out->put2DGather(Vzdata2Di, work.id, (Vzdata2D->getGeom())->getGmap());
            }

            // Reset all classes
            Pdata2D.reset();
            Pdata2Di.reset();
            Shotgeom.reset();
            lmodel.reset();
            modelling.reset();
            if(Precord){
               Pdata2D.reset();
               Pdata2Di.reset();
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

