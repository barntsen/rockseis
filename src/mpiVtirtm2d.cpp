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
   MPIdomaindecomp mpi = MPIdomaindecomp(&argc,&argv);
   if(mpi.getNrank() < 2){
      rs_error("This is a parallel program, it must run with at least 2 processors, use mpirun.");
   }

   if(argc < 2){
      if(mpi.getRank() == 0){
         PRINT_DOC(# MPI 2d vti reverse-time migration configuration file);
         PRINT_DOC(# Domain decomposition parameter);
         PRINT_DOC(        ndomain0 = "1";  # Number of domains along x direction to split the model into);
         PRINT_DOC(        ndomain1 = "1";  # Number of domains along z direction to split the model into);
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(freesurface = "false";  # True if free surface should be on);
         PRINT_DOC(order = "8";  # Order of finite difference stencil);
         PRINT_DOC(lpml = "10"; # Size of pml absorbing boundary (should be larger than order + 5 ));
         PRINT_DOC(source_type = "0"; # Source type 0 - pressure. 1 for Vx. 3 for Vz.);
         PRINT_DOC(snapinc = "4"; # Snap interval in multiples of modelling interval);
         PRINT_DOC(apertx = "0"; # Aperture for local model (source is in the middle));
         PRINT_DOC();
         PRINT_DOC(# Checkpointing parameters);
         PRINT_DOC(snapmethod = "0";  );
         PRINT_DOC(nsnaps = "11";);
         PRINT_DOC(incore = "true";  );
         PRINT_DOC();
         PRINT_DOC(# Migration parameters);
         PRINT_DOC(nhx = "1";);
         PRINT_DOC(nhz = "1";);
         PRINT_DOC();
         PRINT_DOC(# Booleans);
         PRINT_DOC(Pimaging = "true";  # Set these to true if imaging of these events is to be made.);
         PRINT_DOC(Pgather = "false"; # If surface gathers are to be output);
         PRINT_DOC(Simaging = "true";); 
         PRINT_DOC(Sgather = "false"; # If surface gathers are to be output);
         PRINT_DOC();
         PRINT_DOC(# Files);
         PRINT_DOC(        C11 = "c11.rss";);
         PRINT_DOC(        C13 = "c13.rss";);
         PRINT_DOC(        C33 = "c33.rss";);
         PRINT_DOC(        C55 = "c55.rss";);
         PRINT_DOC(Rho = "Rho2d.rss";);
         PRINT_DOC(Wavelet = "Wav2d.rss";);
         PRINT_DOC(Precordfile = "Pshot.rss";);
         PRINT_DOC(Vxrecordfile = "Vxshot.rss";);
         PRINT_DOC(Vzrecordfile = "Vzshot.rss";);
         PRINT_DOC();
         PRINT_DOC(# Output Files);
         PRINT_DOC(Pimagefile = "Pimage2d.rss";);
         PRINT_DOC(Pgatherfile = "Pgather2d.rss";);
         PRINT_DOC(Simagefile = "Simage2d.rss";);
         PRINT_DOC(Sgatherfile = "Sgather2d.rss";);
         PRINT_DOC(Snapfile = "snaps.rss";);
         PRINT_DOC();
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   int lpml;
   bool fs;
   int ndomain0;
   int ndomain1;
   bool incore = false;
   int order;
   int snapinc;
   int nsnaps = 0;
   int snapmethod;
   float apertx;
   int stype;
   int nhx=1, nhz=1;
   std::string Waveletfile;
   std::string C11file;
   std::string C13file;
   std::string C33file;
   std::string C55file;
   std::string Rhofile;
   std::string Snapfile;

   bool Pimaging;
   std::string Pimagefile;
   std::shared_ptr<rockseis::Image2D<float>> pimage;
   bool Pgather;
   std::string Pgatherfile;
   std::shared_ptr<rockseis::Data2D<float>> pgather;

   bool Simaging;
   std::string Simagefile;
   std::shared_ptr<rockseis::Image2D<float>> simage;
   bool Sgather;
   std::string Sgatherfile;
   std::shared_ptr<rockseis::Data2D<float>> sgather;

   std::string Precordfile;
   std::shared_ptr<rockseis::Data2D<float>> Pdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Pdata2Di;

   std::string Vxrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Vxdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Vxdata2Di;

   std::string Vzrecordfile;
   std::shared_ptr<rockseis::Data2D<float>> Vzdata2D;
   std::shared_ptr<rockseis::Data2D<float>> Vzdata2Di;

   // Create a local model class
   std::shared_ptr<rockseis::ModelVti2D<float>> lmodel;

   /* Get parameters from configuration file */
   std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
   if(Inpar->parse(argv[1]) == INPARSE_ERR) 
   {
      rs_error("Parse error on input config file", argv[1]);
   }
   status = false; 
   if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
   if(Inpar->getPar("order", &order) == INPARSE_ERR) status = true;
   if(Inpar->getPar("ndomain0", &ndomain0) == INPARSE_ERR) status = true;
   if(Inpar->getPar("ndomain1", &ndomain1) == INPARSE_ERR) status = true;
   if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("source_type", &stype) == INPARSE_ERR) status = true;
   if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C11", &C11file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C13", &C13file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C33", &C33file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("C55", &C55file) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Pimaging", &Pimaging) == INPARSE_ERR) status = true;
   if(Pimaging){
      if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Simaging", &Simaging) == INPARSE_ERR) status = true;
   if(Simaging){
      if(Inpar->getPar("Simagefile", &Simagefile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Pgather", &Pgather) == INPARSE_ERR) status = true;
   if(Pgather){
      if(Inpar->getPar("Pgatherfile", &Pgatherfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Sgather", &Sgather) == INPARSE_ERR) status = true;
   if(Sgather){
      if(Inpar->getPar("Sgatherfile", &Sgatherfile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Snapfile", &Snapfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Vxrecordfile", &Vxrecordfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Vzrecordfile", &Vzrecordfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
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

   // Setup Domain decomposition
   mpi.setNdomain(ndomain0*ndomain1);
   mpi.splitDomains();

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Vxrecordfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelVti2D<float>> gmodel (new rockseis::ModelVti2D<float>(C11file, C13file, C33file, C55file, Rhofile, lpml, fs));

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

      // Create work queue
      for(long int i=0; i<ngathers; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
         mpi.addWork(work);
      }

      // Perform work in parallel
      mpi.performWork();


      // Image gathers
      if(Pgather){
         std::shared_ptr<rockseis::File> Fimg (new rockseis::File());
         Fimg->input(Pimagefile + "-" + std::to_string(0));
         pgather = std::make_shared<rockseis::Data2D<float>>(Fimg->getN(1),Fimg->getN(3),Fimg->getD(3),Fimg->getO(3));
         pgather->setFile(Pgatherfile);
         pgather->open("o");
         for(long int i=0; i<ngathers; i++) {
            pgather->putImage(Pimagefile + "-" + std::to_string(i));
         }
         pgather->close();
         Fimg->close();
      }
      if(Sgather){
         std::shared_ptr<rockseis::File> Fimg (new rockseis::File());
         Fimg->input(Simagefile + "-" + std::to_string(0));
         sgather = std::make_shared<rockseis::Data2D<float>>(Fimg->getN(1),Fimg->getN(3),Fimg->getD(3),Fimg->getO(3));
         sgather->setFile(Sgatherfile);
         sgather->open("o");
         for(long int i=0; i<ngathers; i++) {
            sgather->putImage(Simagefile + "-" + std::to_string(i));
         }
         sgather->close();
         Fimg->close();
      }

      // Images
      if(Pimaging){
         pimage = std::make_shared<rockseis::Image2D<float>>(Pimagefile, gmodel, nhx, nhz);
         pimage->createEmpty();
         for(long int i=0; i<ngathers; i++) {
            pimage->stackImage(Pimagefile + "-" + std::to_string(i));
            remove_file(Pimagefile + "-" + std::to_string(i));
         }
      }

      if(Simaging){
         simage = std::make_shared<rockseis::Image2D<float>>(Simagefile, gmodel, nhx, nhz);
         simage->createEmpty();
         for(long int i=0; i<ngathers; i++) {
            simage->stackImage(Simagefile + "-" + std::to_string(i));
            remove_file(Simagefile + "-" + std::to_string(i));
         }
      }
   }
   else {
      /* Slave */
      std::shared_ptr<rockseis::RtmVti2D<float>> rtm;
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
            // Do migration
            Sort->readKeymap();
            Sort->readSortmap();

            // Get the shot
            Sort->setDatafile(Vxrecordfile);
            Vxdata2D = Sort->get2DGather(work.id);
            size_t ntr = Vxdata2D->getNtrace();

            Sort->setDatafile(Vzrecordfile);
            Vzdata2D = Sort->get2DGather(work.id);

            Sort->setDatafile(Precordfile);
            Pdata2D = Sort->get2DGather(work.id);

            if(mpi.getDomainrank() == 0){
               // create and write empty gradfiles, with the size of the local model
               lmodel = gmodel->getLocal(Vxdata2D, apertx, SMAP);
               if(Pimaging){
                  pimage = std::make_shared<rockseis::Image2D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, 1, 1);
                  pimage->createEmpty();
                  pimage.reset();
               }
               if(Simaging){
                  simage = std::make_shared<rockseis::Image2D<float>>(Simagefile + "-" + std::to_string(work.id), lmodel, 1, 1);
                  simage->createEmpty();
                  simage.reset();
               }
               lmodel.reset();
               // Synchronize
               mpi.barrier();
            }else{
               // Synchronize
               mpi.barrier();
            }

            lmodel = gmodel->getDomainmodel(Vxdata2D, apertx, SMAP, mpi.getDomainrank(), ndomain0, ndomain1, order);
            (lmodel->getDomain())->setMpi(&mpi);

            // Read wavelet data, set shot coordinates and make a map
            source->read();
            source->copyCoords(Vxdata2D);

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
            source->makeMap(lmodel->getGeom(), SMAP);

            // Interpolate shot
            Pdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            interp->interp(Pdata2D, Pdata2Di);
            Pdata2Di->setField(rockseis::PRESSURE);
            Pdata2Di->makeMap(lmodel->getGeom(), GMAP);

            Vxdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            interp->interp(Vxdata2D, Vxdata2Di);
            Vxdata2Di->setField(rockseis::VX);
            Vxdata2Di->makeMap(lmodel->getGeom(), GMAP);

            Vzdata2Di = std::make_shared<rockseis::Data2D<float>>(ntr, source->getNt(), source->getDt(), 0.0);
            interp->interp(Vzdata2D, Vzdata2Di);
            Vzdata2Di->setField(rockseis::VZ);
            Vzdata2Di->makeMap(lmodel->getGeom(), GMAP);

            rtm = std::make_shared<rockseis::RtmVti2D<float>>(lmodel, source, Pdata2Di, Vxdata2Di, Vzdata2Di, order, snapinc);

            // Setting Image objects
            if(Pimaging){ 
               pimage = std::make_shared<rockseis::Image2D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhz);
               rtm->setPimage(pimage);
            }
            if(Simaging){ 
               simage = std::make_shared<rockseis::Image2D<float>>(Simagefile + "-" + std::to_string(work.id), lmodel, nhx, nhz);
               rtm->setSimage(simage);
            }

            // Setting Snapshot file 
            rtm->setSnapfile(Snapfile + "-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));

            // Setting Snapshot parameters
            rtm->setNcheck(nsnaps);
            rtm->setIncore(incore);

            // Set logfile
            rtm->setLogfile("log.txt-" + std::to_string(work.id)+ "-" + std::to_string(mpi.getDomainrank()));


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
               pimage->stackImage_parallel(Pimagefile + "-" + std::to_string(work.id),(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(2));
            }
            if(Simaging){
               simage->stackImage_parallel(Simagefile + "-" + std::to_string(work.id),(lmodel->getDomain())->getPadl(0),(lmodel->getDomain())->getPadh(0),(lmodel->getDomain())->getPadl(2),(lmodel->getDomain())->getPadh(2));
            }

            // Synchronize
            mpi.barrier();

            // Reset all classes
            Pdata2D.reset();
            Pdata2D.reset();
            Vxdata2D.reset();
            Vxdata2Di.reset();
            Vzdata2D.reset();
            Vzdata2Di.reset();
            lmodel.reset();
            pimage.reset();
            simage.reset();
            rtm.reset();

            // Send result back
            work.status = WORK_FINISHED;
            mpi.sendResult(work);		
         }
      }
   }
}

