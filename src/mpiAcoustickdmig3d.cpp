#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "kdmig.h"
#include "ttable.h"
#include "rays.h"
#include "pml.h"
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
         PRINT_DOC(# MPI 3d acoustic reverse-time migration configuration file);
         PRINT_DOC();
         PRINT_DOC(# Traveltime computation parameters);
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(apertx = "1800"; # Aperture for local model (source is in the middle));
         PRINT_DOC(aperty = "1800"; # Aperture for local model (source is in the middle));
         PRINT_DOC();
         PRINT_DOC(# Migration parameters);
         PRINT_DOC(freqinc = "4"; # Integer frequency interval to sum over);
         PRINT_DOC(minfreq = "100.0"; # Minimum frequency to migrate);
         PRINT_DOC(maxfreq = "100.0"; # Maximum frequency to migrate);
         PRINT_DOC(radius = "50.0"; # Radius of traveltime interpolation);
         PRINT_DOC(nhx = "1";);
         PRINT_DOC(nhy = "1";);
         PRINT_DOC(nhz = "1";);
         PRINT_DOC();
         PRINT_DOC(# Booleans);
         PRINT_DOC(Gather = "false"; # If surface gathers are to be output);
         PRINT_DOC();
         PRINT_DOC(# Input files);
         PRINT_DOC(Vp = "Vp3d.rss";);
         PRINT_DOC(Ttable = "Ttable3d.rss";);
         PRINT_DOC(Precordfile = "Pshots3d.rss";);
         PRINT_DOC();
         PRINT_DOC(# Output files);
         PRINT_DOC(Pimagefile = "Pimage3d.rss";);
         PRINT_DOC(Pgatherfile = "Pgather3d.rss";);
         PRINT_DOC();
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   int lpml = 3;
   float apertx;
   float aperty;
   int nhx=1, nhy=1, nhz=1;
   int freqinc;
   float minfreq;
   float maxfreq;
   float radius;
   bool incore;
   bool homogen;
   std::string Vpfile;
   std::string Ttablefile;
   std::string Pimagefile;
   std::string Precordfile;
   std::shared_ptr<rockseis::Data3D<float>> shot3D;
   std::shared_ptr<rockseis::Image3D<float>> pimage;
   bool Gather;
   std::string Pgatherfile;
   std::shared_ptr<rockseis::Data3D<float>> pgather;
   std::shared_ptr<rockseis::ModelEikonal3D<float>> lmodel;

   std::shared_ptr<rockseis::Ttable3D<float>> ttable;

   /* Get parameters from configuration file */
   std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
   if(Inpar->parse(argv[1]) == INPARSE_ERR) 
   {
      rs_error("Parse error on input config file", argv[1]);
   }
   status = false; 
   if(Inpar->getPar("freqinc", &freqinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("minfreq", &minfreq) == INPARSE_ERR) status = true;
   if(Inpar->getPar("maxfreq", &maxfreq) == INPARSE_ERR) status = true;
   if(Inpar->getPar("radius", &radius) == INPARSE_ERR) status = true;
   if(Inpar->getPar("incore", &incore) == INPARSE_ERR) status = true;
   if(incore){
       if(Inpar->getPar("homogen", &homogen) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
   if(!incore){
       if(Inpar->getPar("Ttable", &Ttablefile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("nhy", &nhy) == INPARSE_ERR) status = true;
   if(Inpar->getPar("nhz", &nhz) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Gather", &Gather) == INPARSE_ERR) status = true;
   if(Gather){
      if(Inpar->getPar("Pgatherfile", &Pgatherfile) == INPARSE_ERR) status = true;
   }

   if(status == true){
      rs_error("Program terminated due to input errors.");
   }

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Precordfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelEikonal3D<float>> gmodel (new rockseis::ModelEikonal3D<float>(Vpfile, lpml));

   // Test for problematic model sampling
   if(incore && !homogen){
       if(gmodel->getDx() != gmodel->getDz() || gmodel->getDx() != gmodel->getDy()){
           rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dy = dz).");
       }
   }
   mpi.setVerbose(false);


   if(mpi.getRank() == 0) {
      // Master
      Sort->createShotmap(Precordfile); 
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
      if(Gather){
         std::shared_ptr<rockseis::File> Fimg (new rockseis::File());
         Fimg->input(Pimagefile + "-" + std::to_string(0));
         pgather = std::make_shared<rockseis::Data3D<float>>(Fimg->getN(1)*Fimg->getN(2),Fimg->getN(3),Fimg->getD(3),Fimg->getO(3));
         pgather->setFile(Pgatherfile);
         pgather->open("o");
         for(long int i=0; i<ngathers; i++) {
            pgather->putImage(Pimagefile + "-" + std::to_string(i));
         }
         pgather->close();
         Fimg->close();
      }

      // Image
      pimage = std::make_shared<rockseis::Image3D<float>>(Pimagefile, gmodel, nhx, nhy, nhz);
      pimage->createEmpty();
      for(long int i=0; i<ngathers; i++) {
         pimage->stackImage(Pimagefile + "-" + std::to_string(i));
         remove_file(Pimagefile + "-" + std::to_string(i));
      }

   }
   else {
      /* Slave */
      std::shared_ptr<rockseis::KdmigAcoustic3D<float>> kdmig;
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
            shot3D = Sort->get3DGather(work.id);

            // Make local model
            lmodel = gmodel->getLocal(shot3D, apertx, aperty, SMAP);
            lmodel->Expand();

            // Map coordinates to model
            shot3D->makeMap(lmodel->getGeom(), SMAP);
            shot3D->makeMap(lmodel->getGeom(), GMAP);

            // Make image class
            pimage = std::make_shared<rockseis::Image3D<float>>(Pimagefile + "-" + std::to_string(work.id), lmodel, nhx, nhy, nhz);
            if(!incore){
                // Create traveltime table class
                ttable = std::make_shared<rockseis::Ttable3D<float>>(Ttablefile);
                ttable->allocTtable();
            }

            // Create imaging class
            kdmig = std::make_shared<rockseis::KdmigAcoustic3D<float>>(lmodel, ttable, shot3D, pimage);

            // Set traveltime parameters
            kdmig->setIncore(incore);
            kdmig->setHomogen(homogen);

            // Set frequency decimation 
            kdmig->setFreqinc(freqinc);

            // Set minimum and maximum frequency to migrate
            kdmig->setMinfreq(minfreq);
            kdmig->setMaxfreq(maxfreq);

            // Set radius of interpolation
            kdmig->setRadius(radius);

            // Set logfile
            //kdmig->setLogfile("log.txt-" + std::to_string(work.id));

            // Run migration
            kdmig->run();

            // Output image
            pimage->write();

            // Reset all classes
            shot3D.reset();
            lmodel.reset();
            pimage.reset();
            kdmig.reset();
            if(!incore){
                ttable.reset();
            }

            // Send result back
            work.status = WORK_FINISHED;
            mpi.sendResult(work);		
         }
      }
   }
}

