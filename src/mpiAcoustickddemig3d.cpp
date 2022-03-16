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
         PRINT_DOC(# MPI 2d acoustic reverse-time migration configuration file);
         PRINT_DOC();
         PRINT_DOC(# Traveltime computation parameters);
         PRINT_DOC();
         PRINT_DOC(incore = "false"; # Whether to compute the traveltimes on the fly. If false then atraveltime table must be built beforehand);
         PRINT_DOC();
         PRINT_DOC(# Modelling parameters);
         PRINT_DOC(apertx = "1800"; # Aperture for local model in x direction (source is in the middle));
         PRINT_DOC(aperty = "1800"; # Aperture for local model in y direction (source is in the middle));
         PRINT_DOC();
         PRINT_DOC(# Migration parameters);
         PRINT_DOC(freqinc = "4"; # Integer frequency interval to sum over);
         PRINT_DOC(minfreq = "100.0"; # Minimum frequency to migrate);
         PRINT_DOC(maxfreq = "100.0"; # Maximum frequency to migrate);
         PRINT_DOC(radius = "50.0"; # Radius of traveltime interpolation);
         PRINT_DOC();
         PRINT_DOC(# Files);
         PRINT_DOC(Vp = "Vp2d.rss";);
         PRINT_DOC(Wavelet = "Wav2d.rss";);
         PRINT_DOC(Ttable = "Ttable2d.rss";);
         PRINT_DOC(Survey = "3DSurvey.rss";);
         PRINT_DOC(Precordfile = "Pshots2d.rss";);
         PRINT_DOC(Pimagefile = "Pimage2d.rss";);
         PRINT_DOC();
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   int lpml = 3;
   float apertx;
   float aperty;
   int freqinc;
   float minfreq;
   float maxfreq;
   float radius;
   bool incore;
   std::string Vpfile;
   std::string Ttablefile;
   std::string Pimagefile;
   std::string Surveyfile;
   std::string Waveletfile;
   std::string Precordfile;
   std::shared_ptr<rockseis::Data3D<float>> shot3D;
   std::shared_ptr<rockseis::Image3D<float>> pimage;
   std::shared_ptr<rockseis::Image3D<float>> lpimage;
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
   if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
   if(!incore){
      if(Inpar->getPar("Ttable", &Ttablefile) == INPARSE_ERR) status = true;
   }
   if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
   if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Pimagefile", &Pimagefile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;

   if(status == true){
      rs_error("Program terminated due to input errors.");
   }

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Surveyfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelEikonal3D<float>> gmodel (new rockseis::ModelEikonal3D<float>(Vpfile, lpml));

   // Create a data class for the source wavelet
   std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Waveletfile));

   // Create an interpolation class
   std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

   // Compute record length in samples
   size_t ntrec; 
   float dtrec;
   ntrec = source->getNt();
   dtrec = source->getDt();


   // Test for problematic model sampling
   if(incore){
      if(gmodel->getDx() != gmodel->getDz()){
         rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dz).");
      }
   }

   if(mpi.getRank() == 0) {
      // Master
      Sort->createShotmap(Surveyfile); 
      Sort->writeKeymap();
      Sort->writeSortmap();

      // Get number of shots
      size_t ngathers =  Sort->getNensemb();

      // Create an empty data file
      Sort->createEmptydataset(Precordfile, ntrec, dtrec, 0.0);

      // Create work queue
      for(long int i=0; i<ngathers; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED,0});
         mpi.addWork(work);
      }

      // Perform work in parallel
      mpi.performWork();
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
            // Do demigration
            Sort->readKeymap();
            Sort->readSortmap();
            Sort->setDatafile(Precordfile);

            // Get the shot
            shot3D = Sort->get3DGather(work.id);

            // Make local model
            lmodel = gmodel->getLocal(shot3D, apertx, aperty, SMAP);
            lmodel->Expand();

            // Map coordinates to model
            shot3D->makeMap(lmodel->getGeom(), SMAP);
            shot3D->makeMap(lmodel->getGeom(), GMAP);

            // Make image class
            pimage = std::make_shared<rockseis::Image3D<float>>(Pimagefile);
            lpimage = pimage->getLocal(shot3D, apertx, aperty, SMAP);

            if(!incore){
                // Create traveltime table class
                ttable = std::make_shared<rockseis::Ttable3D<float>>(Ttablefile);
                ttable->allocTtable();
            }

            // Create imaging class
            kdmig = std::make_shared<rockseis::KdmigAcoustic3D<float>>(lmodel, ttable, shot3D, lpimage);

            // Set traveltime parameters
            kdmig->setIncore(incore);

            // Set frequency decimation 
            kdmig->setFreqinc(1);

            // Set minimum and maximum frequency to migrate
            kdmig->setMinfreq(0.0);
            kdmig->setMaxfreq(125.0);

            // Set radius of interpolation
            kdmig->setRadius(radius);


            // Set logfile
            //kdmig->setLogfile("log.txt-" + std::to_string(work.id));

            // Run demigration
            kdmig->run_demig();

            shot3D->setFile(Precordfile);
            //Sort->put3DGather(shot3D, work.id,(shot3D->getGeom())->getGmap());
            Sort->put3DGather(shot3D, work.id);

            // Reset all classes
            shot3D.reset();
            lmodel.reset();
            pimage.reset();
            lpimage.reset();
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

