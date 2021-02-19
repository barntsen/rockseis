#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "rays.h"
#include "ttable.h"
#include "utils.h"
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
         PRINT_DOC(# MPI 3d acoustic modelling default configuration file);
         PRINT_DOC();
         PRINT_DOC(# Sampling);
         PRINT_DOC(        Souinc = "1";);
         PRINT_DOC(        Recinc = "1";);
         PRINT_DOC();
         PRINT_DOC(# Input files);
         PRINT_DOC(        Vp = "Vp3d.rss";);
         PRINT_DOC(        Survey = "3DSurvey.rss";);
         PRINT_DOC();
         PRINT_DOC(# Output files);
         PRINT_DOC(        Ttablefile = "Ttable3d.rss";);
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   int souinc, recinc;
   int nsoufin, nrecfin;
   int Lpml = 0;
   std::string Surveyfile;
   std::string Vpfile;
   std::string Ttablefile;
   std::shared_ptr<rockseis::Data3D<float>> source;
   std::shared_ptr<rockseis::RaysAcoustic3D<float>> rays;
   std::shared_ptr<rockseis::Ttable3D<float>> Ttable;


   /* Get parameters from configuration file */
   std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
   if(Inpar->parse(argv[1]) == INPARSE_ERR) 
   {
      rs_error("Parse error on input config file", argv[1]);
   }
   status = false; 
   if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Survey", &Surveyfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Souinc", &souinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Recinc", &recinc) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Ttablefile", &Ttablefile) == INPARSE_ERR) status = true;

   if(status == true){
      rs_error("Program terminated due to input errors.");
   }

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Surveyfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelEikonal3D<float>> gmodel (new rockseis::ModelEikonal3D<float>(Vpfile, Lpml));

   // Test for problematic model sampling
   if(gmodel->getDx() != gmodel->getDz() || gmodel->getDx() != gmodel->getDy()){
      rs_error("Input model has different dx and dz values. This is currently not allowed. Interpolate to a unique grid sampling value (i.e dx = dy = dz).");
   }
   // Read and expand global model
   gmodel->readVelocity();
   gmodel->Expand();

   // Create an interpolation class
   std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

   if(mpi.getRank() == 0) {
      // Master

      // Get number of receivers
      Sort->createReceivermap(Surveyfile); 
      size_t nrecgath =  Sort->getNensemb();

      // Get number of shots
      Sort->createShotmap(Surveyfile); 
      size_t nsougath =  Sort->getNensemb();

      Sort->writeKeymap();
      Sort->writeSortmap();

      nsoufin = nsougath/souinc + 1;
      if(nsoufin > nsougath) nsoufin = nsougath;

      nrecfin = nrecgath/recinc + 1;
      if(nrecfin > nrecgath) nrecfin = nrecgath;

      // Create a travel time table class
      size_t ngathers = nsoufin + nrecfin;
      Ttable = std::make_shared<rockseis::Ttable3D<float>> (gmodel, ngathers);
      Ttable->setFilename(Ttablefile);
      Ttable->createEmptyttable();

      /******************             Creating source side     ***********************/
      // Create work queue
      for(long int i=0; i<nsoufin; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
         mpi.addWork(work);
      }

      // Broadcast nsoufin
      MPI_Bcast(&nsoufin, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // Perform work in parallel
      mpi.performWork();

      // Reset mpi 
      mpi.clearWork();
      /******************             Creating receiver side     ***********************/

      // Create new list of positions
      Sort->createReceivermap(Surveyfile); 
      Sort->setReciprocity(true);
      Sort->writeKeymap();
      Sort->writeSortmap();

      for(long int i=0; i<nrecfin; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
         mpi.addWork(work);
      }

      // Perform work in parallel
      mpi.performWork();

   }
   else {
      /* Slave */
      size_t number;
      std::shared_ptr<rockseis::Data3D<float>> Shotgeom;

      // Receive nsoufin from master
      MPI_Bcast(&nsoufin, 1, MPI_INT, 0, MPI_COMM_WORLD);

      /******************             Creating source side     ***********************/
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

            number = work.id*souinc;
            if (number > Sort->getNensemb()-1) number = Sort->getNensemb()-1;
            Shotgeom = Sort->get3DGather(number);

            // Set shot coordinates and make a map
            source = std::make_shared<rockseis::Data3D<float>>(1, 1, 1.0, 0.0);
            source->copyCoords(Shotgeom);

            source->makeMap(gmodel->getGeom(), SMAP);

            // Run modelling 
            rays = std::make_shared<rockseis::RaysAcoustic3D<float>>(gmodel);

            /* initialize traveltime field at source positions */
            rays->insertSource(source, SMAP);
            rays->solve();

            // Create traveltime table
            Ttable = std::make_shared<rockseis::Ttable3D<float>> (Ttablefile);
            Ttable->fetchTtabledata(rays, source, work.id); //Get traveltime data
            Ttable->writeTtable(work.id);
            Ttable.reset();

            // Reset all classes
            source.reset();
            Shotgeom.reset();
            rays.reset();
            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }
      /******************             Creating receiver side     ***********************/
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

            number = work.id*recinc;
            if (number > Sort->getNensemb()-1) number = Sort->getNensemb()-1;
            Shotgeom = Sort->get3DGather(number);

            // Set shot coordinates and make a map
            source = std::make_shared<rockseis::Data3D<float>>(1, 1, 1.0, 0.0);
            source->copyCoords(Shotgeom);

            source->makeMap(gmodel->getGeom(), SMAP);

            // Run modelling 
            rays = std::make_shared<rockseis::RaysAcoustic3D<float>>(gmodel);

            /* initialize traveltime field at source positions */
            rays->insertSource(source, SMAP);
            rays->solve();

            // Create traveltime table
            Ttable = std::make_shared<rockseis::Ttable3D<float>> (Ttablefile);
            Ttable->fetchTtabledata(rays, source, work.id+nsoufin); //Get traveltime data
            Ttable->writeTtable(work.id+nsoufin);
            Ttable.reset();

            // Reset all classes
            source.reset();
            Shotgeom.reset();
            rays.reset();
            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }
   }
}

