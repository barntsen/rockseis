#include <iostream>
#include <fstream>
#include <mpi.h>
#include <memory>
#include <math.h>
#include <inparse.h>
#include "model.h"
#include "image.h"
#include "rays.h"
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
         PRINT_DOC(# MPI 2d acoustic modelling default configuration file);
         PRINT_DOC();
         PRINT_DOC(# Input files);
         PRINT_DOC(        Vp = "Vp2d.rss";);
         PRINT_DOC(        Trecordfile = "Trecord.rss";);
         PRINT_DOC();
         PRINT_DOC(# Output files);
         PRINT_DOC(       Raydensfile = "Rayden.rss";);
      }
      exit(1);
   }
   bool status;
   /* General input parameters */
   std::string Trecordfile;
   std::string Raydensfile;
   std::string Vpfile;

   std::shared_ptr<rockseis::Image2D<float>> rayden;

   // Create a local model class
   std::shared_ptr<rockseis::ModelEikonal2D<float>> lmodel;

   /* Get parameters from configuration file */
   std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
   if(Inpar->parse(argv[1]) == INPARSE_ERR) 
   {
      rs_error("Parse error on input config file", argv[1]);
   }
   status = false; 
   if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Trecordfile", &Trecordfile) == INPARSE_ERR) status = true;
   if(Inpar->getPar("Raydensfile", &Raydensfile) == INPARSE_ERR) status = true;

   if(status == true){
      rs_error("Program terminated due to input errors.");
   }

   // Create a sort class
   std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
   Sort->setDatafile(Trecordfile);

   // Create a global model class
   std::shared_ptr<rockseis::ModelEikonal2D<float>> gmodel (new rockseis::ModelEikonal2D<float>(Vpfile, 10));

   // Create an interpolation class
   std::shared_ptr<rockseis::Interp<float>> interp (new rockseis::Interp<float>(SINC));

   if(mpi.getRank() == 0) {
      // Master
      Sort->createShotmap(Trecordfile); 
      Sort->writeKeymap();
      Sort->writeSortmap();

      // Get number of shots
      size_t ngathers =  Sort->getNensemb();

      // Create work queue
      for(long int i=0; i<ngathers; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
         mpi.addWork(work);
      }

      // Perform work in parallel
      mpi.performWork();

      // Stack gradient
      rayden = std::make_shared<rockseis::Image2D<float>>(Raydensfile, gmodel, 1, 1);
      rayden->createEmpty();

      for(long int i=0; i<ngathers; i++) {
         rayden->stackImage(Raydensfile + "-" + std::to_string(i));
         remove_file(Raydensfile + "-" + std::to_string(i));
      }


   }
   else {
      /* Slave */
      std::shared_ptr<rockseis::Data2D<float>> Shotgeom;
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

            Shotgeom = Sort->get2DGather(work.id);
            size_t ntr = Shotgeom->getNtrace();
            lmodel = gmodel->getLocal(Shotgeom, -1.0*gmodel->getDx(), SMAP);
            lmodel->Expand();

            // Set shot coordinates and make a map
            std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(1, 1, 1.0, 0.0));
            source->copyCoords(Shotgeom);
            source->makeMap(lmodel->getGeom(), SMAP);

            // Run modelling 
            std::shared_ptr<rockseis::RaysAcoustic2D<float>> rays (new rockseis::RaysAcoustic2D<float>(lmodel));

            /* initialize traveltime field at source positions */
            rays->insertSource(source, SMAP);
            rays->solve();

            std::shared_ptr<rockseis::Data2D<float>> Tres2D;
            Tres2D = std::make_shared<rockseis::Data2D<float>>(ntr, 1, 1.0, 0.0);
            Tres2D->copyCoords(Shotgeom);
            Tres2D->makeMap(lmodel->getGeom());

            float *res = Tres2D->getData();
            for (int i=0; i<ntr;i++){
               res[i] = 1;
            }
            rays->insertResiduals(Tres2D, GMAP);
            rays->solve_adj();

            // Creating gradient objects
            rayden = std::make_shared<rockseis::Image2D<float>>(Raydensfile + "-" + std::to_string(work.id), lmodel, 1, 1);

            rayden->allocateImage();
            float *grad = rayden->getImagedata();
            float *lam = rays->getLam();

            int nx, nz;
            int nx_pml, nz_pml;
            int lpml = lmodel->getLpml();
            nx = lmodel->getNx();
            nz = lmodel->getNz();
            nx_pml = lmodel->getNx_pml();
            nz_pml = lmodel->getNz_pml();
            Index Ilam(nx_pml, nz_pml);
            Index Igrad(nx, nz);
            for (int i=0; i<nx; i++){
               for (int j=0; j<nz; j++){
                  if(isnan(lam[Ilam(i+lpml, j+lpml)]) || isinf(lam[Ilam(i+lpml, j+lpml)]) ){
                     grad[Igrad(i,j)] = 0.0;
                  }else{
                     grad[Igrad(i,j)] = lam[Ilam(i+lpml,j+lpml)]/ntr;
                  }

                  if(j == 0 || j == nz-1){
                     grad[Igrad(i,j)] = 0.0;
                  }
                  if(i == 0 || i == nx-1){
                     grad[Igrad(i,j)] = 0.0;
                  }
               }
            }

            // Stack gradient
            rayden->write();


            // Reset all classes
            source.reset();
            rays.reset();
            Shotgeom.reset();
            lmodel.reset();
            rayden.reset();
            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }
   }
}

