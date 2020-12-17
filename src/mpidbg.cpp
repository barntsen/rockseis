#include <iostream>
#include "parallel.h"
#include <mpi.h>

using namespace rockseis;

int main(int argc, char** argv) {
   // Initializing MPI
   MPIdomaindecomp mpi = MPIdomaindecomp(&argc,&argv);
   mpi.setNdomain(2);
   mpi.splitDomains();
   //std::cerr << "Numer of processors: " << mpi.getNrank() << std::endl;
   //std::cerr << "Numer of domains: " << mpi.getNdomain() << std::endl;

   if(mpi.getRank() == 0) {
      // Master

      // Create work queue
      for(long int i=0; i<3; i++) {
         // Work struct
         std::shared_ptr<workModeling_t> work = std::make_shared<workModeling_t>(workModeling_t{i,WORK_NOT_STARTED});
         mpi.addWork(work);
      }

      // Print work queue
      std::cerr << "Work queue before parallelization" << std::endl;
      mpi.printWork();

      // Perform work in parallel
      mpi.performWork();

      // Print work queue
      std::cerr << "Work queue after parallelization" << std::endl;
      mpi.printWork();
   }
   else {
      // Domain Masters
      while(1) {
         if(!mpi.ifActive()){
            std::cerr << "My rank in WORLD: " << mpi.getRank() << " I am not active." << std::endl;
            break;
         }
         workModeling_t work = mpi.receiveWork();


         if(work.MPItag == MPI_TAG_DIE) {
            //mpi.stopDomainSlaves();
            break;
         }

         if(work.MPItag == MPI_TAG_NO_WORK) {
            mpi.sendNoWork(mpi.getMasterComm(), 0);
         }
         else {
            // Do some work
            std::cerr << "Received order to get shot number " << work.id << std::endl;
            std::cerr << "Processing work with id number " << work.id << std::endl;
            work.status = WORK_FINISHED;

            // Send result back
            mpi.sendResult(work);		
         }
      }

   }

   //std::cerr << "My rank in WORLD: " << mpi.getRank() << " My rank in domain: " << mpi.getDomainrank() << " My rank in master is: " << mpi.getMasterrank() << std::endl;
   std::cerr << "My rank in WORLD: " << mpi.getRank() << " I was killed." << std::endl;

}

