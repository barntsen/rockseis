#include <iostream>
#include "parallel.h"
#include <mpi.h>

using namespace rockseis;

int main(int argc, char** argv) {
	// Initializing MPI
	MPImodeling mpi = MPImodeling(&argc,&argv);
	
	if(mpi.getRank() == 0) {
		// Master

		// Create work queue
		for(int i=0; i<20; i++) {
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
		// Slave
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
				work.status = WORK_FINISHED;

				// Send result back
				mpi.sendResult(work);		
			}
		}
	}
}

