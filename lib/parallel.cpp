#include "parallel.h"

namespace rockseis {
// =============== PARENT MPI CLASS =============== //
MPI::MPI() {
	MPI_Init(NULL,NULL);

	// Getting size for rank
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	
	// Getting rank for current processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

MPI::MPI(int *argc, char ***argv) {
	MPI_Init(argc,argv);
	
	// Getting size for rank
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);

	// Getting rank for current processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

MPI::~MPI() {
	MPI_Finalize();
}

int MPI::getNRank() {
	return nrank;
}

int MPI::getRank() {
	return rank;
}

void MPI::stopSlaves() {
	for(int i=0; i<nrank; i++) {
		MPI_Send(0,0,MPI_INT,i,MPI_TAG_DIE,MPI_COMM_WORLD);
	}
}

void MPI::sendNoWork(const int rank) {
	MPI_Send(0,0,MPI_INT,rank,MPI_TAG_NO_WORK,MPI_COMM_WORLD);
}


// =============== MODELING MPI CLASS =============== //
MPImodeling::MPImodeling(): MPI() {
	initTypes();
}

MPImodeling::MPImodeling(int *argc, char ***argv): MPI(argc,argv) {
	initTypes();
}

MPImodeling::~MPImodeling() {
	// Nothing
}

void MPImodeling::initTypes() {
	// Work type
	int count = 3;
	int lengths[3] = {1,1,1};
	MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT};
	MPI_Aint offsets[3];
	offsets[0] = offsetof(workModeling_t,id);
	offsets[1] = offsetof(workModeling_t,status);
	offsets[2] = offsetof(workModeling_t,MPItag);
	MPI_Type_create_struct(count,lengths,offsets,types,&MPIwork);
	MPI_Type_commit(&MPIwork);

	// Result type
	int countR = 4;
	int lengthsR[4] = {1,1,1,1};
	MPI_Datatype typesR[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	MPI_Aint offsetsR[4];
	offsetsR[0] = offsetof(workResult_t,id);
	offsetsR[1] = offsetof(workResult_t,status);
	offsetsR[2] = offsetof(workResult_t,fromRank);
	offsetsR[3] = offsetof(workResult_t,MPItag);
	MPI_Type_create_struct(countR,lengthsR,offsetsR,typesR,&MPIresult);
	MPI_Type_commit(&MPIresult);
}

std::shared_ptr<workModeling_t> MPImodeling::getWork(){
	for(auto const& p: work) {
		if(p->status == WORK_NOT_STARTED) {
			return p;
		}
	}

	return NULL;
}

void MPImodeling::addWork(std::shared_ptr<workModeling_t> _work) {
	work.push_back(_work);
}

void MPImodeling::printWork() {
	std::cerr << "=========================\n";
	std::cerr << "       WORK QUEUE\n";

	if(work.size() == 0) {
		std::cerr << "No work in queue\n";
	}
	else {
		for(auto const& p: work) {
			std::cerr << "id: " << p->id << " status: " << p->status << std::endl;
		}
	}
	std::cerr << "\n\n";
}

void MPImodeling::sendWorkToAll() {
	for(int i=1; i<getNRank(); i++) {
		// Get work from queue
		std::shared_ptr<workModeling_t> work = getWork();
		if(work != NULL) {
			sendWork(work,i);
		}
		else {
			sendNoWork(i);
		}
	}
}

void MPImodeling::sendWork(std::shared_ptr<workModeling_t> work, const int rank) {
	// Changing status on work and sending to slave
	work->status = WORK_RUNNING;
	MPI_Send(work.get(),1,MPIwork,rank,0,MPI_COMM_WORLD);
}

void MPImodeling::performWork() {
	// Sending initial work to all ranks
	sendWorkToAll();

	// Working through queue
	std::shared_ptr<workModeling_t> work = getWork();
	
	while(work != NULL) {
		// Receive work from slaves
		workResult_t result = receiveResult();

		if(result.MPItag != MPI_TAG_NO_WORK) {
			// Check result
			checkResult(result);
		}

		// Sending work to slave
		sendWork(work,result.fromRank);
		
		// Getting new job
		work = getWork();
	}

	// Work queue is finished so collecting results
	for(int i=1; i<getNRank(); ++i) {
		// Receive work from slaves
		workResult_t result = receiveResult();

		if(result.MPItag != MPI_TAG_NO_WORK) {
			// Check result
			checkResult(result);
		}
	}

	// All work is done starting to stop each slave
	stopSlaves();
}


workModeling_t MPImodeling::receiveWork() {
	// Variables
	workModeling_t work;
	MPI_Status status;
	// Receiving
	MPI_Recv(&work,1,MPIwork,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	// Updating struct
	work.MPItag = status.MPI_TAG;

	return work;
}
	
void MPImodeling::sendResult(workModeling_t _work) {
	workResult_t result;

	result.id = _work.id;
	result.status = _work.status;

	MPI_Send(&result,1,MPIresult,0,0,MPI_COMM_WORLD);

}

workResult_t MPImodeling::receiveResult() {
	// Variables
	workResult_t result;
	MPI_Status status;
	// Receiving 
	MPI_Recv(&result,1,MPIresult,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	// Updating struct
	result.fromRank = status.MPI_SOURCE;
	result.MPItag = status.MPI_TAG;

	return result;
}

void MPImodeling::checkResult(workResult_t result) {
	if(result.status != WORK_FINISHED) {
		// Job failed
		work[result.id]->status = WORK_NOT_STARTED;
	}
	else {
		// Job finished successfully
		work[result.id]->status = WORK_FINISHED;
	}
}

}
