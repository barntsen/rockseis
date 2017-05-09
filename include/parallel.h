#ifndef PARALLEL_H
#define PARALLEL_H

// Include statements
#include <vector>
#include <mpi.h>
#include <memory>
#include <iostream>

#define WORK_NOT_STARTED 0
#define WORK_FINISHED 1
#define WORK_RUNNING 2
#define WORK_FAILED 3
#define MPI_TAG_DIE 10
#define MPI_TAG_NO_WORK 11


namespace rockseis {
// =============== MPI WORK STRUCTS =============== //
typedef struct {
	unsigned long id;
	int status;
	int MPItag;
}  workModeling_t;

typedef struct {
	unsigned long int id;
	int status;
	int fromRank;
	int MPItag;
}  workResult_t;


// =============== PARENT MPI CLASS =============== //
/** The parent MPI class
 *
 */
class MPI {
public:
	MPI();					///< Default constructor
	MPI(int *argc, char ***argv);		///< Constructor with input parameters
	~MPI();					///< Destructor

	// Get functions
	int getNrank();				///< Get number of processors
	int getRank();				///< Get rank for current processor

	// Send and receive functions
	virtual void performWork() = 0;
	void stopSlaves();			///< Stop all slaves
	void sendNoWork(const int rank);	///< Send no work tag to rank
private:
	// Variables
	int nrank;	// Number of ranks
	int rank;	// Rank for current rank
	
	// Functions
	virtual void sendWorkToAll() = 0;	///< Send work to all slaves
	virtual void initTypes() = 0;		///< Initialize MPI types for communication
};


// =============== MODELING MPI CLASS =============== //
/** The modeling MPI class
 *
 */
class MPImodeling: public MPI {
public:
	MPImodeling();						///< Default constructor
	MPImodeling(int *argc, char ***argv);			///< Constructor with input parameters
	~MPImodeling();						///< Destructor

	// Work related functions
	void addWork(std::shared_ptr<workModeling_t> _work);	///< Add work
	void printWork();					///< Print work queue (for debugging)

	// Send and receive functions
	void performWork();					///< Perform all work in work queue
	workModeling_t receiveWork();				///< Receive work from slave
	void sendResult(workModeling_t _work);			///< Send result to master
private:
	// Variables
	MPI_Datatype MPIwork;					// MPI type for work
	MPI_Datatype MPIresult;					// MPI type for result
	std::vector<std::shared_ptr<workModeling_t>> work;	// Vector of work pointers
	
	// Functions
	void initTypes();							///< Initialize MPI types
	void sendWork(std::shared_ptr<workModeling_t> work, const int rank);	///< Send work to rank
	void sendWorkToAll();							///< Send work to all ranks
	workResult_t receiveResult();						///< Receive result from slaves
	void checkResult(workResult_t result);					///< Check result and update work queue
	std::shared_ptr<workModeling_t> getWork();				///< Get work that has not started
};

}

#endif //PARALLEL_H
