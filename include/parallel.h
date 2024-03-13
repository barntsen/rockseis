#ifndef PARALLEL_H
#define PARALLEL_H

// Include statements
#include <vector>
#include <mpi.h>
#include <cstddef>
#include <memory>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"

#define WORK_NOT_STARTED 0
#define WORK_RUNNING 1
#define WORK_FINISHED 2
#define WORK_FAILED 3
#define PARALLEL_IO 5
#define MPI_TAG_DIE 10
#define MPI_TAG_NO_WORK 11


namespace rockseis {
// =============== MPI WORK STRUCTS =============== //
typedef struct {
	long int id;
	int status;
	int MPItag;
	time_t start;
	time_t end;
}  workModeling_t;

typedef struct {
	long int id;
	int status;
	int fromRank;
	int MPItag;
}  workModresult_t;

typedef struct {
	long int id;
	long int gath;
	int subgath;
	int nsubgath;
    int gathmaster;
	int status;
	int MPItag;
	time_t start;
	time_t end;
}  workDecomp_t;

typedef struct {
	long int id;
	int status;
	int fromRank;
	int MPItag;
}  workDecompresult_t;



// =============== PARENT MPI CLASS =============== //
/** The parent MPI class
 *
 */
class MPI {
public:
	MPI();					///< Default constructor
	MPI(int *argc, char ***argv);		///< Constructor with input parameters
	virtual ~MPI();					///< Destructor

	// Get functions
	int getNrank();				///< Get number of processors
	int getRank();				///< Get rank for current processor
    char *getName();			///< Get rank for current processor
    std::string getLogfile() { return logfile; }    ///< Get log file name

	// Set functions
    void setLogfile(std::string name) { logfile = name; }    ///< Set log file name
    void clearLogfile() { logfile.clear(); }    ///< Clear log file name
    void setVerbose(bool val) { verbose = val; } ///< Set verbose on
    bool getVerbose() { return verbose; } ///< Set verbose off

    void finalize() { MPI_Finalize(); } ///< Finalize MPI

	// Send and receive functions
	virtual void performWork() = 0;
	void stopSlaves();			///< Stop all slaves
	void sendNoWork(const int rank);	///< Send no work tag to rank
private:
	// Variables
	int nrank;	// Number of ranks
	int rank;	// Rank for current rank
    char name[MPI_MAX_PROCESSOR_NAME]; // Processors name
    std::string logfile; // Logfile
    bool verbose;
	
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
    void clearWork() { work.clear(); } ///< Clears all elements of work
private:
	// Variables
	MPI_Datatype MPIwork;					// MPI type for work
	MPI_Datatype MPIresult;					// MPI type for result
	std::vector<std::shared_ptr<workModeling_t>> work;	// Vector of work pointers
	
	// Functions
	void initTypes();							///< Initialize MPI types
	void sendWork(std::shared_ptr<workModeling_t> work, const int rank);	///< Send work to rank
	void sendWorkToAll();							///< Send work to all ranks
	workModresult_t receiveResult();						///< Receive result from slaves
	workModresult_t receiveResult(const int rank);				///< Receive result from a particular slave
	void checkResult(workModresult_t result);					///< Check result and update work queue
	std::shared_ptr<workModeling_t> getWork();				///< Get work that has not started
	unsigned long int getJobsleft();             ///< Return number of jobs left
};

// =============== MODELING MPI CLASS =============== //
/** The modeling MPI class
 *
 */
class MPIdecomp: public MPI {
public:
	MPIdecomp();						///< Default constructor
	MPIdecomp(int *argc, char ***argv);			///< Constructor with input parameters
	~MPIdecomp();						///< Destructor

	// Work related functions
	void addWork(std::shared_ptr<workDecomp_t> _work);	///< Add work
	void printWork();					///< Print work queue (for debugging)

	// Send and receive functions
	void performWork();					///< Perform all work in work queue
	workDecomp_t receiveWork();				///< Receive work from slave
	void sendResult(workDecomp_t _work);			///< Send result to master
    void clearWork() { work.clear(); } ///< Clears all elements of work
private:
	// Variables
	MPI_Datatype MPIwork;					// MPI type for work
	MPI_Datatype MPIresult;					// MPI type for result
	std::vector<std::shared_ptr<workDecomp_t>> work;	// Vector of work pointers
	
	// Functions
	void initTypes();							///< Initialize MPI types
	void sendWork(std::shared_ptr<workDecomp_t> work, const int rank);	///< Send work to rank
	void sendWorkToAll();							///< Send work to all ranks
	workDecompresult_t receiveResult();						///< Receive result from slaves
	workDecompresult_t receiveResult(const int rank);				///< Receive result from a particular slave
	void checkResult(workDecompresult_t result);					///< Check result and update work queue
	std::shared_ptr<workDecomp_t> getWork();				///< Get work that has not started
	unsigned long int getJobsleft();             ///< Return number of jobs left
};


}

#endif //PARALLEL_H
