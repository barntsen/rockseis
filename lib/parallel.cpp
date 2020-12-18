#include "parallel.h"

namespace rockseis {
// =============== PARENT MPI CLASS =============== //
MPI::MPI() {
	MPI_Init(NULL,NULL);

	// Getting size for rank
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);
	
	// Getting rank for current processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Getting name for current processors
    int length;
    MPI_Get_processor_name(name, &length);
    namelen = length;

    logfile = "mpiqueue.log";
    verbose = true;


}

MPI::MPI(int *argc, char ***argv) {
	MPI_Init(argc,argv);
	
	// Getting size for rank
	MPI_Comm_size(MPI_COMM_WORLD,&nrank);

	// Getting rank for current processors
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Getting name for current processors
    int length;
    MPI_Get_processor_name(name, &length);
    namelen = length;

    logfile = "mpiqueue.log";
    verbose = true;

}


MPI::~MPI() {
	// Nothing
    MPI_Finalize();
}

int MPI::getNrank() {
	return nrank;
}

char *MPI::getName() {
	return name;
}

int MPI::getRank() {
	return rank;
}

void MPI::stopSlaves() {
	for(int i=1; i<nrank; i++) {
		MPI_Send(0,0,MPI_INT,i,MPI_TAG_DIE,MPI_COMM_WORLD);
	}
}

void MPI::stopSlaves(MPI_Comm *Comm) {
   int nslaves;
   MPI_Comm_size(*Comm,&nslaves);
   for(int i=1; i<nslaves; i++) {
      MPI_Send(0,0,MPI_INT,i,MPI_TAG_DIE,*Comm);
   }
}

void MPI::sendNoWork(const int rank) {
   MPI_Send(0,0,MPI_INT,rank,MPI_TAG_NO_WORK,MPI_COMM_WORLD);
}

void MPI::sendNoWork(MPI_Comm *Comm, const int rank) {
   if(*Comm != MPI_COMM_NULL){
      MPI_Send(0,0,MPI_INT,rank,MPI_TAG_NO_WORK,*Comm);
   }
}




// =============== MODELING MPI CLASS =============== //
MPImodeling::MPImodeling(): MPI() {
   initTypes();
}

MPImodeling::MPImodeling(int *argc, char ***argv): MPI(argc,argv) {
   initTypes();
}

MPImodeling::~MPImodeling() {
}

void MPImodeling::initTypes() {
   // Work type
   int count = 5;
   int lengths[5] = {1,1,1,1,1};
   MPI_Datatype types[6] = {MPI_LONG_INT,MPI_INT,MPI_INT,MPI_UNSIGNED_LONG,MPI_UNSIGNED_LONG};
   MPI_Aint offsets[6];
   offsets[0] = offsetof(workModeling_t,id);
   offsets[1] = offsetof(workModeling_t,status);
   offsets[2] = offsetof(workModeling_t,MPItag);
   offsets[3] = offsetof(workModeling_t,start);
   offsets[4] = offsetof(workModeling_t,end);
   MPI_Type_create_struct(count,lengths,offsets,types,&MPIwork);
   MPI_Type_commit(&MPIwork);

   // Result type
   int countR = 4;
   int lengthsR[4] = {1,1,1,1};
   MPI_Datatype typesR[4] = {MPI_LONG_INT,MPI_INT,MPI_INT,MPI_INT};
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

unsigned long int MPImodeling::getJobsleft(){
   int jobs_left = work.size();
   for(auto const& p: work) {
      if(p->status == WORK_FINISHED) {
         jobs_left--;
      }
   }

   return jobs_left;
}

void MPImodeling::addWork(std::shared_ptr<workModeling_t> _work) {
   // Set start and end times to 0
   _work->start=0;
   _work->end=0;
   // Add to queue
   work.push_back(_work);
}

void MPImodeling::printWork() {
   std::string logfilename = this->getLogfile(); 

   // Variables
   char buffer[256],start[256],end[256];
   time_t runtime,now;

   if(!logfilename.empty() && (this->getVerbose() == true)){
      std::ofstream Flog; 
      Flog.open(logfilename);
      if(!Flog.is_open()) rs_error("MPImodeling::printWork: Error writting to MPI log file.");

      // Seeking to beginning of file
      Flog.seekp(0);

      // Title in file
      Flog << "*******************************\n";
      Flog << "*                             *\n";
      Flog << "*    MPI QUEUE STATUS FILE    *\n";
      Flog << "*                             *\n";
      Flog << "*******************************\n";
      Flog << "\n";

      // CPU information
      Flog << "CPU INFORMATION\n";
      // Number of CPUs
      snprintf(buffer,256,"# CPU: %d \n",this->getNrank());
      Flog << buffer;
      Flog << "\n";

      // Job information
      Flog << "JOB INFORMATION\n";
      snprintf(buffer,256,"Jobs: %lu, #Remaining jobs: %lu\n\n",work.size(),this->getJobsleft());
      Flog << buffer;
      Flog << "#JobID: CPUID Status (0=not started, 1=running, 2=finished). \n";
      if(work.size() == 0) {
         Flog << "No work in queue\n";
      }
      else {
         for(auto const& p: work) {
            // Converting time to readable format	
            strcpy(start,ctime(&p->start));
            strcpy(end,ctime(&p->end));
            start[strlen(start)-1] = '\0';
            end[strlen(end)-1] = '\0';

            // Writing beginning of status line
            snprintf(buffer,256,"#%04lu:",p->id);  
            Flog << buffer;
            // Writting rank 
            if(p->MPItag == 0) 
            {
               Flog << "  N/A ";	
            }else{
               snprintf(buffer,256,"  %04d", p->MPItag);
               Flog << buffer;
            }
            snprintf(buffer, 256, "  %d  start: ", p->status);
            Flog << buffer;
            // Start time
            if(p->start == 0) Flog << "N/A";	
            else Flog << start;

            // End time
            Flog << ", end: ";
            if(p->end == 0) 	Flog << "N/A";
            else Flog << end;

            // Runtime
            Flog << ", runtime: ";
            if(p->start == 0 && p->end == 0) {
               Flog << "N/A";
            }
            else {
               // Calculating runtime
               if(p->end <= 0) runtime = time(NULL) - p->start;
               else runtime = p->end - p->start;

               if(runtime < 600) {
                  // Less than 10min, printing out sec and min
                  snprintf(buffer,256,"%ld",runtime);
                  Flog << buffer;
                  Flog << " sec (";
                  snprintf(buffer,256,"%ld",runtime/60);
                  Flog << buffer;
                  Flog << " min)";
               }
               else if(runtime >= 600 && runtime < (3600*12)) {
                  // Less than 12 hours, printing out min and hour
                  snprintf(buffer,256,"%5.2f", (float) runtime/60);
                  Flog << buffer;
                  Flog << " min (";
                  snprintf(buffer,256,"%4.2f",(float) runtime/3600);
                  Flog << buffer;
                  Flog << " h)";
               }
               else {
                  // Long runtimes, printing out hour and days
                  snprintf(buffer,256,"%5.2f",(float) runtime/3600);
                  Flog << buffer;
                  Flog << " h (";
                  snprintf(buffer,256,"%4.2f",(float) runtime/86400);
                  Flog << buffer;
                  Flog << " d)";
               }

               // Printing out running flag if necessary
               if(p->end <= 0) {
                  Flog << " (r)";
               }
            }

            // New line
            Flog << "\n";
         }

         // Details summary at end of file
         now = time(NULL);
         strcpy(start,ctime(&now));
         start[strlen(start)-1] = '\0';
         snprintf(buffer,256,"Jobs: %lu, #Remaining jobs: %lu\n\n",work.size(), this->getJobsleft());
         Flog << buffer;
         Flog << "Last updated: "; Flog << start; Flog << "\n\n";
         Flog.close();
      }

   }
}

void MPImodeling::sendWorkToAll() {
   for(int i=1; i<getNrank(); i++) {
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
   //MPI_Status status;
   // Changing status on work and sending to slave
   work->status = WORK_RUNNING;
   work->start = time(NULL);
   work->MPItag = rank;
   MPI_Send(work.get(),1,MPIwork,rank,0,MPI_COMM_WORLD);
}

void MPImodeling::performWork() {
   // Sending initial work to all ranks
   sendWorkToAll();

   // Working through queue
   std::shared_ptr<workModeling_t> work = getWork();
   // Print initial queue	
   printWork();
   while(work != NULL) {
      // Receive work from slaves
      workResult_t result = receiveResult();

      if(result.MPItag != MPI_TAG_NO_WORK) {
         // Check result
         checkResult(result);
         // Print updated queue
         printWork();
      }

      // Sending work to slave
      sendWork(work,result.fromRank);
      printWork();

      // Getting new job
      work = getWork();
   }

   // Work queue is finished so collecting results
   for(int i=1; i<getNrank(); ++i) {
      // Receive work from slaves
      workResult_t result = receiveResult();

      if(result.MPItag != MPI_TAG_NO_WORK) {
         // Check result
         checkResult(result);
         // Print updated queue
         printWork();
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
   result.fromRank = this->getRank();
   result.MPItag = _work.MPItag;

   MPI_Send(&result,1,MPIresult,0,0,MPI_COMM_WORLD);
}

workResult_t MPImodeling::receiveResult() {
   // Variables
   workResult_t result;
   MPI_Status status;
   // Receiving 
   MPI_Recv(&result,1,MPIresult,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   // Updating struct
   result.MPItag = status.MPI_TAG;

   return result;
}

workResult_t MPImodeling::receiveResult(const int rank) {
   // Variables
   workResult_t result;
   MPI_Status status;
   // Receiving 
   MPI_Recv(&result,1,MPIresult,rank,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
   // Updating struct
   result.MPItag = status.MPI_TAG;

   return result;
}

void MPImodeling::checkResult(workResult_t result) {
   int rank;
   switch(result.status){
      case PARALLEL_IO:
         rank = result.fromRank;
         workResult_t result2;
         result2 = receiveResult(rank);
         if(result2.status != WORK_FINISHED) {
            if(result2.fromRank != 0){
               // Job failed
               work[result2.id]->status = WORK_NOT_STARTED;
            }
         }
         else {
            // Job finished successfully
            work[result2.id]->status = WORK_FINISHED;
            work[result2.id]->end = time(NULL);
            work[result2.id]->MPItag = result2.fromRank;
         }
         break;
      case WORK_FINISHED:
         // Job finished successfully
         work[result.id]->status = WORK_FINISHED;
         work[result.id]->end = time(NULL);
         work[result.id]->MPItag = result.fromRank;
         break;
      default:
         if(result.fromRank != 0){
            // Job failed
            work[result.id]->status = WORK_NOT_STARTED;
         }
         break;
   }
}


// =============== DOMAIN DECOMPOSITION MPI CLASS =============== //
MPIdomaindecomp::MPIdomaindecomp(): MPI() {
   initTypes();

   ndomain = 1;
   masterrank = 0;
   domainrank = 0;
}

MPIdomaindecomp::MPIdomaindecomp(int *argc, char ***argv): MPI(argc,argv) {
   initTypes();


   ndomain = 1;
   masterrank = 0;
   domainrank = 0;
   MPI_COMM_DOMAIN = MPI_COMM_NULL;
   MPI_COMM_MASTERS = MPI_COMM_NULL;
}

MPIdomaindecomp::~MPIdomaindecomp() {
}

void MPIdomaindecomp::initTypes() {
   // Work type
   int count = 5;
   int lengths[5] = {1,1,1,1,1};
   MPI_Datatype types[6] = {MPI_LONG_INT,MPI_INT,MPI_INT,MPI_UNSIGNED_LONG,MPI_UNSIGNED_LONG};
   MPI_Aint offsets[6];
   offsets[0] = offsetof(workModeling_t,id);
   offsets[1] = offsetof(workModeling_t,status);
   offsets[2] = offsetof(workModeling_t,MPItag);
   offsets[3] = offsetof(workModeling_t,start);
   offsets[4] = offsetof(workModeling_t,end);
   MPI_Type_create_struct(count,lengths,offsets,types,&MPIwork);
   MPI_Type_commit(&MPIwork);

   // Result type
   int countR = 4;
   int lengthsR[4] = {1,1,1,1};
   MPI_Datatype typesR[4] = {MPI_LONG_INT,MPI_INT,MPI_INT,MPI_INT};
   MPI_Aint offsetsR[4];
   offsetsR[0] = offsetof(workResult_t,id);
   offsetsR[1] = offsetof(workResult_t,status);
   offsetsR[2] = offsetof(workResult_t,fromRank);
   offsetsR[3] = offsetof(workResult_t,MPItag);
   MPI_Type_create_struct(countR,lengthsR,offsetsR,typesR,&MPIresult);
   MPI_Type_commit(&MPIresult);
}

std::shared_ptr<workModeling_t> MPIdomaindecomp::getWork(){
   for(auto const& p: work) {
      if(p->status == WORK_NOT_STARTED) {
         return p;
      }
   }
   return NULL;
}

unsigned long int MPIdomaindecomp::getJobsleft(){
   int jobs_left = work.size();
   for(auto const& p: work) {
      if(p->status == WORK_FINISHED) {
         jobs_left--;
      }
   }

   return jobs_left;
}

void MPIdomaindecomp::addWork(std::shared_ptr<workModeling_t> _work) {
   // Set start and end times to 0
   _work->start=0;
   _work->end=0;
   // Add to queue
   work.push_back(_work);
}

void MPIdomaindecomp::printWork() {
   std::string logfilename = this->getLogfile(); 

   // Variables
   char buffer[256],start[256],end[256];
   time_t runtime,now;

   if(!logfilename.empty() && (this->getVerbose() == true)){
      std::ofstream Flog; 
      Flog.open(logfilename);
      if(!Flog.is_open()) rs_error("MPIdomaindecomp::printWork: Error writting to MPI log file.");

      // Seeking to beginning of file
      Flog.seekp(0);

      // Title in file
      Flog << "*******************************\n";
      Flog << "*                             *\n";
      Flog << "*    MPI QUEUE STATUS FILE    *\n";
      Flog << "*                             *\n";
      Flog << "*******************************\n";
      Flog << "\n";

      // CPU information
      Flog << "CPU INFORMATION\n";
      // Number of CPUs
      snprintf(buffer,256,"# CPU: %d \n",this->getNrank());
      Flog << buffer;
      Flog << "\n";

      // Job information
      Flog << "JOB INFORMATION\n";
      snprintf(buffer,256,"Jobs: %lu, #Remaining jobs: %lu\n\n",work.size(),this->getJobsleft());
      Flog << buffer;
      Flog << "#JobID: CPUID Status (0=not started, 1=running, 2=finished). \n";
      if(work.size() == 0) {
         Flog << "No work in queue\n";
      }
      else {
         for(auto const& p: work) {
            // Converting time to readable format	
            strcpy(start,ctime(&p->start));
            strcpy(end,ctime(&p->end));
            start[strlen(start)-1] = '\0';
            end[strlen(end)-1] = '\0';

            // Writing beginning of status line
            snprintf(buffer,256,"#%04lu:",p->id);  
            Flog << buffer;
            // Writting rank 
            if(p->MPItag == 0) 
            {
               Flog << "  N/A ";	
            }else{
               snprintf(buffer,256,"  %04d", p->MPItag);
               Flog << buffer;
            }
            snprintf(buffer, 256, "  %d  start: ", p->status);
            Flog << buffer;
            // Start time
            if(p->start == 0) Flog << "N/A";	
            else Flog << start;

            // End time
            Flog << ", end: ";
            if(p->end == 0) 	Flog << "N/A";
            else Flog << end;

            // Runtime
            Flog << ", runtime: ";
            if(p->start == 0 && p->end == 0) {
               Flog << "N/A";
            }
            else {
               // Calculating runtime
               if(p->end <= 0) runtime = time(NULL) - p->start;
               else runtime = p->end - p->start;

               if(runtime < 600) {
                  // Less than 10min, printing out sec and min
                  snprintf(buffer,256,"%ld",runtime);
                  Flog << buffer;
                  Flog << " sec (";
                  snprintf(buffer,256,"%ld",runtime/60);
                  Flog << buffer;
                  Flog << " min)";
               }
               else if(runtime >= 600 && runtime < (3600*12)) {
                  // Less than 12 hours, printing out min and hour
                  snprintf(buffer,256,"%5.2f", (float) runtime/60);
                  Flog << buffer;
                  Flog << " min (";
                  snprintf(buffer,256,"%4.2f",(float) runtime/3600);
                  Flog << buffer;
                  Flog << " h)";
               }
               else {
                  // Long runtimes, printing out hour and days
                  snprintf(buffer,256,"%5.2f",(float) runtime/3600);
                  Flog << buffer;
                  Flog << " h (";
                  snprintf(buffer,256,"%4.2f",(float) runtime/86400);
                  Flog << buffer;
                  Flog << " d)";
               }

               // Printing out running flag if necessary
               if(p->end <= 0) {
                  Flog << " (r)";
               }
            }

            // New line
            Flog << "\n";
         }

         // Details summary at end of file
         now = time(NULL);
         strcpy(start,ctime(&now));
         start[strlen(start)-1] = '\0';
         snprintf(buffer,256,"Jobs: %lu, #Remaining jobs: %lu\n\n",work.size(), this->getJobsleft());
         Flog << buffer;
         Flog << "Last updated: "; Flog << start; Flog << "\n\n";
         Flog.close();
      }
   }
}

void MPIdomaindecomp::sendWorkToAll() {
   for(int i=1; i<getNmaster(); i++) {
      // Get work from queue
      std::shared_ptr<workModeling_t> work = getWork();
      if(work != NULL) {
         sendWork(work,i);
      }
      else {
         sendNoWork(&MPI_COMM_MASTERS,i);
      }
   }
}

void MPIdomaindecomp::sendWork(std::shared_ptr<workModeling_t> work, const int rank) {
   //MPI_Status status;
   // Changing status on work and sending to slave
   if(MPI_COMM_MASTERS != MPI_COMM_NULL){
      work->status = WORK_RUNNING;
      work->start = time(NULL);
      work->MPItag = rank;
      MPI_Send(work.get(),1,MPIwork,rank,0,MPI_COMM_MASTERS);
   }
}

void MPIdomaindecomp::performWork() {
   // Sending initial work to all ranks
   sendWorkToAll();

   // Working through queue
   std::shared_ptr<workModeling_t> work = getWork();
   // Print initial queue	
   printWork();
   while(work != NULL) {
      // Receive work from slaves
      workResult_t result = receiveResult();

      if(result.MPItag != MPI_TAG_NO_WORK) {
         // Check result
         checkResult(result);
         // Print updated queue
         printWork();
      }

      // Sending work to slave
      sendWork(work,result.fromRank);
      printWork();

      // Getting new job
      work = getWork();
   }

   // Work queue is finished so collecting results
   for(int i=1; i<getNmaster(); ++i) {
      // Receive work from slaves
      workResult_t result = receiveResult();

      if(result.MPItag != MPI_TAG_NO_WORK) {
         // Check result
         checkResult(result);
         // Print updated queue
         printWork();
      }
   }

   // All work is done starting to stop each slave
   stopSlaves(&MPI_COMM_MASTERS);
}


workModeling_t MPIdomaindecomp::receiveWork() {
   // Variables
   workModeling_t work;
   MPI_Status status;
   int i;
   if(MPI_COMM_MASTERS != MPI_COMM_NULL){
      // Receiving
      MPI_Recv(&work,1,MPIwork,0,MPI_ANY_TAG,MPI_COMM_MASTERS,&status);
      // Updating struct
      work.MPItag = status.MPI_TAG;
      // Send work to Domain slaves
      for(i=1; i < ndomain; i++){
         MPI_Send(&work,1,MPIwork,i,0,MPI_COMM_DOMAIN);
      }
   }else{
      // Receiving work from Master of Domain
      MPI_Recv(&work,1,MPIwork,0,MPI_ANY_TAG,MPI_COMM_DOMAIN,&status);
   }

   return work;
}

void MPIdomaindecomp::sendResult(workModeling_t _work) {
   workResult_t result;

   if(MPI_COMM_MASTERS != MPI_COMM_NULL){
      result.id = _work.id;
      result.status = _work.status;
      result.fromRank = this->getMasterrank();
      result.MPItag = _work.MPItag;
      MPI_Send(&result,1,MPIresult,0,0,MPI_COMM_MASTERS);
   }
}

workResult_t MPIdomaindecomp::receiveResult() {
   // Variables
   workResult_t result;
   if(MPI_COMM_MASTERS != MPI_COMM_NULL){
      MPI_Status status;
      // Receiving 
      MPI_Recv(&result,1,MPIresult,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_MASTERS,&status);
      // Updating struct
      result.MPItag = status.MPI_TAG;
   }
   return result;
}

workResult_t MPIdomaindecomp::receiveResult(const int rank) {
   workResult_t result;
   if(MPI_COMM_MASTERS != MPI_COMM_NULL){
      // Variables
      MPI_Status status;
      // Receiving 
      MPI_Recv(&result,1,MPIresult,rank,MPI_ANY_TAG,MPI_COMM_MASTERS,&status);
      // Updating struct
      result.MPItag = status.MPI_TAG;
   }
   return result;
}

void MPIdomaindecomp::checkResult(workResult_t result) {
   int rank;
   switch(result.status){
      case PARALLEL_IO:
         rank = result.fromRank;
         workResult_t result2;
         result2 = receiveResult(rank);
         if(result2.status != WORK_FINISHED) {
            if(result2.fromRank != 0){
               // Job failed
               work[result2.id]->status = WORK_NOT_STARTED;
            }
         }
         else {
            // Job finished successfully
            work[result2.id]->status = WORK_FINISHED;
            work[result2.id]->end = time(NULL);
            work[result2.id]->MPItag = result2.fromRank;
         }
         break;
      case WORK_FINISHED:
         // Job finished successfully
         work[result.id]->status = WORK_FINISHED;
         work[result.id]->end = time(NULL);
         work[result.id]->MPItag = result.fromRank;
         break;
      default:
         if(result.fromRank != 0){
            // Job failed
            work[result.id]->status = WORK_NOT_STARTED;
         }
         break;
   }
}

void MPIdomaindecomp::splitDomains() {
   int nrank = this->getNrank();
   int rank = this->getRank();
   bool status;

   if(ndomain > nrank-1) rs_error("MPIdomaindecomp::splitDomains: Number domains cannot exceed number of slaves (number of processors - 1");
   nmaster = (int) ((nrank-1)/ndomain);
   int color;
   int iam = rank;
   int key;
   if(rank > 0){
      if ((iam-1) < ndomain * nmaster) {
         color = (int) ((iam-1)/ndomain); 
         status = MPI_Comm_split(MPI_COMM_WORLD, color, (iam-1), &MPI_COMM_DOMAIN);
         if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
         MPI_Comm_rank(MPI_COMM_DOMAIN, &domainrank);
      } else {
         color = MPI_UNDEFINED;
         status = MPI_Comm_split(MPI_COMM_WORLD, color, 0, &MPI_COMM_DOMAIN);
         if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
         domainrank=-1;
      }
   }else{
      color = MPI_UNDEFINED;
      status = MPI_Comm_split(MPI_COMM_WORLD, color, 0, &MPI_COMM_DOMAIN);
      if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
      domainrank=-1;
   }

   if(rank == 0) {
      nmaster = nmaster+1;
      color=1;
      status=MPI_Comm_split(MPI_COMM_WORLD, color, rank, &MPI_COMM_MASTERS);
      if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
      MPI_Comm_rank(MPI_COMM_MASTERS, &masterrank);
   }else if(domainrank == 0 && (MPI_COMM_DOMAIN != MPI_COMM_NULL)){
      color=1;
      key=(int) ((iam-1)/ndomain+1);
      status=MPI_Comm_split(MPI_COMM_WORLD, color, key, &MPI_COMM_MASTERS);
      if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
      MPI_Comm_rank(MPI_COMM_MASTERS, &masterrank);
   }else{
      color = MPI_UNDEFINED;
      status=MPI_Comm_split(MPI_COMM_WORLD, color, 0, &MPI_COMM_MASTERS);
      if(status != MPI_SUCCESS) rs_error("MPIdomaindecomp::splitDomains: Error splitting the Comms.");
      masterrank=-1;
   }
}



}
