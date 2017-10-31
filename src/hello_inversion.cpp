#include <iostream>
#include <mpi.h>
#include "inversion.h"

using namespace rockseis;

//Global class 
Inversion<float> inv = Inversion<float>();

int main(int argc, char** argv) {

	// Initializing MPI
    std::shared_ptr<MPImodeling> mpi (new MPImodeling(&argc,&argv));
    int task; 

    //MASTER
    if(mpi->getRank() == 0){
        task = RUN_F_GRAD;
        std::cerr << "Rank: " << mpi->getRank() << ": Giving order to run gradient." << std::endl; 
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
        inv.runAcousticfwigrad2d(mpi);

        task = RUN_BS_PROJ;
        std::cerr << "Rank: " << mpi->getRank() << ": Giving order to run b-spline projection." << std::endl; 
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        task = BREAK_LOOP;
        std::cerr << "Rank: " << mpi->getRank() << ": Giving order to break the loop." << std::endl; 
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //SLAVE
    }else{
        bool stop = false;
        while(1)
        {
            MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
            switch(task)
            {
                case RUN_F_GRAD:
                    std::cerr << "Rank: " << mpi->getRank() << ": Running gradient." << std::endl; 
                    inv.runAcousticfwigrad2d(mpi);
                    break;
                case RUN_BS_PROJ:
                    std::cerr << "Rank: " << mpi->getRank() << ": Running B-spline projection." << std::endl; 
                    break;
                case BREAK_LOOP:
                    std::cerr << "Rank: " << mpi->getRank() << ": Leaving the loop." << std::endl; 
                    stop = true;
                    break;
            }
            if(stop){
                break;
            }
        }

        return 0;
    }
}

