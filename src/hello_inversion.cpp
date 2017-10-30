#include <iostream>
#include <mpi.h>
#include "inversion.h"

using namespace rockseis;

int main(int argc, char** argv) {

	// Initializing MPI
    std::shared_ptr<MPImodeling> mpi (new MPImodeling(&argc,&argv));

    Inversion<float> inv = Inversion<float>();
    inv.runAcousticfwigrad2d(mpi);

    return 0;
}

