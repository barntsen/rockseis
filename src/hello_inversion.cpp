#include <iostream>
#include <mpi.h>
#include "inversion.h"

using namespace rockseis;

int main(int argc, char** argv) {
    int status;
    Inversion<float> inv = Inversion<float>(&argc,argv);
    status = inv.runAcousticfwigrad2d();

    std::cerr << status << std::endl;
    return 0;
}

