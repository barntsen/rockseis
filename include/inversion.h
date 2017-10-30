#ifndef INVERSION_H
#define INVERSION_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <time.h>
#include "geometry.h"
#include "utils.h"
#include "file.h"
#include "model.h"
#include "data.h"
#include "waves.h"
#include "der.h"
#include "snap.h"
#include "bspl.h"
#include "interp.h"
#include "parallel.h"
#include "fwi.h"
#include "sort.h"
#include <inparse.h>

namespace rockseis {

// ##### INVERSION CLASS
template<typename T>
class Inversion {
public:
    Inversion(int *argc, char **argv); ///<Constructor
    ~Inversion(); ///<Destructor
    int runAcousticfwigrad2d();

private:
    std::string fwicfg;
    std::shared_ptr<MPImodeling> mpi;
};

}
#endif //INVERSION_H