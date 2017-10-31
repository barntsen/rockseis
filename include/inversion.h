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

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

namespace rockseis {

// ##### INVERSION CLASS
template<typename T>
class Inversion {
public:
    Inversion(); ///<Constructor
    ~Inversion(); ///<Destructor
    void runAcousticfwigrad2d(std::shared_ptr<MPImodeling> mpi);

private:
    std::string fwicfg;
};

}
#endif //INVERSION_H
