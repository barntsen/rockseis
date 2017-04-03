#include "utils.h"

namespace rockseis {
// constructor
Index::Index(int _nx, int _ny, int _nz)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
}
Index::Index(int _nx, int _ny)
{
    nx = _nx;
    ny = _ny;
    nz = 1;
}

// () operator
int Index::operator() (int ix, int iy, int iz)
{
    return (iz*ny*nx + iy*nx + ix);
}
int Index::operator() (int ix, int iy)
{
    return (iy*nx + ix); 
}

// destructor
Index::~Index(){
    /* Do nothing*/
}

// Fatal error message 
void rs_error(std::string msg){
    std::cerr << "Fatal error: " << msg << std::endl;
    exit(1);
}
}


