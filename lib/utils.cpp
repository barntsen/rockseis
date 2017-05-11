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
void rs_error(std::string msg1, std::string msg2){
    std::cerr << "Fatal error: " << msg1 << msg2 << std::endl;
    exit(1);
}
void rs_error(std::string msg1, std::string msg2,  std::string msg3){
    std::cerr << "Fatal error: " << msg1 << msg2 << msg3 << std::endl;
    exit(1);
}

// Fatal error message 
void rs_warning(std::string msg){
    std::cerr << "Warning: " << msg << std::endl;
}
void rs_warning(std::string msg1, std::string msg2){
    std::cerr << "Warning: " << msg1 << msg2 << std::endl;
}
void rs_warning(std::string msg1, std::string msg2,  std::string msg3){
    std::cerr << "Warning: " << msg1 << msg2 << msg3 << std::endl;
}

// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template struct Point2D<int>;
template struct Point2D<float>;
template struct Point3D<int>;
template struct Point3D<float>;

}


