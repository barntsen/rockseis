#include "utils.h"

namespace rockseis {
// constructor
Index::Index(long int _nx, long int _ny)
{
    nx = _nx;
    ny = _ny;
    nz = 1;
    nhx = 1;
    nhy = 1;
    nhz = 1;
}

Index::Index(long int _nx, long int _ny, long int _nz)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
    nhx = 1;
    nhy = 1;
    nhz = 1;
}

Index::Index(long int _nx, long int _ny, long int _nhx, long int _nhy)
{
    nx = _nx;
    ny = _ny;
    nz = 1;
    nhx = _nhx;
    nhy = _nhy;
    nhz = 1;
}

Index::Index(long int _nx, long int _ny, long int _nz, long int _nhx, long int _nhy, long int _nhz)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
    nhx = _nhx;
    nhy = _nhy;
    nhz = _nhz;
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

// Removing file
void remove_file(std::string filename) {
	if(!filename.empty()){
		if( remove( filename.c_str() ) != 0 ){
			rs_error( "remove_file: Error deleting file: ", filename);
		}
	}
}

bool file_exists (std::string filename) {
    if (FILE *file = fopen(filename.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
}

//Timer
double wtime() {
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return ((double)tp.tv_sec + (double)tp.tv_nsec*1.0e-9) ;
}

// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template struct Point2D<int>;
template struct Point2D<float>;
template struct Point3D<int>;
template struct Point3D<float>;

}


