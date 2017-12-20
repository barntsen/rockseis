#ifndef UTILS_H
#define UTILS_H

// Include statements
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#define NHEAD2D 4
#define NHEAD3D 6
#define SMAP 0
#define GMAP 1

#define PRINT_DOC(s) std::cout << (#s) << std::endl;

namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PRESSURE, VX, VY, VZ, SXX, SYY, SZZ, SYZ, SXZ, SXY} rs_field; ///< What kind of data is recorded in file
typedef enum {SOURCE, RECEIVER, CMP} rs_key; ///< Information on how data is sorted
typedef enum {FINISHED, RUNNING, NOT_STARTED, FAILED} rs_status; ///< Status of a process(ex. in modelling, migration, ...)
typedef enum {LINEAR, BSPLINE, SINC} rs_interpmode; ///< Interpolation mode.)
typedef enum {DIFFERENCE, CORRELATION, ADAPTIVE_GAUSS, ADAPTIVE_LINEAR} rs_fwimisfit; ///< Type of misfit for FWI.)

/** The Progress struct
 *
 */
typedef struct{
	clock_t previous, current; ///< Time book keeping
	float persec; ///< Iterations per second 
	char speed[48]; ///< Iterations per second string
	char progress[128]; ///< Progress string
} Progress;

// =============== INDEX CLASS =============== //
/** The Index class.
 * Index for seeking memory in 2D and 3D arrays.
 * */
class Index{
public:
    Index(long int _nx, long int _ny); ///< 2D constructor
    Index(long int _nx, long int _ny, long int _nz); ///< 3D constructor
    Index(long int _nx, long int _ny, long int _nhx, long int _nhy); ///< 4D constructor
    Index(long int _nx, long int _ny, long int _nz, long int _nhx, long int _nhy, long int _nhz); ///< 6D constructor
    ~Index(); ///< Destructor
    long int operator() (long int ix, long int iy) { return (iy*nx + ix); } ///< 1D operator
    long int operator() (long int ix, long int iy, long int iz) { return (iz*ny*nx + iy*nx + ix); } ///< 3D operator 
    long int operator() (long int ix, long int iy, long int ihx, long int ihy) { return (ihy*nhx*ny*nx + ihx*ny*nx + iy*nx + ix); } ///< 4D operator 
    long int operator() (long int ix, long int iy, long int iz, long int ihx, long int ihy, long int ihz) { return (ihz*nhy*nhx*nz*ny*nx + ihy*nhx*nz*ny*nx + ihx*nz*ny*nx + iz*ny*nx + iy*nx + ix); } ///< 6D operator 
private:
    long int nx; 
    long int ny;
    long int nz;
    long int nhx;
    long int nhy;
    long int nhz;
};

// =============== POINT STRUCT =============== //
/** A 2D point struct.
 * Index for seeking memory in 2D and 3D arrays.
 * */
template <typename T> 
struct Point2D { 
	T x; ///< x coordinate
	T y; ///< y coordinate
}; 

/** A 3D point struct.
 * Index for seeking memory in 2D and 3D arrays.
 * */
template <typename T> 
struct Point3D { 
	T x; ///< x coordinate
	T y; ///< y coordinate
	T z; ///< z coordinate
}; 

// =============== Error handling =============== //
/** A fatal error message
 * 
 * */
void rs_error(std::string);
void rs_error(std::string, std::string);
void rs_error(std::string, std::string, std::string);

void rs_warning(std::string);
void rs_warning(std::string, std::string);
void rs_warning(std::string, std::string, std::string);

void remove_file(std::string);

}
#endif //UTILS_H
