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

#define PRINT_DOC(s) std::cerr << (#s) << std::endl;

namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PRESSURE, VX, VY, VZ, SXX, SYY, SZZ, SYZ, SXZ, SXY} rs_field; ///< What kind of data is recorded in file
typedef enum {SOURCE, RECEIVER, CMP} rs_key; ///< Information on how data is sorted
typedef enum {GENERIC, REGULAR, DATA2D, DATA3D, SNAPSHOT, KEYMAP, SORTMAP} rs_datatype; ///< Information about the file content (Regular model, 2D data, 3D data, etc.)
typedef enum {FINISHED, RUNNING, NOT_STARTED, FAILED} rs_status; ///< Status of a process(ex. in modelling, migration, ...)
typedef enum {LINEAR, BSPLINE, SINC} rs_interpmode; ///< Interpolation mode.)

// =============== INDEX CLASS =============== //
/** The Index class.
 * Index for seeking memory in 2D and 3D arrays.
 * */
class Index{
public:
    Index(int _nx, int _ny); ///< 2D constructor
    Index(int _nx, int _ny, int _nz); ///< 3D constructor
    ~Index(); ///< Destructor
    int operator() (int ix, int iy); ///< 2D operator
    int operator() (int ix, int iy, int iz); ///< 3D operator 
private:
    int nx; 
    int ny;
    int nz;
    
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

}
#endif //UTILS_H
