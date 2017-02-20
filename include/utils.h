#ifndef UTILS_H
#define UTILS_H

// Include statements
#include <vector>
#include <string>

namespace rockseis {
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

// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template struct Point2D<int>;
template struct Point2D<float>;
template struct Point3D<int>;
template struct Point3D<float>;
}
#endif //UTILS_H
