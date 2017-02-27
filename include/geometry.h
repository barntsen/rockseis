#ifndef GEOMETRY_H
#define GEOMETRY_H

// Include statements
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <memory>
#include "utils.h"

#define MAXDIMS 9


namespace rockseis {
/// =============== ABSTRACT GEOMETRY CLASS =============== //
/** The geometry class
 *
 */
template<typename T>
class Geometry{
public:
    Geometry(); 		///< Constructor
    virtual ~Geometry();       	///< Destructor

    // Get functions
    size_t getN(int dim) { if( (dim <= MAXDIMS) && (dim > 0) ) return n[dim-1]; else return 0; }	///< Get dimension size
    T getD(int dim) { if( (dim <= MAXDIMS) && (dim > 0) ) return d[dim-1]; else return 0.; }	///< Get sampling interval
    T getO(int dim) { if( (dim <= MAXDIMS) && (dim > 0) ) return o[dim-1]; else return 0.; }	///< Get origin

    // Set functions
    void setN(int dim, size_t val); 	///< Set dimension size
    void setD(int dim, T val) { if( (dim <= MAXDIMS) && (dim > 0) ) d[dim-1] = val; }	///< Set sampling interval
    void setO(int dim, T val) { if( (dim <= MAXDIMS) && (dim > 0) ) o[dim-1] = val; }	///< Set origin

    // Print geometry functions
    void print(); ///< Print geometry information 

    // Clear geometry 
    void clear(); ///< Zeroes out the geometry

    // Compare geometry
    bool compare(std::shared_ptr<Geometry<T>> other);  ///< Returns 0 if non-zero entries in two geometries are equal

private:
    size_t n[MAXDIMS];	// Dimension sizes
    T d[MAXDIMS];	// Sampling interval
    T o[MAXDIMS];	// Origin 
};


// =============== 2D DATA GEOMETRY CLASS =============== //
/** The 2D data geometry class
 *
 */
template<typename T>
class Geometry2D: public Geometry<T> {
public:
    Geometry2D(size_t _n); 	///< Constructor
    virtual ~Geometry2D();       	///< Destructor
    /** Creates integer map from data coordinates and model geometry.
     * */
    void makeMap(std::shared_ptr<Geometry<T>> _geom); 

    // Get functions
    rockseis::Point2D<T> *getScoords() { return scoords; } ///< Get scoords
    rockseis::Point2D<T> *getGcoords() { return gcoords; } ///< Get gcoords
    rockseis::Point2D<int> *getSmap() { return smap; } ///< Get smap
    rockseis::Point2D<int> *getGmap() { return gmap; } ///< Get gmap

private:
    rockseis::Point2D<T> *scoords;             ///< Coordinates (x,z) of the source for the traces
    rockseis::Point2D<T> *gcoords;             ///< Coordinates (x,z) of the receivers for the traces
    rockseis::Point2D<int> *smap;              ///< Map of the coordinates (i,j) to the regular geometry grid of the source coordinates
    rockseis::Point2D<int> *gmap;              ///< Map of the coordinates (i,j) to the regular geometry grid of the receiver coordinates
};

// =============== 3D DATA GEOMETRY CLASS =============== //
/** The 3D data geometry class
 *
 */
template<typename T>
class Geometry3D: public Geometry<T> {
public:
    Geometry3D(int _n); 	///< Constructor
    virtual ~Geometry3D();       	///< Destructor
    /** Creates integer map from data coordinates and model geometry.
     * */
    void makeMap(std::shared_ptr<Geometry<T>> _geom);

    // Get functions
    rockseis::Point3D<T> *getScoords() { return scoords; } ///< Get scoords
    rockseis::Point3D<T> *getGcoords() { return gcoords; } ///< Get gcoords
    rockseis::Point3D<int> *getSmap() { return smap; } ///< Get smap
    rockseis::Point3D<int> *getGmap() { return gmap; } ///< Get gmap

private:
    rockseis::Point3D<T> *scoords;             ///< Coordinates (x,y,z) of the sources for the traces
    rockseis::Point3D<T> *gcoords;             ///< Coordinates (x,y,z) of the receivers for the traces
    rockseis::Point3D<int> *smap;              ///< Map of the coordinates (i,j,k) to the regular geometry grid of the source coordinates
    rockseis::Point3D<int> *gmap;              ///< Map of the coordinates (i,j,k) to the regular geometry grid of the receiver coordinates
};

}

#endif //GEOMETRY_H
