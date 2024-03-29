#ifndef GEOMETRY_H
#define GEOMETRY_H

// Include statements
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <cmath>
#include <iostream>
#include <memory>
#include "utils.h"

#define MAXDIMS 9
#define CTHRES 1e-2
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define CLOSE(x, y, d) ((ABS(x - y) <= (d)*ABS(y)) ? 1 : 0)

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
    size_t getNtot(); ///< Get total number of samples

    // Set functions
    void setN(int dim, size_t val); 	///< Set dimension size
    void setD(int dim, T val) { if( (dim <= MAXDIMS) && (dim > 0) ) d[dim-1] = val; }	///< Set sampling interval
    void setO(int dim, T val) { if( (dim <= MAXDIMS) && (dim > 0) ) o[dim-1] = val; }	///< Set origin

    // Floor function (to avoid calling std::floor which is unstable in Centos)
    int mapfloor( const T x ) { return x >= 0 ? (int) x : (int) x - 1; } ///< Rounds to the lowest integer

    // Ceil function (to avoid calling std::floor which is unstable in Centos)
    int mapceil( const T x ) { return x >= 0 ? (int) (x + 1) : (int) (x); } ///< Rounds to the highest integer

    // Print geometry functions
    void print(); ///< Print geometry information 

    // Clear geometry 
    void clear(); ///< Zeroes out the geometry

    // Compare geometry
    bool compare(std::shared_ptr<Geometry<T>> other);  ///< Returns 0 if non-zero entries in two geometries are equal

    bool getVerbose() { return verbose; } ///< If to print warnings 
    void setVerbose(bool val) { verbose = val; } ///< Switch to print warnings 

private:
    size_t n[MAXDIMS];	// Dimension sizes
    T d[MAXDIMS];	// Sampling interval
    T o[MAXDIMS];	// Origin 
    bool verbose;
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
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padhx, int padhy, T shiftx, T shifty); 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padhx, int padhy) {makeMap(_geom,map,padlx,padly,padhx,padhy, 0, 0);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, int padx, int pady) {makeMap(_geom,SMAP,padx,pady,padx,pady,0,0); makeMap(_geom,GMAP,padx,pady,padx,pady,0,0);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom) {makeMap(_geom,SMAP,0,0,0,0,0,0); makeMap(_geom,GMAP,0,0,0,0,0,0);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map) {makeMap(_geom,map,0,0,0,0,0,0);} 
    void copySmap2Gmap();
    void copyGmap2Smap();

    // Get functions
    rockseis::Point2D<T> *getScoords() { return scoords; } ///< Get scoords
    rockseis::Point2D<T> *getGcoords() { return gcoords; } ///< Get gcoords
    rockseis::Point2D<int> *getSmap() { return smap; } ///< Get smap
    rockseis::Point2D<int> *getGmap() { return gmap; } ///< Get gmap
    rockseis::Point2D<T> *getSshift() { return sshift; } ///< Get sshift
    rockseis::Point2D<T> *getGshift() { return gshift; } ///< Get gshift

private:
    rockseis::Point2D<T> *scoords;             ///< Coordinates (x,z) of the source for the traces
    rockseis::Point2D<T> *gcoords;             ///< Coordinates (x,z) of the receivers for the traces
    rockseis::Point2D<int> *smap;              ///< Map of the coordinates (i,j) to the regular geometry grid of the source coordinates
    rockseis::Point2D<int> *gmap;              ///< Map of the coordinates (i,j) to the regular geometry grid of the receiver coordinates
    rockseis::Point2D<T> *sshift;              ///< Value of error in source coordinate location between actual geometry and the regular geometry grid of the mapped coordinates
    rockseis::Point2D<T> *gshift;              ///< Value of error in receiver coordinate location between actual geometry and the regular geometry grid of the mapped coordinates
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
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padlz, int padhx, int padhy, int padhz, T shiftx, T shifty, T shiftz); 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padlz, int padhx, int padhy, int padhz) {makeMap(_geom,map,padlx,padly,padlz,padhx,padhy,padhz,0,0,0);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, int padx, int pady, int padz) {makeMap(_geom,SMAP,padx,pady,padz,padx,pady,padz); makeMap(_geom,GMAP,padx,pady,padz,padx,pady,padz);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom) {makeMap(_geom,SMAP,0,0,0,0,0,0); makeMap(_geom,GMAP,0,0,0,0,0,0);} 
    void makeMap(std::shared_ptr<Geometry<T>> _geom, bool map) {makeMap(_geom,map,0,0,0,0,0,0);} 
    void copySmap2Gmap();
    void copyGmap2Smap();

    // Get functions
    rockseis::Point3D<T> *getScoords() { return scoords; } ///< Get scoords
    rockseis::Point3D<T> *getGcoords() { return gcoords; } ///< Get gcoords
    rockseis::Point3D<int> *getSmap() { return smap; } ///< Get smap
    rockseis::Point3D<int> *getGmap() { return gmap; } ///< Get gmap
    rockseis::Point3D<T> *getSshift() { return sshift; } ///< Get sshift
    rockseis::Point3D<T> *getGshift() { return gshift; } ///< Get gshift

private:
    rockseis::Point3D<T> *scoords;             ///< Coordinates (x,y,z) of the sources for the traces
    rockseis::Point3D<T> *gcoords;             ///< Coordinates (x,y,z) of the receivers for the traces
    rockseis::Point3D<int> *smap;              ///< Map of the coordinates (i,j,k) to the regular geometry grid of the source coordinates
    rockseis::Point3D<int> *gmap;              ///< Map of the coordinates (i,j,k) to the regular geometry grid of the receiver coordinates
    rockseis::Point3D<T> *sshift;              ///< Value of error in source coordinate location between actual geometry and the regular geometry grid of the mapped coordinates
    rockseis::Point3D<T> *gshift;              ///< Value of error in receiver coordinate location between actual geometry and the regular geometry grid of the mapped coordinates
};

}

#endif //GEOMETRY_H
