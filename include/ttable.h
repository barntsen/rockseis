#ifndef TTABLE_H
#define TTABLE_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "data.h"
#include "file.h"
#include "rays.h"
#include "kdtree.h"

#define k2D(i,j) ((j)*nx +(i))
#define k2Di(i,j) ((j)*nx_i +(i))
#define k3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define k3Di(i,j,k) ((k)*nx_i*ny_i + (j)*nx_i + (i))

#define TTABLE_OK 1
#define TTABLE_ERR 0

#define MAXINTERP 5


namespace rockseis {

// =============== ABSTRACT TTABLE CLASS =============== //
/** The abstract ttable class
 *
 */
template<typename T>
class Ttable {
public:
    Ttable(); ///<Constructor
    Ttable(const int dim, const int ntable); ///<Constructor
    virtual ~Ttable();	///< Destructor

    // Set functions
    void setDim(const int _dim) { dim = _dim; } 	///< Set the dimension
    void setNx(const int _nx) { geometry->setN(1, _nx); }///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }///< Set Nz
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set lpml
    void setNtable(const int nt) { geometry->setN(4, nt); } ///< Set number of tables
    void setFilename(std::string name) { filename = name; } ///< Set filename
    void setRadius(T val) { radius = val; } ///< Set interpolation radius for position
    void setAllocated(bool val) { allocated = val; }

    // Get functions
    int getDim() { return dim; }		///< Get dimension
    int getNx() { return geometry->getN(1); }	///< Get Nx
    int getNx_pml() { return geometry->getN(1) + 2*lpml ; }	///< Get Nx padded
    int getNy() { return geometry->getN(2); }	///< Get Ny
    int getNy_pml() { return geometry->getN(2) + 2*lpml ; }	///< Get Ny padded
    int getNz() { return geometry->getN(3); }	///< Get Nz
    int getNz_pml() { return geometry->getN(3) + 2*lpml ; }	///< Get Nz padded
    int getLpml() { return lpml; }		///< Get lpml
    T getDx() { return geometry->getD(1); }	///< Get Dx
    T getDy() { return geometry->getD(2); }	///< Get Dy
    T getDz() { return geometry->getD(3); }	///< Get Dz
    T getOx() { return geometry->getO(1); }	///< Get Ox
    T getOy() { return geometry->getO(2); }	///< Get Oy
    T getOz() { return geometry->getO(3); }	///< Get Oz
    T getRadius() { return radius; } ///< Get interpolation radius for position
    kdtree* getKdtree() { return ptree; }
    std::shared_ptr<Geometry3D<T>> getGeom() { return geometry; } ///< Get data geometry
    int getNtable() { return geometry->getN(4); }	///< Get number of tables
    bool getAllocated() {return allocated;}
    std::string getFilename() { return filename; } ///< Get filename
    void makeGeom(const int ntable) { if(!geomset) { geometry = std::make_shared<Geometry3D<T>>(ntable);  geomset = true; }} ///< Make geometry


    // Floor function (to avoid calling std::floor which is unstable in Centos)
    int tabfloor( const T x ) { return x >= 0 ? (int) x : (int) x - 1; } ///< Rounds to the lowest integer

   private:
    std::string filename; ///< filename
    std::shared_ptr<Geometry3D<T>> geometry; ///< 3D geometry
    kdtree *ptree; ///< Tree structure to search for positions and neighbours
    int dim; 
    T radius;
    int lpml; 
    bool allocated; ///< flag to see if data is allocated
    bool geomset;
};

// =============== TTABLE 2D CLASS =============== //
/** The abstract ttable class
 *
 */
template<typename T>
class Ttable2D: public Ttable<T> {
public:
    Ttable2D(); ///<Constructor
    Ttable2D(std::shared_ptr<ModelEikonal2D<T>> model, int _ntable); ///<Constructor
    Ttable2D(std::string tablefile); ///<Constructor
    ~Ttable2D();	///< Destructor

    //File functions
    void createEmptyttable(); ///< Create empty dataset
    void fetchTtabledata(std::shared_ptr<RaysAcoustic2D<T>> rays, std::shared_ptr<Data2D<T>> source, const size_t number);
    void putTtabledata(std::shared_ptr<RaysAcoustic2D<T>> rays);
    void writeTtable(const size_t it); ///< Write travel time table
    void readTtable(const size_t it); ///< Read travel time table

    // Memory functions
    void allocTtable(); ///< Allocate data in travel time table
    void setData(T *ptr) {data = ptr; } ///< Set data pointer
    T* getData() {return data;} ///< Get data pointer
    T* getWrk() {return wrk;} ///< Get wrk pointer
    T* getDist() {return interpdist;} ///< Get interp dist array pointer
    size_t* getPch() {return interpch;} ///< Get pch pointer

    // Memory functions
    void freeTtable(); ///< Free data in travel time table

    //Insert source 
    void insertSource(std::shared_ptr<Data2D<T>> source, bool maptype, int traceno);

    //Interpolation function
    void interpTtable(std::shared_ptr<Ttable2D<T>> ttablei, T rad); ///< Interpolate to another ttable

   private:
    T interpoint2D[4];
    T interpoint1D[2];
    T *data;
    T *wrk;
    T *interpdist;
    size_t *interpch;

};

// =============== TTABLE 3D CLASS =============== //
/** The abstract ttable class
 *
 */
template<typename T>
class Ttable3D: public Ttable<T> {
public:
    Ttable3D(); ///<Constructor
    Ttable3D(std::shared_ptr<ModelEikonal3D<T>> model, int _ntable); ///<Constructor
    Ttable3D(std::string tablefile); ///<Constructor
    ~Ttable3D();	///< Destructor

    void createEmptyttable(); ///< Create empty dataset
    void fetchTtabledata(std::shared_ptr<RaysAcoustic3D<T>> rays, std::shared_ptr<Data3D<T>> source, const size_t number);
    void putTtabledata(std::shared_ptr<RaysAcoustic3D<T>> rays);
    void writeTtable(const size_t it); ///< Write travel time table
    void readTtable(const size_t it); ///< Read travel time table

    // Memory functions
    void allocTtable(); ///< Allocate data in travel time table
    void setData(T *ptr) {data = ptr; } ///< Set data pointer
    T* getData() {return data;} ///< Get data pointer
    T* getWrk() {return wrk;} ///< Get wrk pointer
    T* getDist() {return interpdist;} ///< Get interp dist array pointer
    size_t* getPch() {return interpch;} ///< Get pch pointer

    // Memory functions
    void freeTtable(); ///< Free data in travel time table

    //Insert source 
    void insertSource(std::shared_ptr<Data3D<T>> source, bool maptype, int traceno);

    //Interpolation function
    void interpTtable(std::shared_ptr<Ttable3D<T>> ttablei, T rad); ///< Interpolate to another ttable

   private:
    T interpoint3D[8];
    T interpoint2D[4];
    T interpoint1D[2];
    T *data;
    T *wrk;
    T *interpdist;
    size_t *interpch;
};



}
#endif //TTABLE_H
