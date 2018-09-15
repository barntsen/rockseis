#ifndef RAYS_H
#define RAYS_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "model.h"
#include "geometry.h"
#include "utils.h"
#include "data.h"
#include "file.h"
#include "snap.h"
#include "interp.h"
#include <algorithm>

#define RAYS_OK 1;
#define RAYS_ERR 0;

#define LANC_SIZE 3
#define LANC(x,a) (this->sinc(x)*this->sinc((x)/a))
#define TTNORM 1e-10
#define EPS_ADJ 1.0e-12
#define SQ(x) ((x)*(x))
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#define TMAX 99.0
#define MAXITER 30

namespace rockseis {

/** The abstract waves class
 *
 */
template<typename T>
class Rays {
public:
    Rays();		///< Constructor
    Rays(const int _dim, const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz);	///< Constructor
    virtual ~Rays();	///< Destructor
    
    // Get functions
    int getDim() { return dim; }		///< Get dimension
    int getNx() { return geometry->getN(1); }	///< Get Nx
    int getNx_pml() { return geometry->getN(1) + 2*lpml ; }	///< Get Nx padded
    int getNy() { return geometry->getN(2); }	///< Get Ny
    int getNy_pml() { return geometry->getN(2) + 2*lpml ; }	///< Get Ny padded
    int getNz() { return geometry->getN(3); }	///< Get Nz
    int getNz_pml() { return geometry->getN(3) + 2*lpml ; }	///< Get Nz padded
    T getDx() { return geometry->getD(1); }	///< Get Dx
    T getDy() { return geometry->getD(2); }	///< Get Dy
    T getDz() { return geometry->getD(3); }	///< Get Dz
    T getOx() { return geometry->getO(1); }	///< Get Ox
    T getOy() { return geometry->getO(2); }	///< Get Oy
    T getOz() { return geometry->getO(3); }	///< Get Oz
    int getLpml() { return lpml; }		///< Get lpml

    //Interpolation function
    T sinc(T x) { if(x == 0) return 1; else return std::sin(PI*x)/(PI*x); } 

        // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }///< Set Nz
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setDim(const int _dim) { dim = _dim; }		///< Set the dimension
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set lpml

private:
    int dim;
    std::shared_ptr<Geometry<T>> geometry; // regular geometry
    int lpml;

};

/** The 2D Acoustic RAYS class
 *
 */
template<typename T>
class RaysAcoustic2D: public Rays<T> {
public:
    RaysAcoustic2D();					///< Constructor
    RaysAcoustic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz);	///< Constructor
    RaysAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model);	///< Constructor
    ~RaysAcoustic2D();					///< Destructor
    
    // Eikonal solver
    T norm1(T *TT, T *TTold);
    void sweep(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy); ///< Compute first arrival traveltimes
    void sweep_adj(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy);  ///< Compute adjoint eikonal equation
    void solve();
    void solve_adj();

    // Get functions
    T *getTT() { return TT; } 
    T *getLam() { return lam; } 
    bool *getRecmask() { return recmask; } 

    // Insert source functions
    void insertSource(std::shared_ptr<Data2D<T>> source, bool maptype); ///< Insert source position
    void insertResiduals(std::shared_ptr<Data2D<T>> source, bool maptype); ///< Insert source position
    void createRecmask(std::shared_ptr<Data2D<T>> source, bool maptype); ///< Create a mask for adjoint eikonal computation

    // Record data at receivers functions
    void recordData(std::shared_ptr<Data2D<T>> data, bool maptype); ///< Record traveltime at receivers 

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    T *TT; // Traveltime
    T *lam; // Adjoint state traveltime
    bool *recmask; // Boolean array indicating where there are receivers
};

/** The 3D Acoustic RAYS class
 *
 */
template<typename T>
class RaysAcoustic3D: public Rays<T> {
public:
    RaysAcoustic3D();					///< Constructor
    RaysAcoustic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz);	///< Constructor
    RaysAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model);	///< Constructor
    ~RaysAcoustic3D();					///< Destructor
    
    // Eikonal solver
    T norm1(T *TT, T *TTold);
    void sweep(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int nz1, int nz2, int ndz); ///< Compute first arrival traveltimes
    void sweep_adj(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int nz1, int nz2, int ndz); ///< Compute first arrival traveltimes
    void solve();
    void solve_adj();

    // Get functions
    T *getTT() { return TT; } 
    T *getLam() { return lam; } 
    bool *getRecmask() { return recmask; } 

    // Insert source functions
    void insertSource(std::shared_ptr<Data3D<T>> source, bool maptype); ///< Insert source position
    void insertResiduals(std::shared_ptr<Data3D<T>> source, bool maptype); ///< Insert source position
    void createRecmask(std::shared_ptr<Data3D<T>> source, bool maptype); ///< Create a mask for adjoint eikonal computation

    // Record data at receivers functions
    void recordData(std::shared_ptr<Data3D<T>> data, bool maptype); ///< Record traveltime at receivers 

private:
    std::shared_ptr<ModelAcoustic3D<T>> model;
    T *TT; // Traveltime
    T *lam; // Adjoint state traveltime
    bool *recmask; // Boolean array indicating where there are receivers
};


}
#endif //RAYS_H
