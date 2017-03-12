#ifndef SNAP_H
#define SNAP_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "waves.h"
#include "utils.h"
#include "file.h"

namespace rockseis {

// =============== ABSTRACT SNAP CLASS =============== //
/** The abstract snap class
 *
 */
template<typename T>
class Snap {
public:
    Snap(); ///< Constructor
    virtual ~Snap();	///< Destructor
    
    // Get functions
    int getNx() { return geometry->getN(1); }		///< Get Nx
    int getNy() { return geometry->getN(2); }		///< Get Ny
    int getNz() { return geometry->getN(3); }		///< Get Nz
    int getLpml() { return lpml; }		///< Get Lpml
    int getNx_pml() { return geometry->getN(1) + 2 * lpml; }	///< Nx_pml = Nx + 2*lpml 
    int getNy_pml() { return geometry->getN(2) + 2 * lpml; }	///< Ny_pml = Ny + 2*lpml 
    int getNz_pml() { return geometry->getN(3) + 2 * lpml; }	///< Nz_pml = Nz + 2*lpml 
    T getDx() { return geometry->getD(1); }		///< Get Dx
    T getDy() { return geometry->getD(2); }		///< Get Dy
    T getDz() { return geometry->getD(3); }		///< Get Dz
    T getOx() { return geometry->getO(1); }		///< Get Ox
    T getOy() { return geometry->getO(2); }		///< Get Oy
    T getOz() { return geometry->getO(3); }		///< Get Oz
    int getEnd_diff() { return end_diff; }		///< Get End_diff
    int getSnapit() { return snapit; }		///< Return current snap number
    std::shared_ptr<Geometry<T>> getGeom() { return geometry; } ///< Get geometry
    
    // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }	///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }	///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }	///< Set Nz
    void setNt(const int _nt) { geometry->setN(4, _nt); }	///< Set Nt
    void setNt_mod(const int _nt) { nt_mod = _nt; }	///< Set Nt_mod
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set size of padding boundary for PML
    void setNx_pml() { nx_pml = geometry->getN(1) + 2 * lpml; }	///< Nx_pml = Nx + 2*lpml 
    void setNy_pml() { ny_pml = geometry->getN(2) + 2 * lpml; }	///< Ny_pml = Ny + 2*lpml 
    void setNz_pml() { nz_pml = geometry->getN(3) + 2 * lpml; }	///< Nz_pml = Nz + 2*lpml 
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setDt(const T _dt) { geometry->setD(4, _dt); }	///< Set Dt
    void setDt_mod(const T _dt) { dt_mod = _dt; }	///< Set Dt_mod
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setOt(const T _ot) { geometry->setO(4, _ot); }	///< Set Ot
    void setOt_mod(const T _ot) { ot_mod = _ot; }	///< Set Ot_mod
    void setEnd_diff(const int it) { end_diff = it; } ///< Set discrete time difference between end in snap and end in modelling
    void setSnapit(const int it) { snapit = it; } ///< Set current snap number
    
    // I/O functions
    virtual void readSnap() = 0;	///< Read snap (virtual)
    virtual void writeSnap() = 0;	///< write snap (virtual)

private:
    int lpml;
    std::shared_ptr<Geometry<T>> geometry; 
    int nx_pml;
    int ny_pml;
    int nz_pml;
    int nt_mod;
    T dt_mod;
    T ot_mod;
    int end_diff;
    int snapit; 
};

// =============== 2D ACOUSTIC SNAP CLASS =============== //
/** The 2D acoustic snap class
 *
 */
template<typename T>
class SnapAcoustic2D: public Snap<T> {
public:
    SnapAcoustic2D(std::shared_ptr<WavesAcoustic2D<T>> waves, const int snap_inc);	///< Constructor
    ~SnapAcoustic2D();	///< Destructor
    
    // I/O functions
    void readSnap();	 ///< Read a snap from file
    void writeSnap(); ///< Write a snap to file
    // Get functions
    std::string getAxfile() { return Axfile; } ///<
    std::string getAzfile() { return Azfile; } ///<
    std::string getPfile() { return Pfile; } ///<
    
private:
    std::string Axfile; ///< Filename to x-acceleration snap
    std::string Azfile; ///< Filename to z-acceleration snap
    std::string Pfile; ///< Filename to pressure snap

};

// =============== 3D ACOUSTIC SNAP CLASS =============== //
/** The 3D acoustic snap class
 *
 */
template<typename T>
class SnapAcoustic3D: public Snap<T> {
public:
    SnapAcoustic3D(std::shared_ptr<WavesAcoustic3D<T>> waves, const int snap_inc);	///< Constructor
    ~SnapAcoustic3D();	///< Destructor
    
    // I/O functions
    void readSnap();	///< Read a snap from file
    void writeSnap();	///< NOT IMPLEMENTED
    // Get functions
    std::string getAxfile() { return Axfile; } ///<
    std::string getAyfile() { return Ayfile; } ///<
    std::string getAzfile() { return Azfile; } ///<
    std::string getPfile() { return Pfile; } ///<

private:
    std::string Axfile; ///< Filename to x-acceleration snap
    std::string Ayfile; ///< Filename to y-acceleration snap
    std::string Azfile; ///< Filename to z-acceleration snap
    std::string Pfile; ///< Filename to pressure snap
};

// =============== 2D ELASTIC SNAP CLASS =============== //
/** The 2D elastic snap class
 *
 */
template<typename T>
class SnapElastic2D: public Snap<T> {
public:
    SnapElastic2D(std::shared_ptr<WavesElastic2D<T>> waves, const int snap_inc);	///< Constructor
    ~SnapElastic2D();	///< Destructor
    
    // I/O functions
    void readSnap();	///< Read a snap from file
    void writeSnap();	///< Write a snap to file

    // Get functions
    std::string getVxfile() { return Vxfile; } ///<
    std::string getVzfile() { return Vzfile; } ///<
    std::string getSxxfile() { return Sxxfile; } ///<
    std::string getSzzfile() { return Szzfile; } ///<
    std::string getSxzfile() { return Sxzfile; } ///<

private:
    std::string Sxxfile; ///< Filename to normal stress in x direction
    std::string Szzfile; ///< Filename to normal stress in z direction
    std::string Sxzfile; ///< Filename to shear stress 
    std::string Vxfile; ///< Filename to x-velocity snap
    std::string Vzfile; ///< Filename to z-velocity snap
};

// =============== 3D ELASTIC SNAP CLASS =============== //
/** The 3D elastic snap class
 *
 */
template<typename T>
class SnapElastic3D: public Snap<T> {
public:
    SnapElastic3D(std::shared_ptr<WavesElastic3D<T>> waves, const int snap_inc);	///< Constructor
    ~SnapElastic3D();	///< Destructor
    
    // I/O functions
    void readSnap();	///< Read a snap to file
    void writeSnap();	///< Write a snap to file

    // Get functions
    std::string getSxxfile() { return Sxxfile; }
    std::string getSyyfile() { return Syyfile; }
    std::string getSzzfile() { return Szzfile; }
    std::string getSxzfile() { return Sxzfile; }
    std::string getSxyfile() { return Sxyfile; }
    std::string getSyzfile() { return Syzfile; }
    std::string getVxfile() { return Vxfile; }
    std::string getVyfile() { return Vyfile; }
    std::string getVzfile() { return Vzfile; }

private:
    std::string Vxfile; ///< Filename to x-velocity snap
    std::string Vyfile; ///< Filename to y-velocity snap
    std::string Vzfile; ///< Filename to z-velocity snap
    std::string Sxxfile; ///< Filename to normal stress in x direction
    std::string Syyfile; ///< Filename to normal stress in y direction
    std::string Szzfile; ///< Filename to normal stress in z direction
    std::string Sxzfile; ///< Filename to shear stress 
    std::string Syzfile; ///< Filename to shear stress 
    std::string Sxyfile; ///< Filename to shear stress 
};

}
#endif //SNAP_H
