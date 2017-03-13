#ifndef SNAP_H
#define SNAP_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "model.h"
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
    int getNt() { return geometry->getN(4); }		///< Get Nt
    int getNt_mod() { return nt_mod; }		///< Get modelling nt
    int getLpml() { return lpml; }		///< Get Lpml
    int getNx_pml() { return geometry->getN(1) + 2 * lpml; }	///< Nx_pml = Nx + 2*lpml 
    int getNy_pml() { return geometry->getN(2) + 2 * lpml; }	///< Ny_pml = Ny + 2*lpml 
    int getNz_pml() { return geometry->getN(3) + 2 * lpml; }	///< Nz_pml = Nz + 2*lpml 
    T getDx() { return geometry->getD(1); }		///< Get Dx
    T getDy() { return geometry->getD(2); }		///< Get Dy
    T getDz() { return geometry->getD(3); }		///< Get Dz
    T getDt() { return geometry->getD(4); }		///< Get Dt
    T getDt_mod() { return dt_mod; }		///< Get modelling dt
    T getOx() { return geometry->getO(1); }		///< Get Ox
    T getOy() { return geometry->getO(2); }		///< Get Oy
    T getOz() { return geometry->getO(3); }		///< Get Oz
    T getOt() { return geometry->getO(4); }		///< Get Oz
    T getOt_mod() { return ot_mod; }		///< Get modelling ot
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
    SnapAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, const int _nt_mod, T _dt_mod, T _ot_mod, const int snap_inc);	///< Constructor
    ~SnapAcoustic2D();	///< Destructor
    
    // I/O functions
    void readSnap();	 ///< Read a snap from file
    void writeSnap(); ///< Write a snap to file
    // Get functions
    std::string getAxfile() { return Axfile; } ///< Get Ax snap file
    std::string getAzfile() { return Azfile; } ///< Get Az snap file
    std::string getPfile() { return Pfile; } ///< Get P snap file
    bool getP() { return P; } ///< Get P flag
    bool getAx() { return Ax; } ///< Get Ax flag
    bool getAz() { return Az; } ///< Get Az flag

    // Put functions
    void putAzfile(std::string filename); ///< Set Az snap file ond open file for writting
    void putAxfile(std::string filename); ///< Set Ax snap file ond open file for writting
    void putPfile(std::string filename); ///< Set P snap file ond open file for writting
    
private:
    std::string Axfile; 
    std::string Azfile; 
    std::string Pfile; 
    std::shared_ptr<rockseis::File> FAx;
    std::shared_ptr<rockseis::File> FAz;
    std::shared_ptr<rockseis::File> FP;
    bool Ax, Az, P;
};

// =============== 3D ACOUSTIC SNAP CLASS =============== //
/** The 3D acoustic snap class
 *
 */
template<typename T>
class SnapAcoustic3D: public Snap<T> {
public:
    SnapAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, const int snap_inc);	///< Constructor
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
    std::string Axfile; 
    std::string Ayfile; 
    std::string Azfile; 
    std::string Pfile; 
    bool Ax, Ay, Az, P; 
};

// =============== 2D ELASTIC SNAP CLASS =============== //
/** The 2D elastic snap class
 *
 */
template<typename T>
class SnapElastic2D: public Snap<T> {
public:
    SnapElastic2D(std::shared_ptr<ModelElastic2D<T>> model, const int snap_inc);	///< Constructor
    ~SnapElastic2D();	///< Destructor
    
    // I/O functions
    void readSnap();	///< Read a snap from file
    void writeSnap();	///< Write a snap to file

    // Get functions
    std::string getVxfile() { return Vxfile; } ///< Get Vx snap file
    std::string getVzfile() { return Vzfile; } ///< Get Vx snap file
    std::string getSxxfile() { return Sxxfile; } ///< Get Sxx snap file
    std::string getSzzfile() { return Szzfile; } ///< Get Szz snap file
    std::string getSxzfile() { return Sxzfile; } ///< Get Sxz snap file

private:
    std::string Sxxfile; 
    std::string Szzfile; 
    std::string Sxzfile; 
    std::string Vxfile; 
    std::string Vzfile; 
    bool Vx, Vz, Sxx, Szz, Sxz; 
};

// =============== 3D ELASTIC SNAP CLASS =============== //
/** The 3D elastic snap class
 *
 */
template<typename T>
class SnapElastic3D: public Snap<T> {
public:
    SnapElastic3D(std::shared_ptr<ModelElastic3D<T>> model, const int snap_inc);	///< Constructor
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
    std::string Vxfile; 
    std::string Vyfile; 
    std::string Vzfile; 
    std::string Sxxfile;
    std::string Syyfile;
    std::string Szzfile;
    std::string Sxzfile;
    std::string Syzfile;
    std::string Sxyfile;
    bool Vx, Vy, Vz, Sxx, Syy, Szz, Sxz, Syz, Sxy; 
};

}
#endif //SNAP_H
