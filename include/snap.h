#ifndef SNAP_H
#define SNAP_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "data.h"
#include "file.h"
#include "waves.h"

#define SNAP_OK 1
#define SNAP_ERR 0

#define NPTR 3

namespace rockseis {

// =============== ABSTRACT SNAP CLASS =============== //
/** The abstract snap class
 *
 */
template<typename T>
class Snap {
public:
    Snap(); ///<Constructor
    virtual ~Snap();	///< Destructor

    // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }///< Set Nz
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set lpml
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setSnapnt(const int nt) { geometry->setN(4, nt); } ///< Set snapnt
    void setSnapdt(const T dt) { geometry->setD(4, dt); } ///< Set snapdt
    void setSnapot(const T ot) { geometry->setO(4, ot); } ///< Set snapot
    void setSnapit(const int it) { snapit = it; } ///< Set snapit
    void setEnddiff(const int it) { enddiff = it; } ///< Set enddiff
    void setSnapinc(const int inc); ///< Set snapinc
    void setData(T *ptr, int i) { if(i >= 0 && i < NPTR) data[i] = ptr; } ///< Set data pointer

    // Get functions
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
    T getDt() { return geometry->getD(4); }	///< Get Dt
    T getOx() { return geometry->getO(1); }	///< Get Ox
    T getOy() { return geometry->getO(2); }	///< Get Oy
    T getOz() { return geometry->getO(3); }	///< Get Oz
    T getOt() { return geometry->getO(4); }	///< Get Ot
    int getSnapnt() { return geometry->getN(4); }	///< Get Snap nt
    T getSnapdt() { return geometry->getD(4); }	///< Get Snap dt
    T getSnapot() { return geometry->getO(4); }	///< Get Snap ot
    int getSnapit() { return snapit; } ///< Get snapit
    int getEnddiff() { return enddiff; } ///< Get enddiff
    int getSnapinc() { return snapinc; } ///< Get snapinc
    T* getData(int i) { if(i >= 0 && i < 3) return data[i]; } ///< Get data pointer

    // Memory functions
    void allocSnap(); ///< Allocate data in snapshot
    void freeSnap(); ///< Free data in snapshot

    //File functions
    bool openSnap(std::string filename, char flag); ///< Open a snapshot for reading, writting or appending
    void closeSnap(); ///< Close snapshot file 
    void writeSnap(const int it); ///< Write snapshot
    void readSnap(const int it); ///< Read snapshot

   private:
    std::string filename; ///< filename
    std::shared_ptr<rockseis::File> Fp; ///< File handle
    std::shared_ptr<Geometry<T>> geometry; // regular geometry
    int lpml; // PML boundary size
    bool open; ///< flag to see if file is open
    rs_field field; ///< enum indicating which field to snap
    T *data[NPTR];
    bool allocated[NPTR]; ///< flag to see if data is allocated
    int enddiff;
    int snapit; 
    int snapinc;
};
}
#endif //SNAP_H
