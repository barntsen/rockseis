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

#define NPTR 9

namespace rockseis {

// =============== ABSTRACT SNAP CLASS =============== //
/** The abstract snap class
 *
 */
template<typename T>
class Snapshot {
public:
    Snapshot(); ///<Constructor
    virtual ~Snapshot();	///< Destructor

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
    void setSnapinc(const int inc) { snapinc = inc; } ///< Set snapinc
    void setData(T *ptr, int i) { if(i >= 0 && i < NPTR) data[i] = ptr; } ///< Set data pointer
    void setDim(const int _dim) { dim = _dim; }	///< Set lpml

    // Get functions
    int getNx() { return geometry->getN(1); }	///< Get Nx
    int getNx_pml() { return geometry->getN(1) + 2*lpml ; }	///< Get Nx padded
    int getNy() { return geometry->getN(2); }	///< Get Ny
    int getNy_pml() { return geometry->getN(2) + 2*lpml ; }	///< Get Ny padded
    int getNz() { return geometry->getN(3); }	///< Get Nz
    int getNz_pml() { return geometry->getN(3) + 2*lpml ; }	///< Get Nz padded
    int getLpml() { return lpml; }		///< Get lpml
    int getDim() { return dim; }		///< Get dim
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
    bool getOpen() { return open; }
    bool getAllocated(int i) { if(i >= 0 && i < 3) return allocated[i]; else {rs_error("Eror in getAllocated, index is out of bounds."); return 0;} }
    T* getData(int i) { if(i >= 0 && i < 3) return data[i]; else return nullptr;} ///< Get data pointer
    std::shared_ptr<rockseis::File> getFp() { return Fp; } ///< Get pointer to file handle

    // Memory functions
    void allocSnap(const int i); ///< Allocate data in snapshot
    void freeSnaps(); ///< Free data in snapshot

    //File functions
    bool openSnap(std::string filename, char flag); ///< Open a snapshot for reading, writting or appending
    bool openEdge(std::string filename, char flag); ///< Open a edge snapshot for reading, writting or appending
    void closeSnap(); ///< Close snapshot file 
    void removeSnap(); ///< Delete snapshot file

   private:
    std::string filename; ///< filename
    std::shared_ptr<rockseis::File> Fp; ///< File handle
    std::shared_ptr<Geometry<T>> geometry; // regular geometry
    int dim;
    int lpml; // PML boundary size
    bool open; ///< flag to see if file is open
    rs_field field; ///< enum indicating which field to snap
    T *data[NPTR];
    bool allocated[NPTR]; ///< flag to see if data is allocated
    int enddiff;
    int snapit; 
    int snapinc;
};


template<typename T>
class Snapshot2D: public Snapshot<T> {
public:
    Snapshot2D(std::shared_ptr<ModelEikonal2D<T>> model, int snapinc); 	///< Constructor
    Snapshot2D(std::shared_ptr<WavesAcoustic2D<T>> waves, int snapinc); 	///< Constructor
    Snapshot2D(std::shared_ptr<WavesElastic2D<T>> waves, int snapinc); 	///< Constructor
    Snapshot2D(std::shared_ptr<WavesElastic2D_DS<T>> waves, int snapinc); 	///< Constructor
    ~Snapshot2D();       	///< Destructor

    void writeSnap(const int it); ///< Write snapshot
    void readSnap(const int it); ///< Read snapshot

    void writeEdge(const int it); ///< Write edges
    void readEdge(const int it); ///< Read edges
};

template<typename T>
class Snapshot3D: public Snapshot<T> {
public:
    Snapshot3D(std::shared_ptr<ModelEikonal3D<T>> model, int snapinc); 	///< Constructor
    Snapshot3D(std::shared_ptr<WavesAcoustic3D<T>> waves, int snapinc); 	///< Constructor
    Snapshot3D(std::shared_ptr<WavesElastic3D<T>> waves, int snapinc); 	///< Constructor
    Snapshot3D(std::shared_ptr<WavesElastic3D_DS<T>> waves, int snapinc); 	///< Constructor
    ~Snapshot3D();       	///< Destructor

    void writeSnap(const int it); ///< Write snapshot
    void readSnap(const int it); ///< Read snapshot

    void writeEdge(const int it); ///< Write edges
    void readEdge(const int it); ///< Read edges
};

}
#endif //SNAP_H
