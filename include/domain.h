#ifndef DOMAIN_H
#define DOMAIN_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"

#define DOMAIN_OK 1
#define DOMAIN_ERR 0

namespace rockseis {

// =============== ABSTRACT DOMAIN CLASS =============== //
/** The abstract snap class
 *
 */
template<typename T>
class Domain {
public:
    Domain(); ///<Constructor
    virtual ~Domain();  ///< Destructor

    // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }///< Set Nz
    void setLpml(const int _lpml) { lpml = _lpml; }     ///< Set lpml
    void setDim(const int _dim) { dim = _dim; } ///< Set dim
    void setDx(const T _dx) { geometry->setD(1, _dx); } ///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); } ///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); } ///< Set Dz
    void setOx(const T _ox) { geometry->setO(1, _ox); } ///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); } ///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); } ///< Set Oz
    void setLow(const int val) { low = val; }   ///< Set rank of low side 
    void setHigh(const int val) { high = val; } ///< Set rank of high side 
    void setPadl(const int _padl) { padl = _padl; }     ///< Set padl
    void setPadh(const int _padh) { padh = _padh; }     ///< Set padh
    void setNd(const int _nd) { nd=_nd; }///< Set number of domains
    void setD(const int _d) { d=_d; }///< Set domain number

    // Get functions
    int getNx() { return geometry->getN(1); }   ///< Get Nx
    int getNx_pml() { return geometry->getN(1) + 2*lpml ; }     ///< Get Nx padded
    int getNy() { return geometry->getN(2); }   ///< Get Ny
    int getNy_pml() { return geometry->getN(2) + 2*lpml ; }     ///< Get Ny padded
    int getNz() { return geometry->getN(3); }   ///< Get Nz
    int getNz_pml() { return geometry->getN(3) + 2*lpml ; }     ///< Get Nz padded
    int getLpml() { return lpml; }              ///< Get lpml
    int getDim() { return dim; }                ///< Get dim
    int getLow() { return low; }                ///< Get low rank 
    int getHigh() { return high; }              ///< Get high rank 
    int getPadl() { return padl; }              ///< Get padl
    int getPadh() { return padh; }              ///< Get padh
    int getNd() { return nd; } ///< GetNd
    int getD() { return d; } /// Get d
    T getDx() { return geometry->getD(1); }     ///< Get Dx
    T getDy() { return geometry->getD(2); }     ///< Get Dy
    T getDz() { return geometry->getD(3); }     ///< Get Dz
    T getOx() { return geometry->getO(1); }     ///< Get Ox
    T getOy() { return geometry->getO(2); }     ///< Get Oy
    T getOz() { return geometry->getO(3); }     ///< Get Oz

   private:
    std::shared_ptr<Geometry<T>> geometry; // regular geometry
    int lpml; // PML boundary size
    int low,high;
    int dim;
    int padl, padh;
    int nd;
    int d;
};

}
#endif //DOMAIN_H

