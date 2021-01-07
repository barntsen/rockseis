#ifndef DOMAIN_H
#define DOMAIN_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
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
    void setLpad(const int _lpad) { lpad = _lpad; }     ///< Set lpad
    void setDim(const int _dim) { dim = _dim; } ///< Set dim
    void setLow(const int val) { low = val; }   ///< Set rank of low side 
    void setHigh(const int val) { high = val; } ///< Set rank of high side 
    void setPadl(const int _padl) { padl = _padl; }     ///< Set padl
    void setPadh(const int _padh) { padh = _padh; }     ///< Set padh
    void setNx_orig(const int val) { nx_orig = val; }     ///< Set nx_orig
    void setNy_orig(const int val) { ny_orig = val; }     ///< Set ny_orig
    void setNz_orig(const int val) { nz_orig = val; }     ///< Set nz_orig
    void setNx_pad(const int val) { nx_pad = val; }     ///< Set nx_pad
    void setNy_pad(const int val) { ny_pad = val; }     ///< Set ny_pad
    void setNz_pad(const int val) { nz_pad = val; }     ///< Set nz_pad
    void setNd(const int _nd) { nd=_nd; }///< Set number of domains
    void setD(const int _d) { d=_d; }///< Set domain number

    void setIx0(const int _ix0) { ix0 = _ix0; }     ///< Set ix0 (origin in global model)
    void setIy0(const int _iy0) { iy0 = _iy0; }     ///< Set iy0 (origin in global model)
    void setIz0(const int _iz0) { iz0 = _iz0; }     ///< Set iz0 (origin in global model)

    // Get functions
    int getNx_orig() { return nx_orig; }     ///< Get original model Nx 
    int getNy_orig() { return ny_orig ; }     ///< Get original model
    int getNz_orig() { return nz_orig ; }     ///< Get original model
    int getNx_pad() { return nx_pad; }     ///< Get Nx padded
    int getNy_pad() { return ny_pad ; }     ///< Get Ny padded
    int getNz_pad() { return nz_pad ; }     ///< Get Nz padded
    int getLpad() { return lpad; }              ///< Get lpad
    int getDim() { return dim; }                ///< Get dim
    int getLow() { return low; }                ///< Get low rank 
    int getHigh() { return high; }              ///< Get high rank 
    int getPadl() { return padl; }              ///< Get padl
    int getPadh() { return padh; }              ///< Get padh
    int getNd() { return nd; } ///< GetNd
    int getD() { return d; } /// Get d
    int getIx0() { return ix0; } /// Get ix0
    int getIy0() { return iy0; } /// Get iy0
    int getIz0() { return iz0; } /// Get iz0
    bool getStatus() { return status; }
    void setupDomain(const int nx, const int ny, const int nz, const int d, const int nd, const int order);

   private:
    int nx_orig, ny_orig, nz_orig;
    int nx_pad, ny_pad, nz_pad;
    int lpad; // PML boundary size
    int low,high;
    int dim;
    int padl, padh;
    int nd;
    int d;
    int ix0, iy0, iz0; // Global indexes for the origin of domain 
    bool status; // 1 for on, 0 for off. 
};

}
#endif //DOMAIN_H

