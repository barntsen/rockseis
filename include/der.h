#ifndef DER_H
#define DER_H

#include <stdio.h>
#include <stdint.h>
#include "utils.h"

namespace rockseis {
/* Der class */
template<typename T>
class Der {
public:
    Der();
    Der(const int _nx, const int _ny, const int _nz, const T _dx, const T _dy, const T _dz, const int _order);
    ~Der();
    /* Get functions */
    int getNx() { return nx; }		///< Get Nx
    int getNy() { return ny; }		///< Get Ny
    int getNz() { return nz; }		///< Get Nz
    int getOrder() { return order; } 	///< Get Order
    T getDx() { return dx; }		///< Get Dx
    T getDy() { return dy; }		///< Get Dy
    T getDz() { return dz; }		///< Get Dz
    T* getDf() { return df; }           ///< Get derivative function
    
    // Set functions
    void setNx(const int _nx) { nx = _nx; }	///< Set Nx
    void setNy(const int _ny) { ny = _ny; }	///< Set Ny
    void setNz(const int _nz) { nz = _nz; }	///< Set Nz
    void setDx(const T _dx) { dx = _dx; }	///< Set Dx
    void setDy(const T _dy) { dy = _dy; }	///< Set Dy
    void setDz(const T _dz) { dz = _dz; }	///< Set Dz
    void setOrder(const int _order);  ///< Reset the order of the FD approximation
    
    /* Derivative functions */
    //Forward
    void ddx_fw(T* f); ///< Forward derivative in x
    void ddy_fw(T* f); ///< Forward derivative in y
    void ddz_fw(T* f); ///< Forward derivative in z
    //Backward
    void ddx_bw(T* f); ///< Backward derivative in x
    void ddy_bw(T* f); ///< Backward derivative in y
    void ddz_bw(T* f); ///< Backward derivative in z
private:
    int nx;
    int ny;
    int nz;
    int order;
    T dx;
    T dy;
    T dz;
    T *coeffs;  // Vector of size order containing the coefficients of the FD stencil.
    T *df; // Array to store the derivative
};


}
#endif //DER_H
