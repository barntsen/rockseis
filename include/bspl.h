#ifndef BSPL_H
#define BSPL_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "utils.h"
extern "C" {
#include "bsplc.h"
}

#define BSPL_OK 1
#define BSPL_ERR 0

namespace rockseis {

/* Bspl class */
class Bspl {
public:
    Bspl();
    Bspl(const int _nx, const int _ny, const int _nz, const float _dx, const float _dy, const float _dz, const double _dtx, const double _dty, const double _dtz, const int _order, const int _dim);
    ~Bspl();
    /* Get functions */
    int getNx() { return nx; }		///< Get Nx
    int getNy() { return ny; }		///< Get Ny
    int getNz() { return nz; }		///< Get Nz
    int getNc() { return nc; }		///< Get Nc
    int getOrder() { return order; } 	///< Get Order
    float getDx() { return dx; }		///< Get Dx
    float getDy() { return dy; }		///< Get Dy
    float getDz() { return dz; }		///< Get Dz
    double getDtx() { return dtx; }		///< Get Dtx
    double getDty() { return dty; }		///< Get Dty
    double getDtz() { return dtz; }		///< Get Dtz
    float* getMod() { return mod; }           ///< Get evaluatedt spline
    float* getSpline();           ///< Get spline coefficients
    
    // Set functions
    void setNx(const int _nx) { nx = _nx; }	///< Set Nx
    void setNy(const int _ny) { ny = _ny; }	///< Set Ny
    void setNz(const int _nz) { nz = _nz; }	///< Set Nz
    void setDx(const float _dx) { dx = _dx; }	///< Set Dx
    void setDy(const float _dy) { dy = _dy; }	///< Set Dy
    void setDz(const float _dz) { dz = _dz; }	///< Set Dz
    
    /* Bspline evaluation function */
    void bisp();

private:
    int nx;
    int ny;
    int nz;
    int order;
    float dx;
    float dy;
    float dz;
    double dtx;
    double dty;
    double dtz;
    int nc;
    int dim;
    bspl_1dspline spline1;
    bspl_2dspline spline2;
    bspl_3dspline spline3;
    float *mod; // Array to store the evaluated spline
    bool allocated;
};

}
#endif //BSPL_H
