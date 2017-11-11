#ifndef BSPL_H
#define BSPL_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "utils.h"

#define BSPL_OK 1
#define BSPL_ERR 0

namespace rockseis {

/* Bspl class */
template<typename T>
class Bspl {
public:
    Bspl();
    ~Bspl();
   
    /* Bspline evaluation function */
    void bspl(T *h, T *t, int k, T x, int l);
};

/* Bspl1D class */
template<typename T>
class Bspl1D: public Bspl<T> {
public:
    Bspl1D();
    Bspl1D(const int _nx, const T _dx, const T _dtx, const int _order);
    ~Bspl1D();

    T* getSpline() { return c; }           ///< Get spline coefficients
    T* getMod() { return mod; }           ///< Get evaluated spline model
    int getNc() { return nc; }           ///< Get spline coefficients

    /* Bspline evaluation function */
    void bisp();

private:
	T *c;		// B-spline coefficient array
	int kx;	 		// B-spline order 
	T *tx;    	 	// Knots in the x direction 
	T *x;    	 	// Evaluation array
	int mx;             // Length of evaluated model 
	int ntx; 		// Length of the Knot array
	int nc;                 // Number of coefficients
	T dtx;		// Size of the knot spacing
	T *mod;     // Model for evaulating spline
};

/* Bspl2D class */
template<typename T>
class Bspl2D: public Bspl<T> {
public:
    Bspl2D();
    Bspl2D(const int _nx, const int _nz, const T _dx, const T _dz, const T _dtx, const T _dtz, const int _kx, const int _kz);
    ~Bspl2D();

    T* getSpline() { return c; }           ///< Get spline coefficients
    T* getMod() { return mod; }           ///< Get evaluated spline model
    int getNc() { return nc; }           ///< Get spline coefficients

    /* Bspline evaluation function */
    void bisp();

private:
	T *c;		// B-spline coefficient array
	T **wx;		// Temporary array to store B-splines
	T **wz;	  	// Temporary array to store B-splines
	int kx;	 		// B-spline order 
	int kz;	 		// B-spline order 
	T *tx;    	 	// Knots in the x direction 
	T *tz;    	 	// Knots in the y direction 
	int *lx;    	 	// indexing(book keeper) array
	int *lz;    	 	// indexing(book keeper) array
	T *x;    	 	// Evaluation array
	T *z;    	 	// Evaluation array
	int mx, mz;             // Length of evaluated model 
	int ntx, ntz; 		// Length of the Knot array
	int nc;                 // Number of coefficients
	T dtx,dtz;		// Size of the knot spacing
	T *mod;     // Model for evaulating spline
};

/* Bspl3D class */
template<typename T>
class Bspl3D: public Bspl<T> {
public:
    Bspl3D();
    Bspl3D(const int _nx, const int _ny, const int _nz, const T _dx, const T _dy,const T _dz, const T _dtx, const T _dty, const T _dtz, const int _kx, const int _ky, const int _kz);
    ~Bspl3D();

    T* getSpline() { return c; }           ///< Get spline coefficients
    T* getMod() { return mod; }           ///< Get evaluated spline model
    int getNc() { return nc; }           ///< Get spline coefficients

    /* Bspline evaluation function */
    void bisp();

private:
	T *c;		// B-spline coefficient array
	T **wx;		// Temporary array to store B-splines
	T **wy;		// Temporary array to store B-splines
	T **wz;	  	// Temporary array to store B-splines
	int kx;	 		// B-spline order 
	int ky;	 		// B-spline order 
	int kz;	 		// B-spline order 
	T *tx;    	 	// Knots in the x direction 
	T *ty;    	 	// Knots in the y direction 
	T *tz;    	 	// Knots in the y direction 
	int *lx;    	 	// indexing(book keeper) array
	int *ly;    	 	// indeying(book keeper) array
	int *lz;    	 	// indexing(book keeper) array
	T *x;    	 	// Evaluation array
	T *y;    	 	// Evaluation array
	T *z;    	 	// Evaluation array
	int mx, my,mz;             // Length of evaluated model 
	int ntx,nty,ntz; 		// Length of the Knot array
	int nc;                 // Number of coefficients
	T dtx,dty,dtz;		// Size of the knot spacing
	T *mod;     // Model for evaulating spline
};



}

#endif //BSPL_H
