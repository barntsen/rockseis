/* 
   Tools for evaluating functions from B-spline coefficients.

   Author: Wiktor Weibull, wiktor.weibull@ntnu.no, NTNU, 2013
*/

/*
  Copyright (C) 2013 Norwegian University of Science and Technology 
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* Including libraries */
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct bspl_spline1 *bspl_1dspline;
/*^*/

struct bspl_spline1{
	float *c;		// B-spline coefficient array
	int kx;	 		// B-spline order 
	float *tx;    	 	// Knots in the x direction 
	float *x;    	 	// Evaluation array
	int mx;             // Length of evaluated model 
	int ntx; 		// Length of the Knot array
	int nc;                 // Number of coefficients
	float dtx;		// Size of the knot spacing
};
/*^*/

typedef struct bspl_spline2 *bspl_2dspline;
/*^*/

struct bspl_spline2{
	float *c;		// B-spline coefficient array
	float **wx;		// Temporary array to store B-splines
	float **wz;	  	// Temporary array to store B-splines
	int kx;	 		// B-spline order 
	int kz;	 		// B-spline order 
	float *tx;    	 	// Knots in the x direction 
	float *tz;    	 	// Knots in the y direction 
	int *lx;    	 	// indexing(book keeper) array
	int *lz;    	 	// indexing(book keeper) array
	float *x;    	 	// Evaluation array
	float *z;    	 	// Evaluation array
	int mx, mz;             // Length of evaluated model 
	int ntx, ntz; 		// Length of the Knot array
	int nc;                 // Number of coefficients
	float dtx,dtz;		// Size of the knot spacing
};
/*^*/

typedef struct bspl_spline3 *bspl_3dspline;
/*^*/

struct bspl_spline3{
	float *c;		// B-spline coefficient array
	float **wx;		// Temporary array to store B-splines
	float **wy;	  	// Temporary array to store B-splines
	float **wz;	  	// Temporary array to store B-splines
	int kx;	 		// B-spline order 
	int ky;	 		// B-spline order 
	int kz;	 		// B-spline order 
	float *tx;    	 	// Knots in the x direction 
	float *ty;    	 	// Knots in the x direction 
	float *tz;    	 	// Knots in the y direction 
	int *lx;    	 	// indexing(book keeper) array
	int *ly;    	 	// indexing(book keeper) array
	int *lz;    	 	// indexing(book keeper) array
	float *x;    	 	// Evaluation array
	float *y;    	 	// Evaluation array
	float *z;    	 	// Evaluation array
	int mx, my, mz;          // Length of evaluated model 
	int ntx, nty, ntz; 	// Length of the Knot array
	int nc;                 // Number of coefficients
	float dtx, dty, dtz;	// Size of the knot spacing
};
/*^*/


bspl_1dspline bspl_1dspline_init(int nx, float dx, double dtx, int kx); 
/*<Initializing the 1D B-spline struct.>*/

bspl_2dspline bspl_2dspline_init(int nx,int nz, float dx, float dz, double dtx, double dtz, int kx, int kz); 
/*<Initializing the 2D B-spline struct.>*/

bspl_3dspline bspl_3dspline_init(int nx, int ny, int nz, float dx, float dy, float dz, double dtx, double dty, double dtz, int kx, int ky, int kz); 
/*<Initializing the 2D B-spline struct.>*/

void bspl_1dspline_free(bspl_1dspline spline);
/*<Free bspl_1dspline structure from memory.>*/

void bspl_2dspline_free(bspl_2dspline spline);
/*<Free bspl_2dspline structure from memory.>*/

void bspl_3dspline_free(bspl_3dspline spline);
/*<Free bspl_3dspline structure from memory.>*/

void bspl(float *h, float *t, int k, float x, int l);
/*<evaluates the (k+1) non-zero b-splines of degree k at t(l) <= x < t(l+1) using the stable recurrence relation of de Boor and Cox.>*/

void bspl_bisp1d(float *mod, bspl_1dspline spline);
/*<Evaluates the spline coefficients in spline->c over z.>*/

void bspl_bisp2d(float *mod, bspl_2dspline spline);
/*<Evaluates the spline coefficients in spline->c over z.>*/

void bspl_bisp3d(float *mod, bspl_3dspline spline);
/*<Evaluates the spline coefficients in spline->c over z.>*/

