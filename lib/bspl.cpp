#include "bspl.h"

namespace rockseis {
/* Constructors*/
Bspl::Bspl()
{
    rs_error("Bspl::Bspl(): Default constructor is not implemented.");
}

Bspl::Bspl(const int _nx, const int _ny, const int _nz, const float _dx, const float _dy, const float _dz, const double _dtx, const double _dty, const double _dtz, const int _order, const int _dim)
{
nx = _nx;
ny = _ny;
nz = _nz;

dx = _dx;
dy = _dy;
dz = _dz;

dtx = _dtx;
dty = _dty;
dtz = _dtz;
dim = _dim;

    switch (dim){
        case 1:
            spline1 = bspl_1dspline_init(_nx, _dx, _dtx, _order);
            nc = spline1->nc;
            break;
        case 2:
            spline2 = bspl_2dspline_init(_nx, _nz, _dx, _dz, _dtx, _dtz, _order, _order);
            nc = spline2->nc;
            break;
        case 3:
            spline3 = bspl_3dspline_init(_nx, _ny, _nz, _dx, _dy, _dz, _dtx, _dty, _dtz, _order, _order, _order);
            nc = spline3->nc;
            break;
        default:
            rs_error("Bspl::Bspl(): Dimension must be 1, 2 or 3");
            break;
    }
    allocated = false;
}

/* Destructor*/
Bspl::~Bspl(){
    /* Free allocated variables */
    switch (dim){
        case 1:
            bspl_1dspline_free(spline1);
            break;
        case 2:
            bspl_2dspline_free(spline2);
            break;
        case 3:
            bspl_3dspline_free(spline3);
            break;
        default:
            rs_error("Bspl::Bspl(): Dimension must be 1, 2 or 3");
            break;
    }
    if(this->allocated){
        free(mod);
        allocated = false;
    }
}

void Bspl::bisp()
{
    switch (dim){
        case 1:
            if(!allocated){
                mod = (float *) calloc(nx, sizeof(float));
                allocated = true;
            }
			bspl_bisp1d(mod, spline1); // Evaluate spline for this coefficient
            break;
        case 2:
            if(!allocated){
                mod = (float *) calloc(nx*nz, sizeof(float));
                allocated = true;
            }
            bspl_bisp2d(mod, spline2); // Evaluate spline for this coefficient
            break;
        case 3:
            if(!allocated){
                mod = (float *) calloc(nx*ny*nz, sizeof(float));
                allocated = true;
            }
			bspl_bisp3d(mod, spline3); // Evaluate spline for this coefficient
            break;
        default:
            rs_error("Bspl::Bspl(): Dimension must be 1, 2 or 3");
            break;
    }
}

float * Bspl::getSpline(){
    switch (dim){
        case 1:
            return spline1->c;
            break;
        case 2:
            return spline2->c;
            break;
        case 3:
            return spline3->c;
            break;
        default:
            rs_error("Bspl::Bspl(): Dimension must be 1, 2 or 3");
            break;
    }
    return NULL;
}

}
