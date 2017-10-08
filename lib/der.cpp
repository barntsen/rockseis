#include "der.h"

#define w11  1.0029f
/*^*/

#define w21  1.188401e+0
#define w22 -7.046382e-2
/*^*/

#define w31  1.216624e+0
#define w32 -9.197724e-2
#define w33  1.300041e-2
/*^*/

#define w41  1.230862e+0
#define w42 -1.034123e-1
#define w43  2.011671e-2
#define w44 -3.245760e-3
/*^*/

#define w51  1.239407e+0
#define w52 -1.105315e-1
#define w53  2.496329e-2
#define w54 -5.804879e-3
#define w55  9.358680e-4
/*^*/

#define w61  1.245095e+0
#define w62 -1.153979e-1
#define w63  2.848442e-2
#define w64 -7.899473e-3
#define w65  1.898222e-3
#define w66 -2.935304e-4
/*^*/

#define w71  1.249150e+0
#define w72 -1.189375e-1
#define w73  3.116091e-2
#define w74 -9.623004e-3
#define w75  2.812069e-3
#define w76 -6.655106e-4
#define w77  9.728191e-5
/*^*/

#define w81  1.252186e+0
#define w82 -1.216289e-1
#define w83  3.326513e-2
#define w84 -1.105777e-2
#define w85  3.647865e-3
#define w86 -1.066144e-3
#define w87  2.437491e-4
#define w88 -3.351959e-5
/*^*/

#define CARRAYSIZE 36 // Array of coefficients
#define ind(i,j,k) ((k)*ny*nx + (j)*nx + (i))
#define cind(i) (((i)*(i+1))/2)

namespace rockseis {
/* Constructors*/
template <typename T>
Der<T>::Der()
{
    nx = 1;
    ny = 1;
    nz = 1;
    dx = 1.;
    dy = 1.;
    dz = 1.;
    order = 1;
    coeffs = (T *) calloc(CARRAYSIZE,sizeof(T));
    coeffs[0] = w11;
    df = (T *) calloc(nx*ny*nz,sizeof(T));
}

template <typename T>
Der<T>::Der(const int _nx, const int _ny, const int _nz, const T _dx, const T _dy, const T _dz, const int _order)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
    dx = _dx;
    dy = _dy;
    dz = _dz;
    order = _order;
    
    /* Check for possibility of interger overflow */
    long int lnx, lny, lnz;
    lnx = nx;
    lny = ny;
    lnz = nz;
    if((lnx*lny*lnz - 1) != ind(nx-1, ny-1, nz-1)) rs_error("Der::Der: The model size is beyond the size this program can model.");
    
    if(order < 1) order = 1;
    if(order > 8) order = 8;
    coeffs = (T *) calloc(CARRAYSIZE,sizeof(T));
    coeffs[0] = w11;
    coeffs[1] = w21;
    coeffs[2] = w22;
    coeffs[3] = w31;
    coeffs[4] = w32;
    coeffs[5] = w33;
    coeffs[6] = w41;
    coeffs[7] = w42;
    coeffs[8] = w43;
    coeffs[9] = w44;
    coeffs[10] = w51;
    coeffs[11] = w52;
    coeffs[12] = w53;
    coeffs[13] = w54;
    coeffs[14] = w55;
    coeffs[15] = w61;
    coeffs[16] = w62;
    coeffs[17] = w63;
    coeffs[18] = w64;
    coeffs[19] = w65;
    coeffs[20] = w66;
    coeffs[21] = w71;
    coeffs[22] = w72;
    coeffs[23] = w73;
    coeffs[24] = w74;
    coeffs[25] = w75;
    coeffs[26] = w76;
    coeffs[27] = w77;
    coeffs[28] = w81;
    coeffs[29] = w82;
    coeffs[30] = w83;
    coeffs[31] = w84;
    coeffs[32] = w85;
    coeffs[33] = w86;
    coeffs[34] = w87;
    coeffs[35] = w88;

    // Allocating the derivative array
    df = (T *) calloc(nx*ny*nz,sizeof(T));
}

/* Destructor*/
template <typename T>
Der<T>::~Der(){
    /* Free allocated variables */
    free(df);
    free(coeffs);
}

template <typename T>
void Der<T>::setOrder(const int _order)
{
    order=_order;
}

/* Derivative functions */
// Forward derivatives
template <typename T>
void Der<T>::ddx_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);

            for(int ix=0; ix < order-1; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                c = coeffs + cind(ix);
                df1[0] = (f1[ind(+1, 0, 0)]- f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<ix+1; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)] - f1[ind(-io, 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }

            c = coeffs + cind(order-1);
            for(int ix=order-1; ix < nx-order-1; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(+1, 0, 0)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)]-f1[ind(-io, 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }

            for(int ix=nx-order-1; ix < nx-1; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                c = coeffs + cind(nx-ix-2);
                df1[0] = (f1[ind(+1, 0, 0)]- f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<nx-ix-1; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)] - f1[ind(-io, 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }
            df1 = df2 + ind(nx-1,0,0);
            df1[0] = 0.0;
        }
    }
}

template <typename T>
void Der<T>::ddy_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;

    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < order-1; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            c = coeffs + cind(iy);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, +1, 0)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<iy+1; io++){
                    df1[0] += (f1[ind(0, io+1, 0)]-f1[ind(0, -io, 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }
        c = coeffs + cind(order-1);
        for(int iy=order-1; iy < ny-order-1; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, +1, 0)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, io+1, 0)]-f1[ind(0, -io, 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }
        for(int iy=ny-order-1; iy < ny-1; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            c = coeffs + cind(ny-iy-2);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, +1, 0)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<ny-iy-1; io++){
                    df1[0] += (f1[ind(0, io+1, 0)]-f1[ind(0, -io, 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }
        df2 = df3 + ind(0,ny-1,0);
        for(int ix=0; ix < nx; ix++){
            df1 = df2 + ind(ix,0,0);
            df1[0] = 0.0;
        }
    }
}

template <typename T>
void Der<T>::ddz_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;
    
    // Compute derivative
    for(int iz=0; iz < order-1; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        c = coeffs + cind(iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<iz+1; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }

    c = coeffs + cind(order-1);
    for(int iz=order-1; iz < nz-order-1; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }

    for(int iz=nz-order-1; iz < nz-1; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        c = coeffs + cind(nz-iz-2);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*c[0];
                for(int io=1; io<nz-iz-1; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }

    df3 = df + ind(0,0,nz-1);
    for(int iy=0; iy < ny; iy++){
	    df2 = df3 + ind(0,iy,0);
	    f2 = f3 + ind(0,iy,0);
	    for(int ix=0; ix < nx; ix++){
		    df1 = df2 + ind(ix,0,0);
		    f1 = f2 + ind(ix,0,0);
		    df1[0] = 0.0;
	    }
    }
}


// Backward derivatives
template <typename T>
void Der<T>::ddx_bw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            df1 = df2;
            df1[0] = 0.0;
            for(int ix=1; ix < order; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                c = coeffs + cind(ix-1);
                df1[0] = (f1[ind(0, 0, 0)]- f1[ind(-1, 0, 0)])*c[0];
                for(int io=1; io<ix; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }

            c = coeffs + cind(order-1);
            for(int ix=order; ix < nx-order; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(-1, 0, 0)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }

            for(int ix=nx-order; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                c = coeffs + cind(nx-ix-1);
                df1[0] = (f1[ind(0, 0, 0)]- f1[ind(-1, 0, 0)])*c[0];
                for(int io=1; io<nx-ix; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*c[io];
                }
                df1[0] /= dx;
            }
        }
    }
}


template <typename T>
void Der<T>::ddy_bw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        df2 = df3;
        for(int ix=0; ix < nx; ix++){
            df1 = df2 + ind(ix,0,0);
            df1[0] = 0.0;
        }
        for(int iy=1; iy < order; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            c = coeffs + cind(iy-1);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, -1, 0)])*c[0];
                for(int io=1; io<iy; io++){
                    df1[0] += (f1[ind(0, io, 0)]-f1[ind(0, -(io+1), 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }

           
        c = coeffs + cind(order-1);
        for(int iy=order; iy < ny-order; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, -1, 0)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, io, 0)]-f1[ind(0, -(io+1), 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }
        for(int iy=ny-order; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            c = coeffs + cind(ny-iy-1);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, -1, 0)])*c[0];
                for(int io=1; io<ny-iy; io++){
                    df1[0] += (f1[ind(0, io, 0)]-f1[ind(0, -(io+1), 0)])*c[io];
                }
                df1[0] /= dy;
            }
        }
    }
}

template <typename T>
void Der<T>::ddz_bw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* c;
   
    // Compute derivative
    df3 = df;
    for(int iy=0; iy < ny; iy++){
	    df2 = df3 + ind(0,iy,0);
	    for(int ix=0; ix < nx; ix++){
		    df1 = df2 + ind(ix,0,0);
		    df1[0] = 0.0;
	    }
    }

    for(int iz=1; iz < order; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        c = coeffs + cind(iz-1);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*c[0];
                for(int io=1; io<iz; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }
 
    c = coeffs + cind(order-1);
    for(int iz=order; iz < nz-order; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*c[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }

    for(int iz=nz-order; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        c = coeffs + cind(nz-iz-1);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*c[0];
                for(int io=1; io<nz-iz; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*c[io];
                }
                df1[0] /= dz;
            }
        }
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Der<float>;
template class Der<double>;
}
