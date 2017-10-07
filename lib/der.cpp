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

#define ind(i,j,k) ((k)*ny*nx + (j)*nx + (i))

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
    coeffs = (T *) calloc(order,sizeof(T));
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
    coeffs = (T *) malloc(order*sizeof(T));
    switch(order){
        case 1:
            coeffs[0] = w11;
            break;
        case 2:
            coeffs[0] = w21;
            coeffs[1] = w22;
            break;
        case 3:
            coeffs[0] = w31;
            coeffs[1] = w32;
            coeffs[2] = w33;
            break;
        case 4:
            coeffs[0] = w41;
            coeffs[1] = w42;
            coeffs[2] = w43;
            coeffs[3] = w44;
            break;
        case 5:
            coeffs[0] = w51;
            coeffs[1] = w52;
            coeffs[2] = w53;
            coeffs[3] = w54;
            coeffs[4] = w55;
            break;
        case 6:
            coeffs[0] = w61;
            coeffs[1] = w62;
            coeffs[2] = w63;
            coeffs[3] = w64;
            coeffs[4] = w65;
            coeffs[5] = w66;
            break;
        case 7:
            coeffs[0] = w71;
            coeffs[1] = w72;
            coeffs[2] = w73;
            coeffs[3] = w74;
            coeffs[4] = w75;
            coeffs[5] = w76;
            coeffs[6] = w77;
            break;
        case 8:
            coeffs[0] = w81;
            coeffs[1] = w82;
            coeffs[2] = w83;
            coeffs[3] = w84;
            coeffs[4] = w85;
            coeffs[5] = w86;
            coeffs[6] = w87;
            coeffs[7] = w88;
            break;
    }
    
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
    // Reallocate and reset coefficients
    free(coeffs);
    coeffs = (T *) malloc(order*sizeof(T));
    switch(order){
        case 1:
            coeffs[0] = w11;
            break;
        case 2:
            coeffs[0] = w21;
            coeffs[1] = w22;
            break;
        case 3:
            coeffs[0] = w31;
            coeffs[1] = w32;
            coeffs[2] = w33;
            break;
        case 4:
            coeffs[0] = w41;
            coeffs[1] = w42;
            coeffs[2] = w43;
            coeffs[3] = w44;
            break;
        case 5:
            coeffs[0] = w51;
            coeffs[1] = w52;
            coeffs[2] = w53;
            coeffs[3] = w54;
            coeffs[4] = w55;
            break;
        case 6:
            coeffs[0] = w61;
            coeffs[1] = w62;
            coeffs[2] = w63;
            coeffs[3] = w64;
            coeffs[4] = w65;
            coeffs[5] = w66;
            break;
        case 7:
            coeffs[0] = w71;
            coeffs[1] = w72;
            coeffs[2] = w73;
            coeffs[3] = w74;
            coeffs[4] = w75;
            coeffs[5] = w76;
            coeffs[6] = w77;
            break;
        case 8:
            coeffs[0] = w81;
            coeffs[1] = w82;
            coeffs[2] = w83;
            coeffs[3] = w84;
            coeffs[4] = w85;
            coeffs[5] = w86;
            coeffs[6] = w87;
            coeffs[7] = w88;
            break;
    }
}

/* Derivative functions */
// Forward derivatives
template <typename T>
void Der<T>::ddx_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    
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
                df1[0] = (f1[ind(+1, 0, 0)]- f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<ix; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)] - f1[ind(-io, 0, 0)])*coeffs[io];
                }
                df1[0] /= dx;
            }

            for(int ix=order-1; ix < nx-order; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(+1, 0, 0)]-f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)]-f1[ind(-io, 0, 0)])*coeffs[io];
                }
                df1[0] /= dx;
            }

            for(int ix=nx-order; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(+1, 0, 0)]- f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<nx-ix; io++){
                    df1[0] += (f1[ind(io+1, 0, 0)] - f1[ind(-io, 0, 0)])*coeffs[io];
                }
                df1[0] /= dx;
            }
        }
    }
}

template <typename T>
void Der<T>::ddy_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    
    // Clear derivative
    for(int iz=0; iz < nz; iz++){
        for(int iy=0; iy < order+1; iy++){
            for(int ix=0; ix < nx; ix++){
                df[ind(ix,iy,iz)] = 0.0;
            }
        }
        for(int iy=ny-order-1; iy < ny; iy++){
            for(int ix=0; ix < nx; ix++){
                df[ind(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=order-1; iy < ny-order; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, +1, 0)]-f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, io+1, 0)]-f1[ind(0, -io, 0)])*coeffs[io];
                }
                df1[0] /= dy;
            }
        }
        
    }
}

template <typename T>
void Der<T>::ddz_fw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    
    // Compute derivative
    for(int iz=0; iz < order-1; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<iz; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*coeffs[io];
                }
                df1[0] /= dz;
            }
        }
    }
    
    for(int iz=order-1; iz < nz-order; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*coeffs[io];
                }
                df1[0] /= dz;
            }
        }
    }

    for(int iz=nz-order; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, +1)]-f1[ind(0, 0, 0)])*coeffs[0];
                for(int io=1; io<nz-iz; io++){
                    df1[0] += (f1[ind(0, 0, io+1)]-f1[ind(0, 0, -io)])*coeffs[io];
                }
                df1[0] /= dz;
            }
        }
    }
}


// Backward derivatives
template <typename T>
void Der<T>::ddx_bw(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);

            for(int ix=0; ix < order; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]- f1[ind(-1, 0, 0)])*coeffs[0];
                for(int io=1; io<ix; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*coeffs[io];
                }
                df1[0] /= dx;
            }

            for(int ix=order; ix < nx-order-1; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(-1, 0, 0)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*coeffs[io];
                }
                df1[0] /= dx;
            }

            for(int ix=nx-order-1; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]- f1[ind(-1, 0, 0)])*coeffs[0];
                for(int io=1; io<nx-ix; io++){
                    df1[0] += (f1[ind(io, 0, 0)]-f1[ind(-(io+1), 0, 0)])*coeffs[io];
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
    
    // Clear derivative
    for(int iz=0; iz < nz; iz++){
        for(int iy=0; iy < order + 1; iy++){
            for(int ix=0; ix < nx; ix++){
                df[ind(ix,iy,iz)] = 0.0;
            }
        }
        for(int iy=ny-order-1; iy < ny; iy++){
            for(int ix=0; ix < nx; ix++){
                df[ind(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    
    // Compute derivative
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=order; iy < ny-order-1; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, -1, 0)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, io, 0)]-f1[ind(0, -(io+1), 0)])*coeffs[io];
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
   
    // Compute derivative
    for(int iz=0; iz < order; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*coeffs[0];
                for(int io=1; io<iz; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*coeffs[io];
                }
                df1[0] /= dz;
            }
        }
    }
 
    for(int iz=order; iz < nz-order-1; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*coeffs[0];
                for(int io=1; io<order; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*coeffs[io];
                }
                df1[0] /= dz;
            }
        }
    }

    for(int iz=nz-order-1; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                f1 = f2 + ind(ix,0,0);
                df1[0] = (f1[ind(0, 0, 0)]-f1[ind(0, 0, -1)])*coeffs[0];
                for(int io=1; io<nz-iz; io++){
                    df1[0] += (f1[ind(0, 0, io)]-f1[ind(0, 0, -(io+1))])*coeffs[io];
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
