#include "hilbert.h"

#define CARRAYSIZE 36 // Array of coefficients
#define ind(i,j,k) ((k)*ny*nx + (j)*nx + (i))
#define cind(i) (((i)*(i+1))/2)
#define MAX(x,y) (((x)>=(y))?(x):(y))

namespace rockseis {
/* Constructors*/
template <typename T>
Hilbert<T>::Hilbert()
{
    nx = 1;
    ny = 1;
    nz = 1;
    df = (T *) calloc(nx*ny*nz,sizeof(T));
}

template <typename T>
Hilbert<T>::Hilbert(const int _nx, const int _ny, const int _nz)
{
    nx = _nx;
    ny = _ny;
    nz = _nz;
    
    /* Check for possibility of interger overflow */
    long int lnx, lny, lnz;
    lnx = nx;
    lny = ny;
    lnz = nz;
    if((lnx*lny*lnz - 1) != ind(nx-1, ny-1, nz-1)) rs_error("Hilbert::Hilbert: The model size is beyond the size this program can handle.");
    
       // Allocating the hilbert transform array
    df = (T *) calloc(nx*ny*nz,sizeof(T));

    // Allocating Fourier transform
    long int nmax = 0;
    nmax = MAX(nx, MAX(ny,nz));
    fft1d = std::make_shared<Fft<T>>(nmax);
}

/* Destructor*/
template <typename T>
Hilbert<T>::~Hilbert(){
    /* Free allocated variables */
    free(df);
}

template<typename T>
void Hilbert<T>::Hilbert1D()
/*< Compute the Hilbert transform of a trace */
{
	unsigned long i;

	/* Compute size of complex array */
	unsigned long int nf=fft1d->getNfft();
	unsigned long int nfs=nf/2;
	
	/* Apply forward fourier transform */
	fft1d->fft1d(1);

	T *cdata;
	cdata = fft1d->getData();
	/* Zero out negative frequencies */
	for(i=nfs+1; i<nf; i++)
	{
        cdata[2*i] = 0.0;
        cdata[2*i+1] = 0.0;
	}

	/* Apply backward fourier transform */
	fft1d->fft1d(-1);
}

/* Hilbert transform functions */
// Forward derivatives
template <typename T>
void Hilbert<T>::hilbertx(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* cdata = fft1d->getData();
    long int nf = fft1d->getNfft();

    // Compute Hilbert transform
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        for(int iy=0; iy < ny; iy++){
            df2 = df3 + ind(0,iy,0);
            f2 = f3 + ind(0,iy,0);

            for(int ix=0; ix < nx; ix++){
                f1 = f2 + ind(ix,0,0);
                cdata[2*ix] =  f1[ind(0, 0, 0)];
                cdata[2*ix+1] =  0.0;
            }
            for(int ix=nx; ix < nf; ix++){
                cdata[2*ix] =  0.0;
                cdata[2*ix+1] =  0.0;
            }

            // Transform cdata
            Hilbert1D();

            /* Get imaginary part (Hilbert transform*0.5) */
            for(int ix=0; ix < nx; ix++){
                df1 = df2 + ind(ix,0,0);
                df1[0] = 2.0*cdata[2*ix+1]/nf;
            }
        }
    }
}

template <typename T>
void Hilbert<T>::hilberty(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* cdata = fft1d->getData();
    long int nf = fft1d->getNfft();

    // Compute Hilbert transform
    for(int iz=0; iz < nz; iz++){
        df3 = df + ind(0,0,iz);
        f3 = f + ind(0,0,iz);
        int iy = 0;
        df2 = df3 + ind(0,iy,0);
        f2 = f3 + ind(0,iy,0);
        for(int ix=0; ix < nx; ix++){
            df1 = df2 + ind(ix,0,0);
            f1 = f2 + ind(ix,0,0);
            for(iy=0; iy < ny; iy++){
                cdata[2*iy] =  f1[ind(0, iy, 0)];
                cdata[2*iy+1] =  0.0;
            }
            for(iy=ny; iy < nf; iy++){
                cdata[2*iy] =  0.0;
                cdata[2*iy+1] =  0.0;
            }

            // Transform cdata
            Hilbert1D();

            /* Get imaginary part (Hilbert transform*0.5) */
            for(iy=0; iy < ny; iy++){
                df1[ind(0, iy, 0)] = 2.0*cdata[2*iy+1]/nf;
            }
        }
    }
}

template <typename T>
void Hilbert<T>::hilbertz(T *f){
    T *df1, *f1;
    T* df2, *f2;
    T* df3, *f3;
    T* cdata = fft1d->getData();
    long int nf = fft1d->getNfft();
    
    // Compute derivative
    int iz=0;
    df3 = df + ind(0,0,iz);
    f3 = f + ind(0,0,iz);
    for(int iy=0; iy < ny; iy++){
        df2 = df3 + ind(0,iy,0);
        f2 = f3 + ind(0,iy,0);
        for(int ix=0; ix < nx; ix++){
            df1 = df2 + ind(ix,0,0);
            f1 = f2 + ind(ix,0,0);
            for(iz=0; iz<nz; iz++){
                cdata[2*iz] =  f1[ind(0, 0, iz)];
                cdata[2*iz+1] =  0.0;
            }
            for(iz=nz; iz < nf; iz++){
                cdata[2*iz] =  0.0;
                cdata[2*iz+1] =  0.0;
            }
            // Transform cdata
            Hilbert1D();

            /* Get imaginary part (Hilbert transform*0.5) */
            for(iz=0; iz < nz; iz++){
                df1[ind(0, 0, iz)] = 2.0*cdata[2*iz+1]/nf;
            }
        }
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Hilbert<float>;
template class Hilbert<double>;
}
