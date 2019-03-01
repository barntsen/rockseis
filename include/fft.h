#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "utils.h"

#define TWOPI	(2.0*PI)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define CEIL(x)  ((x) >= 0 ? (int) (x + 1) : (int) (x))

namespace rockseis {
/* Fft class */
template<typename T>
class Fft {
public:
    Fft();
    Fft(const unsigned long _n);
    ~Fft();

    /* Get functions */
    unsigned long getNfft() { return nfft[0]; }		///< Get Nfft for 1D
    unsigned long getNfft(int i) { if(i<3) return nfft[i]; else return 0; }		///< Get Nfft for higher dimensions
    unsigned long getN() { return n[0]; }		///< Get n
    unsigned long getN(int i) { if(i<3) return n[i]; else return 0; }		///< Get Nfft for higher dimensions
    T* getData() { return data; }           ///< Get FFT data

    /* FFT functions */
    void fft1d(int isign); ///< Wrapper to fast fourier trasnform in one dimension

private:
    unsigned long n[3];
    unsigned long nfft[3];
    int ndim;
    T *data; // Array to store the fourier trasnform
    void four1(T *data, unsigned long nn, int isign); ///< Fast fourier trasnform in one dimension
    void fourn(T *data, unsigned long *nn, int ndim, int isign); ///< Fast fourier transform in n dimensions
};


}
#endif //FFT_H
