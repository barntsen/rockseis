#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "utils.h"

#define TWOPI	(2.0*PI)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr


namespace rockseis {
/* Fft class */
template<typename T>
class Fft {
public:
    Fft();
    Fft(const int _n);
    ~Fft();

    /* Get functions */
    int getNfft() { return nfft; }		///< Get Nfft 
    int getN() { return n; }		///< Get n
    T* getData() { return data; }           ///< Get FFT data


    /* FFT functions */
    void fft1d(int isign); ///< Fast fourier trasnform in one dimension
    

private:
    int n;
    int nfft;
    T *data; // Array to store the fourier trasnform
    void four1(T *data, unsigned long nn, int isign);
};


}
#endif //FFT_H
