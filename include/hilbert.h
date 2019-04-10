#ifndef HILBERT_H
#define HILBERT_H

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <memory>
#include "utils.h"
#include "fft.h"

namespace rockseis {
/* Hilbert class */
template<typename T>
class Hilbert {
public:
    Hilbert();
    Hilbert(const int _nx, const int _ny, const int _nz);
    ~Hilbert();
    /* Get functions */
    int getNx() { return nx; }		///< Get Nx
    int getNy() { return ny; }		///< Get Ny
    int getNz() { return nz; }		///< Get Nz
    T* getDf() { return df; }           ///< Get Hilbert transform data
    
    // Set functions
    void setNx(const int _nx) { nx = _nx; }	///< Set Nx
    void setNy(const int _ny) { ny = _ny; }	///< Set Ny
    void setNz(const int _nz) { nz = _nz; }	///< Set Nz
    
    /* Hilbert transform functions */
    //Directions
    void hilbertx(T* f); ///< Hilbert transform over x
    void hilberty(T* f); ///< Hilbert transform over x
    void hilbertz(T* f); ///< Hilbert transform over x

    // Hilbert transform
    void Hilbert1D();
    
private:
    int nx;
    int ny;
    int nz;
    T *df; // Array to store the Hilbert
    std::shared_ptr<Fft<T>> fft1d; // Fourier transform object
};


}
#endif //HILBERT_H
