#include <iostream>
#include <memory>
#include "fft.h"

int main(int argc, char* argv[])
{
	int i;
	int Nx;
	int NFFT;
	double *x;
	double *X;

	/* generate a ramp with 10 numbers */
	Nx = 8;
	printf("Nx = %d\n", Nx);
	x = (double *) malloc(Nx * sizeof(double));
	for(i=0; i<Nx; i++)
	{
		x[i] = i;
	}
    x[0] = 0.;
    x[1] = 0.;
    x[2] = 1.;
    x[3] = 0.;
    x[4] = 2.;
    x[5] = 0.;
    x[6] = 3.;
    x[7] = 0.;

	/* calculate NFFT as the next higher power of 2 >= Nx */
    std::shared_ptr<rockseis::Fft<double>> fftcfg (new rockseis::Fft<double>(Nx));
    NFFT = fftcfg->getNfft();
	printf("NFFT = %d\n", NFFT);
    X = fftcfg->getData();

	/* Storing x(n) in a complex array to make it work with four1. 
	This is needed even though x(n) is purely real in this case. */
	for(i=0; i<Nx; i++)
	{
		X[2*i] = x[i];
		X[2*i+1] = 0.0;
	}
	/* pad the remainder of the array with zeros (0 + 0 j) */
	for(i=Nx; i<NFFT; i++)
	{
		X[2*i] = 0.0;
		X[2*i+1] = 0.0;
	}

	printf("\nInput complex sequence (padded to next highest power of 2):\n");
	for(i=0; i<NFFT; i++)
	{
		printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i], X[2*i+1]);
	}

	/* calculate FFT */
	fftcfg->fft1d(1);

	printf("\nFFT:\n");
	for(i=0; i<NFFT; i++)
	{
		printf("X[%d] = (%.2f + j %.2f)\n", i, X[2*i], X[2*i+1]);
	}

	/* calculate IFFT */
	fftcfg->fft1d(-1);

	/* normalize the IFFT */
	for(i=0; i<NFFT; i++)
	{
		X[2*i] /= NFFT;
		X[2*i+1] /= NFFT;
	}

	printf("\nComplex sequence reconstructed by IFFT:\n");
	for(i=0; i<NFFT; i++)
	{
		printf("x[%d] = (%.2f + j %.2f)\n", i, X[2*i], X[2*i+1]);
	}

	getchar();
}
