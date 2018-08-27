#include "fft.h"
namespace rockseis {
/* Constructors*/
template <typename T>
Fft<T>::Fft(int _n)
{
    n = _n;
    nfft = (int)pow(2.0, ceil(log((double)n)/log(2.0))); 
    data = (T *) calloc(2*nfft,sizeof(T));
}

/* Destructor*/
template <typename T>
Fft<T>::~Fft(){
    /* Free allocated variables */
    free(data);
}


template <typename T>
void Fft<T>::four1(T *data, unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	T wtemp,wr,wpr,wpi,wi,theta;
	T tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(TWOPI/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

template <typename T>
void Fft<T>::fft1d(int isign)
{
    four1(this->data-1, this->nfft, isign);
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Fft<float>;
template class Fft<double>;
}
