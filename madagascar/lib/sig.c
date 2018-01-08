/* 
   Signal processing functions

Author: Wiktor W. Weibull, wiktor.weibull@ntnu.no, NTNU, 2013
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

#include "sig.h" 

/* Functions */ 
/*------------------------------------------------------------*/
float sinc(const float x)
/*< Sinc function >*/
{
	if(x==0.0){
		return (1.0);
	}
        return (sin(PI*x)/(PI*x));
}

void utils_padmodel1d(float *model, float *padded, const int n1, const int pad1)
	/*<Pad input 1d model>*/
{
	int ix;

	for(ix=pad1;ix<n1-pad1; ix++){
		padded[ix]=model[ix-pad1];	
	}

	for(ix=0; ix<pad1; ix++){
		padded[ix]=padded[pad1];
		padded[n1-pad1+ix]=padded[n1-pad1-1];
	}

}

void utils_padmodel2d(float **model, float **padded, const int n1, const int n2, const int pad1, const int pad2)
	/*<Pad input 2d model>*/
{
	int ix,iy;

	for(ix=pad1;ix<n1-pad1; ix++){
		for(iy=pad2;iy<n2-pad2; iy++){
			padded[ix][iy]=model[ix-pad1][iy-pad2];	
		}
	}

	for(ix=0; ix<pad1; ix++){
		for(iy=0; iy<n2; iy++){
			padded[ix][iy]=padded[pad1][iy];
			padded[n1-pad1+ix][iy]=padded[n1-pad1-1][iy];
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<pad2; iy++){
			padded[ix][iy]=padded[ix][pad2];
			padded[ix][n2-pad2+iy]=padded[ix][n2-pad2-1];
		}
	}

}

void utils_padmodel3d(float ***model, float ***padded, const int n1, const int n2, const int n3, const int pad1, const int pad2,const int pad3) 
	/*<Pads input 3d model>*/
{
	int ix,iy,iz;

	for(ix=pad1;ix<n1-pad1; ix++){
		for(iy=pad2;iy<n2-pad2; iy++){
			for(iz=pad3;iz<n3-pad3; iz++){
				padded[ix][iy][iz]=model[ix-pad1][iy-pad2][iz-pad3];	
			}
		}
	}

	for(ix=0; ix<pad1; ix++){
		for(iy=0; iy<n2; iy++){
			for(iz=0; iz<n3; iz++){
				padded[ix][iy][iz]=padded[pad1][iy][iz];
				padded[n1-pad1+ix][iy][iz]=padded[n1-pad1-1][iy][iz];
			}
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<pad2; iy++){
			for(iz=0; iz<n3; iz++){
				padded[ix][iy][iz]=padded[ix][pad2][iz];
				padded[ix][n2-pad2+iy][iz]=padded[ix][n2-pad2-1][iz];
			}
		}
	}
	for(ix=0; ix<n1; ix++){
		for(iy=0; iy<n2; iy++){
			for(iz=0; iz<pad3; iz++){
				padded[ix][iy][iz]=padded[ix][iy][pad3];
				padded[ix][iy][n3-pad3+iz]=padded[ix][iy][n3-pad3-1];
			}
		}
	}

}


sig_cfg_fft1d sig_fft1d_init(int n1_)
	/*< initialize 1d fft >*/
{
	sig_cfg_fft1d fft;
	fft = (sig_cfg_fft1d) sf_alloc(1,sizeof(*fft));

	fft->n1 = n1_;

	fft->ctmp1 = (kiss_fft_cpx*) sf_complexalloc(fft->n1);

	fft->forw1 = kiss_fft_alloc(fft->n1,0,NULL,NULL);
	fft->invs1 = kiss_fft_alloc(fft->n1,1,NULL,NULL);

	if (NULL == fft->forw1 || NULL == fft->invs1) 
		sf_error("%s: KISS FFT allocation error",__FILE__);

	fft->fftscale = 1./sqrtf(fft->n1);

	return fft;
}

void sig_fft1d_close(sig_cfg_fft1d fft)
	/*< close  fft >*/
{
	free (fft->ctmp1);

	free (fft->forw1);
	free (fft->invs1);
}

/*------------------------------------------------------------*/
void sig_fft1d(bool inv          /* inverse/forward flag */, 
		kiss_fft_cpx *pp /* [1...n1] */,
		sig_cfg_fft1d fft) 
/*< apply 1-D FFT >*/
{
	int i1;

	if (inv) {

		/* IFT */
		kiss_fft_stride(fft->invs1,
				pp,
				fft->ctmp1,
				1);

		for (i1=0; i1<fft->n1; i1++) {
			pp[i1] = fft->ctmp1[i1];
		}

		/* scaling */
		for (i1=0; i1<fft->n1; i1++) {
			pp[i1] = sf_crmul(pp[i1],fft->fftscale);
		}

	} else {

		/* scaling */
		for (i1=0; i1<fft->n1; i1++) {
			pp[i1] = sf_crmul(pp[i1],fft->fftscale);
		}

		/* FFT */
		kiss_fft_stride(fft->forw1,
				pp,
				fft->ctmp1, 
				1);

		for (i1=0; i1<fft->n1; i1++) {
			pp[i1] = fft->ctmp1[i1];
		}
	}
}



/*------------------------------------------------------------*/
sig_cfg_fft2d sig_fft2d_init(int n1_, int n2_)
	/*< initialize  2d fft >*/
{
	/*------------------------------------------------------------*/
	sig_cfg_fft2d fft;
	fft = (sig_cfg_fft2d) sf_alloc(1,sizeof(*fft));

	fft->n1 = n1_;
	fft->n2 = n2_;

	fft->ctmp1 = (kiss_fft_cpx*) sf_complexalloc(fft->n1);
	fft->ctmp2 = (kiss_fft_cpx*) sf_complexalloc(fft->n2);

	fft->forw1 = kiss_fft_alloc(fft->n1,0,NULL,NULL);
	fft->invs1 = kiss_fft_alloc(fft->n1,1,NULL,NULL);

	fft->forw2 = kiss_fft_alloc(fft->n2,0,NULL,NULL);
	fft->invs2 = kiss_fft_alloc(fft->n2,1,NULL,NULL);

	if (NULL == fft->forw2 || NULL == fft->invs2 || 
			NULL == fft->forw1 || NULL == fft->invs1) 
		sf_error("%s: KISS FFT allocation error",__FILE__);

	fft->fftscale = 1./sqrtf(fft->n1*fft->n2);

	return fft;
}

/*------------------------------------------------------------*/
void sig_fft2d_close(sig_cfg_fft2d fft)
	/*< close  fft >*/
{
	free (fft->ctmp1);
	free (fft->ctmp2);

	free (fft->forw1);
	free (fft->invs1);

	free (fft->forw2);
	free (fft->invs2);
}

/*------------------------------------------------------------*/
void sig_fft2d(bool inv          /* inverse/forward flag */, 
		kiss_fft_cpx **pp /* [1...n2][1...n1] */,
		sig_cfg_fft2d fft) 
/*< apply 2-D FFT >*/
{
	int i1,i2;

	if (inv) {

		/* IFT 1 */
		for (i2=0; i2 < fft->n2; i2++) {
			kiss_fft_stride(fft->invs1,
					pp[i2],
					fft->ctmp1,
					1);

			for (i1=0; i1<fft->n1; i1++) {
				pp[i2][i1] = fft->ctmp1[i1];
			}
		}

		/* IFT 2 */
		for (i1=0; i1 < fft->n1; i1++) {
			kiss_fft_stride(fft->invs2,
					pp[0]+i1,
					fft->ctmp2,
					fft->n1);

			for (i2=0; i2<fft->n2; i2++) {
				pp[i2][i1] = fft->ctmp2[i2];
			}
		}

		/* scaling */
		for     (i2=0; i2<fft->n2; i2++) {
			for (i1=0; i1<fft->n1; i1++) {
				pp[i2][i1] = sf_crmul(pp[i2][i1],fft->fftscale);
			}
		}

	} else {

		/* scaling */
		for     (i2=0; i2<fft->n2; i2++) {
			for (i1=0; i1<fft->n1; i1++) {
				pp[i2][i1] = sf_crmul(pp[i2][i1],fft->fftscale);
			}
		}

		/* FFT 2 */
		for (i1=0; i1 < fft->n1; i1++) {
			kiss_fft_stride(fft->forw2,
					pp[0]+i1,
					fft->ctmp2, 
					fft->n1);

			for (i2=0; i2<fft->n2; i2++) {
				pp[i2][i1] = fft->ctmp2[i2];
			}
		}

		/* FFT 1 */
		for (i2=0; i2 < fft->n2; i2++) {
			kiss_fft_stride(fft->forw1,
					pp[i2],
					fft->ctmp1,
					1);

			for (i1=0; i1<fft->n1; i1++) {
				pp[i2][i1] = fft->ctmp1[i1];
			}
		}

	}

}

void sig_taper(float *data, int n, int ntaper)
	/*<Apply cosine taper on both ends of a trace>*/
{
	int i;
	float ifloat;
	float filter;
	for (i=0; i<ntaper; i++){
		ifloat=(float) i;
		filter = 0.5*(1.0-cosf(PI*ifloat/(ntaper-1)));
		filter=filter*filter;
		data[i] *= filter;
		data[n-1-i] *= filter;
	}
}

void sig_smoother1d(float *data, int n1, int pad, float perc1)
	/*<Apply 2D gaussian low-pass filter in the wavenumber domain>*/
{
	int L1,L2;
	float *H1;
	float *wrk;
	kiss_fft_cpx *cwrk;
	float *tmp;
	int fn1, i1; 
	float ifloat;
	sig_cfg_fft1d cfg; 

	/* Find next fast even (x2) dimension */
	fn1=2*kiss_fft_next_fast_size(((n1+2*pad)+1)/2);

	if(fn1%2) fn1++; // Make it even

	// Initialize fft struct
	cfg= sig_fft1d_init(fn1);

	cwrk= (kiss_fft_cpx*) sf_complexalloc(fn1);
	for (i1=0; i1<fn1; i1++){
		cwrk[i1].r= (kiss_fft_scalar) 0.0;
		cwrk[i1].i= (kiss_fft_scalar) 0.0;
	}

	//Pad initial model to avoid edge effects
	wrk=sf_floatalloc(n1+2*pad);
	utils_padmodel1d(data, wrk, n1+2*pad, pad);

	//Compute the filter in dimension 1
	H1=sf_floatalloc(fn1); // Allocate memory

	L1=floor(perc1*fn1);
	L2=ceil((fn1-L1)/2);
	L1=fn1-2*L2;
	for (i1=0;i1<fn1; i1++){
		H1[i1]=0.0;
	}
	for (i1=L1+L2;i1<fn1;i1++){
		H1[i1]=0.0;
	}
	for (i1=0; i1<L1; i1++){
		ifloat= (float) i1;
		H1[L2+i1]=0.5*(1.0-cosf(2.0*PI*ifloat/(L1-1)));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn1);
	for(i1=0;i1<fn1/2;i1++){
		tmp[i1]=H1[(fn1/2)+i1];
	}
	for(i1=fn1/2;i1<fn1;i1++){
		tmp[i1]=H1[i1-(fn1/2)];
	}
	for (i1=0;i1<fn1;i1++){
		H1[i1]=tmp[i1];
	}
	free(tmp);

	// Copy data into cwrk array 
	for (i1=0; i1<n1+2*pad; i1++){
		cwrk[i1].r = (kiss_fft_scalar ) wrk[i1];
		cwrk[i1].i= (kiss_fft_scalar ) 0.0;
	}

	/* Apply forward fourier transform */
	sig_fft1d(0,cwrk,cfg);

	/* Apply filter */
	for (i1=0; i1<fn1; i1++){
		cwrk[i1]=sf_crmul(cwrk[i1], H1[i1]);
	}

	/* Apply backward fourier transform */
	sig_fft1d(1,cwrk,cfg);

	/* Take the real part */
	for (i1=0; i1<n1; i1++){
		data[i1]=(float) cwrk[i1+pad].r;
	}

	/* Free memory */
	sig_fft1d_close(cfg);
	free(cwrk);
	free(wrk);
	free(H1);
}


void sig_smoother2d(float **data, int n1, int n2, int pad, float perc1, float perc2)
	/*<Apply 2D gaussian low-pass filter in the wavenumber domain>*/
{
	int L1,L2;
	float *H1, *H2;
	float **wrk;
	kiss_fft_cpx **cwrk;
	float *tmp;
	int fn1, fn2, i1, j1; 
	float ifloat;
	sig_cfg_fft2d cfg; 

	/* Find next fast even (x2) dimension */
	fn1=2*kiss_fft_next_fast_size(((n1+2*pad)+1)/2);
	fn2=2*kiss_fft_next_fast_size(((n2+2*pad)+1)/2);

	if(fn1%2) fn1++; // Make it even
	if(fn2%2) fn2++; // Make it even

	// Initialize fft struct
	cfg= sig_fft2d_init(fn1,fn2);

	cwrk= (kiss_fft_cpx**) sf_complexalloc2(fn1,fn2);
	for (i1=0; i1<fn2; i1++){
		for (j1=0; j1<fn1; j1++){
			cwrk[i1][j1].r= (kiss_fft_scalar) 0.0;
			cwrk[i1][j1].i= (kiss_fft_scalar) 0.0;
		}
	}

	//Pad initial model to avoid edge effects
	wrk=sf_floatalloc2(n1+2*pad,n2+2*pad);
	utils_padmodel2d(data, wrk, n2+2*pad, n1+2*pad, pad, pad);

	//Compute the filter in dimension 1
	H1=sf_floatalloc(fn1); // Allocate memory

	L1=floor(perc1*fn1);
	L2=ceil((fn1-L1)/2);
	L1=fn1-2*L2;
	for (i1=0;i1<fn1; i1++){
		H1[i1]=0.0;
	}
	for (i1=L1+L2;i1<fn1;i1++){
		H1[i1]=0.0;
	}
	for (i1=0; i1<L1; i1++){
		ifloat= (float) i1;
		H1[L2+i1]=0.5*(1.0-cosf(2.0*PI*ifloat/(L1-1)));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn1);
	for(i1=0;i1<fn1/2;i1++){
		tmp[i1]=H1[(fn1/2)+i1];
	}
	for(i1=fn1/2;i1<fn1;i1++){
		tmp[i1]=H1[i1-(fn1/2)];
	}
	for (i1=0;i1<fn1;i1++){
		H1[i1]=tmp[i1];
	}
	free(tmp);

	//Compute the filter in dimension 2
	H2=sf_floatalloc(fn2); // Allocate memory

	L1=floor(perc2*fn2);
	L2=ceil((fn2-L1)/2);
	L1=fn2-2*L2;
	for (i1=0;i1<fn2; i1++){
		H2[i1]=0.0;
	}
	for (i1=L1+L2;i1<fn2;i1++){
		H2[i1]=0.0;
	}
	for (i1=0; i1<L1; i1++){
		ifloat= (float) i1;
		H2[L2+i1]=0.5*(1.0-cosf(2.0*PI*ifloat/(L1-1)));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn2);
	for(i1=0;i1<fn2/2;i1++){
		tmp[i1]=H2[(fn2/2)+i1];
	}
	for(i1=fn2/2;i1<fn2;i1++){
		tmp[i1]=H2[i1-(fn2/2)];
	}
	for (i1=0;i1<fn2;i1++){
		H2[i1]=tmp[i1];
	}
	free(tmp); 

	// Copy data into cwrk array 
	for (i1=0; i1<n2+2*pad; i1++){
		for (j1=0; j1<n1+2*pad; j1++){
			cwrk[i1][j1].r = (kiss_fft_scalar ) wrk[i1][j1];
			cwrk[i1][j1].i= (kiss_fft_scalar ) 0.0;
		}
	}

	/* Apply forward fourier transform */
	sig_fft2d(0,cwrk,cfg);

	/* Apply filter */
	for (i1=0; i1<fn2; i1++){
		for (j1=0; j1<fn1; j1++){
			cwrk[i1][j1]=sf_crmul(cwrk[i1][j1], H1[j1]*H2[i1]);
		}
	}

	/* Apply backward fourier transform */
	sig_fft2d(1,cwrk,cfg);

	/* Take the real part */
	for (i1=0; i1<n2; i1++){
		for (j1=0; j1<n1; j1++){
			data[i1][j1]=(float) cwrk[i1+pad][j1+pad].r;
		}
	}

	/* Free memory */
	sig_fft2d_close(cfg);
	free(cwrk[0]);
	free(cwrk);
	free(wrk[0]);
	free(wrk);
	free(H1);
	free(H2);

}

/*************************
 *        3d             *
 *************************/

/*------------------------------------------------------------*/
sig_cfg_fft3d sig_fft3d_init(int n1_, int n2_, int n3_)
	/*< initialize  3d fft >*/
{
	/*------------------------------------------------------------*/
	sig_cfg_fft3d fft;
	fft = (sig_cfg_fft3d) sf_alloc(1,sizeof(*fft));

	fft->n1 = n1_;
	fft->n2 = n2_;
	fft->n3 = n3_;

	fft->ctmp1 = (kiss_fft_cpx*) sf_complexalloc(fft->n1);
	fft->ctmp2 = (kiss_fft_cpx*) sf_complexalloc(fft->n2);
	fft->ctmp3 = (kiss_fft_cpx*) sf_complexalloc(fft->n3);

	fft->forw1 = kiss_fft_alloc(fft->n1,0,NULL,NULL);
	fft->invs1 = kiss_fft_alloc(fft->n1,1,NULL,NULL);

	fft->forw2 = kiss_fft_alloc(fft->n2,0,NULL,NULL);
	fft->invs2 = kiss_fft_alloc(fft->n2,1,NULL,NULL);

	fft->forw3 = kiss_fft_alloc(fft->n3,0,NULL,NULL);
	fft->invs3 = kiss_fft_alloc(fft->n3,1,NULL,NULL);

	if (NULL == fft->forw3 || NULL == fft->invs3 || 
			NULL == fft->forw2 || NULL == fft->invs2 || 
			NULL == fft->forw1 || NULL == fft->invs1) 
		sf_error("%s: KISS FFT allocation error",__FILE__);

	fft->fftscale = 1./sqrtf(fft->n1*fft->n2*fft->n3);

	return fft;
}

/*------------------------------------------------------------*/
void sig_fft3d_close(sig_cfg_fft3d fft)
	/*< close  fft >*/
{
	free (fft->ctmp1);
	free (fft->ctmp2);
	free (fft->ctmp3);

	free (fft->forw1);
	free (fft->invs1);

	free (fft->forw2);
	free (fft->invs2);

	free (fft->forw3);
	free (fft->invs3);
}

/*------------------------------------------------------------*/
void sig_fft3d(bool inv          /* inverse/forward flag */, 
		kiss_fft_cpx ***pp /* [1...n3][1...n2][1...n1] */,
		sig_cfg_fft3d fft) 
/*< apply 3-D FFT >*/
{
	int i1,i2,i3;

	if (inv) {

		/* IFT 1 */
		for(i3=0; i3<fft->n3; i3++){
			for (i2=0; i2 < fft->n2; i2++) {
				kiss_fft_stride(fft->invs1,
						pp[i3][i2],
						fft->ctmp1,
						1);

				for (i1=0; i1<fft->n1; i1++) {
					pp[i3][i2][i1] = fft->ctmp1[i1];
				}
			}
		}

		/* IFT 2 */
		for(i3=0; i3<fft->n3; i3++){
			for (i1=0; i1 < fft->n1; i1++) {
				kiss_fft_stride(fft->invs2,
						pp[i3][0]+i1,
						fft->ctmp2,
						fft->n1);

				for (i2=0; i2<fft->n2; i2++) {
					pp[i3][i2][i1] = fft->ctmp2[i2];
				}
			}
		}

		/* IFT 3 */
		for(i2=0; i2<fft->n2; i2++){
			for (i1=0; i1 < fft->n1; i1++) {
				kiss_fft_stride(fft->invs3,
						pp[0][0]+fft->n1*i2+i1,
						fft->ctmp3,
						fft->n1*fft->n2);

				for (i3=0; i3<fft->n3; i3++) {
					pp[i3][i2][i1] = fft->ctmp3[i3];
				}
			}
		}

		/* scaling */
		for     (i3=0; i3<fft->n3; i3++) {
			for     (i2=0; i2<fft->n2; i2++) {
				for (i1=0; i1<fft->n1; i1++) {
					pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->fftscale);
				}
			}
		}

	} else {

		/* scaling */
		for     (i3=0; i3<fft->n3; i3++) {
			for     (i2=0; i2<fft->n2; i2++) {
				for (i1=0; i1<fft->n1; i1++) {
					pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fft->fftscale);
				}
			}
		}

		/* FFT 3 */
		for(i2=0; i2<fft->n2; i2++){
			for (i1=0; i1 < fft->n1; i1++) {
				kiss_fft_stride(fft->forw3,
						pp[0][0]+fft->n1*i2+i1,
						fft->ctmp3,
						fft->n1*fft->n2);

				for (i3=0; i3<fft->n3; i3++) {
					pp[i3][i2][i1] = fft->ctmp3[i3];
				}
			}
		}
		/* FFT 2 */
		for(i3=0; i3<fft->n3; i3++){
			for (i1=0; i1 < fft->n1; i1++) {
				kiss_fft_stride(fft->forw2,
						pp[i3][0]+i1,
						fft->ctmp2,
						fft->n1);

				for (i2=0; i2<fft->n2; i2++) {
					pp[i3][i2][i1] = fft->ctmp2[i2];
				}
			}
		}

		/* FFT 1 */
		for(i3=0; i3<fft->n3; i3++){
			for (i2=0; i2 < fft->n2; i2++) {
				kiss_fft_stride(fft->forw1,
						pp[i3][i2],
						fft->ctmp1,
						1);

				for (i1=0; i1<fft->n1; i1++) {
					pp[i3][i2][i1] = fft->ctmp1[i1];
				}
			}
		}
	}

}

void sig_smoother3d(float ***data, int n1, int n2, int n3, int pad, float perc1, float perc2, float perc3)
	/*<Apply 3D gaussian low-pass filter in the frequency domain>*/
{
	int L1, L2;
	float *H1, *H2, *H3;
	float ***wrk;
	kiss_fft_cpx ***cwrk;
	float *tmp;
	int fn1, fn2, fn3, i, j, k; 
	sig_cfg_fft3d cfg; 

	/* Find next fast even (x2) dimension */
	fn1=2*kiss_fft_next_fast_size(((n1+2*pad)+1)/2);
	fn2=2*kiss_fft_next_fast_size(((n2+2*pad)+1)/2);
	fn3=2*kiss_fft_next_fast_size(((n3+2*pad)+1)/2);

	if(fn1%2) fn1++; // Make it even
	if(fn2%2) fn2++; // Make it even
	if(fn3%2) fn3++; // Make it even


	//Pad initial model to avoid edge effects
	wrk=sf_floatalloc3(n1+2*pad,n2+2*pad, n3+2*pad);
	utils_padmodel3d(data, wrk, n3+2*pad, n2+2*pad, n1+2*pad, pad, pad, pad);

	//Compute the filter in dimension 1
	H1=sf_floatalloc(fn1); // Allocate memory

	L1=floor(perc1*fn1);
	L2=ceil((fn1-L1)/2);
	L1=fn1-2*L2;
	for (i=0;i<fn1; i++){
		H1[i]=0.0;
	}
	for (i=L1+L2;i<fn1;i++){
		H1[i]=0.0;
	}
	for (i=0; i<L1; i++){
		H1[L2+i]=0.5*(1.0-cos(2.0*PI*((float) i/(L1-1))));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn1);
	for(i=0;i<fn1/2;i++){
		tmp[i]=H1[(fn1/2)+i];
	}
	for(i=fn1/2;i<fn1;i++){
		tmp[i]=H1[i-(fn1/2)];
	}
	for (i=0;i<fn1;i++){
		H1[i]=tmp[i];
	}
	free(tmp);

	//Compute the filter in dimension 2
	H2=sf_floatalloc(fn2); // Allocate memory

	L1=floor(perc2*fn2);
	L2=ceil((fn2-L1)/2);
	L1=fn2-2*L2;
	for (i=0;i<fn2; i++){
		H2[i]=0.0;
	}
	for (i=L1+L2;i<fn2;i++){
		H2[i]=0.0;
	}
	for (i=0; i<L1; i++){
		H2[L2+i]=0.5*(1.0-cos(2.0*PI*((float) i/(L1-1))));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn2);
	for(i=0;i<fn2/2;i++){
		tmp[i]=H2[(fn2/2)+i];
	}
	for(i=fn2/2;i<fn2;i++){
		tmp[i]=H2[i-(fn2/2)];
	}
	for (i=0;i<fn2;i++){
		H2[i]=tmp[i];
	}
	free(tmp); 

	//Compute the filter in dimension 3
	H3=sf_floatalloc(fn3); // Allocate memory

	L1=floor(perc3*fn3);
	L2=ceil((fn3-L1)/2);
	L1=fn3-2*L2;
	for (i=0;i<fn3; i++){
		H3[i]=0.0;
	}
	for (i=L1+L2;i<fn3;i++){
		H3[i]=0.0;
	}
	for (i=0; i<L1; i++){
		H3[L2+i]=0.5*(1.0-cos(2.0*PI*((float) i/(L1-1))));
	}

	// Transpose filter to correct position after fft
	tmp=sf_floatalloc(fn3);
	for(i=0;i<fn3/2;i++){
		tmp[i]=H3[(fn3/2)+i];
	}
	for(i=fn3/2;i<fn3;i++){
		tmp[i]=H3[i-(fn3/2)];
	}
	for (i=0;i<fn3;i++){
		H3[i]=tmp[i];
	}
	free(tmp);


	// Initialize fft struct
	cfg= sig_fft3d_init(fn1,fn2,fn3);

	// Allocate complex wrk array
	cwrk= (kiss_fft_cpx***) sf_complexalloc3(fn1,fn2,fn3);
	for (i=0; i<fn3; i++){
		for (j=0; j<fn2; j++){
			for (k=0; k<fn1; k++){
				cwrk[i][j][k].r=0.0;
				cwrk[i][j][k].i=0.0;
			}
		}
	}

	// Copy data into cwrk array 
	for (i=0; i<n3+2*pad; i++){
		for (j=0; j<n2+2*pad; j++){
			for (k=0; k<n1+2*pad; k++){
				cwrk[i][j][k].r = wrk[i][j][k];
				cwrk[i][j][k].i= 0.0;
			}
		}
	}
	// Free work array
	free(wrk[0][0]);
	free(wrk[0]);
	free(wrk);

	/* Apply forward fourier transform */
	sig_fft3d(0,cwrk,cfg);

	/* Apply filter */
	for (i=0; i<fn3; i++){
		for (j=0; j<fn2; j++){
			for (k=0; k<fn1; k++){
				cwrk[i][j][k]=sf_crmul(cwrk[i][j][k], H1[k]*H2[j]*H3[i]);
			}
		}
	}

	/* Apply backward fourier transform */
	sig_fft3d(1,cwrk,cfg);

	/* Take the real part */
	for (i=0; i<n3; i++){
		for (j=0; j<n2; j++){
			for (k=0; k<n1; k++){
				data[i][j][k]=cwrk[i+pad][j+pad][k+pad].r;
			}
		}
	}

	/* Free memory */
	sig_fft3d_close(cfg);
	free(cwrk[0][0]);
	free(cwrk[0]);
	free(cwrk);
	free(H1);
	free(H2);
	free(H3);
}

void sig_gauss(float *pulse, float f0, float t0, float dt, int nt, int order)
	/*< Creates 0th, 1th or 2nd order derivative of Gaussian wavelet >*/
{
	int i;
	float t;
	float s=-PI*PI*f0*f0;

	for(i=0; i<nt; i++){
		t = i*dt - t0;
		switch (order){
			case 0:
				pulse[i]=t*expf(s*t*t);
				break;
			case 1:
				pulse[i]=(1.0+2.0*s*t*t)*expf(s*t*t);
				break;
			case 2:
				pulse[i]=-2.0*s*t*((1+2*s*t*t)+2.0)*expf(s*t*t);
				break;
			default:
				pulse[i]=-2.0*s*t*((1+2*s*t*t)+2.0)*expf(s*t*t);
				break;
		}
	}
}

void sig_1Dsinc_lanczos(float *op, int N, float s)
	/*< Creates Sinc-lanczos shift operator>*/
{
	int Nc=(N-1)/2;   // Center index of the operator
	int i;
	float arg;

	for(i=0; i<N; i++){
		arg=(Nc-(i-s));
		op[i]=sinc(arg)*sinc(arg/Nc);
	}
}

void sig_tukey(float *pulse, float f0, float f1, float f2, float f3, float t0, double dt, int nt)
/*< Creates wavelet with flat pass band>*/
{
	int i;
	double f,df;
	int nf,nfs;
	sig_cfg_fft1d cfg; 

	/* Compute size of complex array */
	nf=2*kiss_fft_next_fast_size(nt);
	cfg= sig_fft1d_init(nf);

	nfs=nf/2 + 1;
	df=(1.0/(2.0*dt))/nfs;
	kiss_fft_cpx *W;
	W = (kiss_fft_cpx *) calloc(nf, sizeof(kiss_fft_cpx));

	/* Compute wavelet spectrum  */
	for(i=0; i<nfs; i++)
	{
		f = i*df;
		if(f < f0) W[i].r = 0.0;
		if(f>= f0 && f < f1) W[i].r = 0.5*(1.0 - cosf(PI*(f-f0)/(f1-f0))); 
		if(f >= f1 && f <= f2 ) W[i].r = 1;
		if(f > f2)  W[i].r = 0.5*(1.0 - cosf(PI*(f3-f)/(f3-f2))); 
		if(f > f3) W[i].r = 0;
	}

	/* Apply time delay */
	kiss_fft_cpx shift;
	for(i=0; i<nfs; i++)
	{
		f = i*df;
		shift.r = cosf(-2.0*PI*f*t0);
		shift.i = sinf(-2.0*PI*f*t0);
		W[i] = sf_cmul(W[i], shift);
	}

	/* Reconstruct negative part of frequency spectrum */
	for(i=nfs-2; i>0; i--)
	{
		W[nfs-2 - i + nfs].r = W[i].r;  
		W[nfs-2 - i + nfs].i = -1.0*W[i].i;  
	}


	/* Apply backward fourier transform */
	sig_fft1d(1,W,cfg);
	for(i=0; i<nt; i++)
	{
		pulse[i] = W[i].r;
	}
	free(W);
}

void sig_filt(float *pulse, float f0, float f1, float f2, float f3, double dt, int nt)
/*< Filters with a hamming window>*/
{
	int i;
	double f,df;
	int nf,nfs;
	sig_cfg_fft1d cfg; 

	/* Compute size of complex array */
	nf=2*kiss_fft_next_fast_size(nt);
	cfg= sig_fft1d_init(nf);

	nfs=nf/2 + 1;
	df=(1.0/(2.0*dt))/nfs;
	kiss_fft_cpx *W, *cdata;
	W = (kiss_fft_cpx *) calloc(nf, sizeof(kiss_fft_cpx));
	cdata = (kiss_fft_cpx *) calloc(nf, sizeof(kiss_fft_cpx));
    for(i=0; i<nt; i++){
        cdata[i].r = pulse[i];
    }

	/* Apply forward fourier transform */
	sig_fft1d(0,cdata,cfg);

	/* Compute window spectrum  */
	for(i=0; i<nfs; i++)
	{
		f = i*df;
		if(f < f0) W[i].r = 0.0;
		if(f>= f0 && f < f1) W[i].r = 0.5*(1.0 - cosf(PI*(f-f0)/(f1-f0))); 
		if(f >= f1 && f <= f2 ) W[i].r = 1;
		if(f > f2)  W[i].r = 0.5*(1.0 - cosf(PI*(f3-f)/(f3-f2))); 
		if(f > f3) W[i].r = 0;
	}

	/* Reconstruct negative part of frequency spectrum */
	for(i=nfs-2; i>0; i--)
	{
		W[nfs-2 - i + nfs].r = W[i].r;  
		W[nfs-2 - i + nfs].i = -1.0*W[i].i;  
	}

    // Apply filter
	for(i=0; i<nf; i++)
    {
		W[i] = sf_cmul(W[i], cdata[i]);
    }

	/* Apply backward fourier transform */
	sig_fft1d(1,W,cfg);
	for(i=0; i<nt; i++)
	{
		pulse[i] = W[i].r;
	}
	free(W);
	free(cdata);
}


