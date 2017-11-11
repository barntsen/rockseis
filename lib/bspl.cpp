#include "bspl.h"

namespace rockseis {
/* Constructors*/
template<typename T>
Bspl<T>::Bspl()
{
// Do nothing
}

template<typename T>
Bspl<T>::~Bspl()
{
// Do nothing
}

template<typename T>
void Bspl<T>::bspl(T *h, T *t, int k, T x, int l)
/*<evaluates the (k+1) non-zero b-splines of degree k at t(l) <= x < t(l+1) using the stable recurrence relation of de Boor and Cox.>*/
{
	T hh[6];
	int li, lj;
	T f;
	T one = 1.0; /* 0.1e+01; */ 
	
	h[0] = one;
	int i, j;
	for (j=1; j<=k; j++){ 
	    for (i=1; i<=j; i++){
	            hh[i-1] = h[i-1];
	    }
	    h[0] = 0.0;
	    for(i=1; i<=j; i++){
		li = l+i;
		lj = li-j;
		f = hh[i-1]/(t[li-1]-t[lj-1]);
		h[i-1] = h[i-1]+f*(t[li-1]-x);
		h[i] = f*(x-t[lj-1]);
	    }
	}
}


/* Bspl1D class */
template<typename T>
Bspl1D<T>::Bspl1D(const int _nx, const T _dx, const T _dtx, const int _kx) : Bspl<T>()
{
	/* Setting size values for the evaluation of the spline */
	mx=_nx;
	dtx=_dtx;
	kx=_kx;
	
	/* Compute size of knot array */
	double xmax=mx*_dx;
	double xmin=-_dx;

	ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(mx%2==0){
		if(ntx%2 != 0){
			ntx++;
			dtx=(xmax-xmin)/(ntx-2*kx-1);
		}
	}else{
		if(ntx%2 == 0){
			ntx++;
			dtx=(xmax-xmin)/(ntx-2*kx-1);
		}
	}

	nc = (ntx-kx-1);
	c = (T *) calloc(nc,sizeof(T));
	tx = (T *) calloc(ntx,sizeof(T));
	x = (T *) calloc(mx,sizeof(T));
	mod = (T *) calloc(mx,sizeof(T));

	/* Computing the evaluation arrays */
	int i;
	for (i=0; i<mx; i++){
		x[i]=i*_dx;
	}

	/* Computing the Knotts array */
	for (i=0; i<kx+1; i++){
		tx[i]=xmin;
	}

	for (i=kx+1; i<ntx-kx-1; i++){
		tx[i]=xmin+(i-kx)*dtx;
	}

	for (i=ntx-kx-1; i<ntx; i++){
		tx[i]=xmax;
	}
}

template<typename T>
Bspl1D<T>::~Bspl1D()
{
	free(c);
	free(tx);
	free(x);
	free(mod);
}

template<typename T>
void Bspl1D<T>::bisp()
/*<Evaluates the spline coefficients in c over z.>*/
{

	int i;
	int i1, l, l1;
	T arg;
	T h[6];
    for(i1=0; i1<6; i1++) h[1]=0.0;
	T sp = 0.0;

	int nx=this->ntx;
	int kx=this->kx;
	int mx=this->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	T tb = this->tx[kx1-1];
	T te = this->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = this->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
        if(arg>this->tx[l]){
			while(arg>this->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		this->bspl(h,this->tx,kx,arg,l);
		l1 = l-kx1;
		sp = 0.0;
		for(i1=1; i1<=kx1; i1++){
			l1 = l1+1;
			sp += h[i1-1]*this->c[l1-1];
		}  
		mod[i-1] = sp;
	}
}

/* Bspl2D class */
template<typename T>
Bspl2D<T>::Bspl2D(const int _nx, const int _nz, const T _dx, const T _dz, const T _dtx, const T _dtz, const int _kx, const int _kz)
{
	/* Setting size values for the evaluation of the spline */
	this->mx=_nx;
	this->mz=_nz;
	this->dtx=_dtx;
	this->dtz=_dtz;
	this->kx=_kx;
	this->kz=_kz;
	
	/* Compute size of knot array */
	double xmax=mx*_dx;
	double xmin=-_dx;
	double zmax=mz*_dz;
	double zmin=-_dz;

	this->ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(this->mx%2==0){
		if(this->ntx%2 != 0){
			this->ntx++;
			this->dtx=(xmax-xmin)/(this->ntx-2*kx-1);
		}
	}else{
		if(this->ntx%2 == 0){
			this->ntx++;
			this->dtx=(xmax-xmin)/(this->ntx-2*kx-1);
		}
	}
	this->ntz=rint((zmax-zmin)/dtz) +(2*kz) + 1;
	/* To preserve symmetry, ntz and mz must be both even or both odd. Thus we need to find the closest dtz that ensures this condition. */
	if(this->mz%2==0){
		if(this->ntz%2 != 0){
			this->ntz++;
			this->dtz=(zmax-zmin)/(this->ntz-2*kz-1);
		}
	}else{
		if(this->ntz%2 == 0){
			this->ntz++;
			this->dtz=(zmax-zmin)/(this->ntz-2*kz-1);
		}
	}

	int i;
	this->wx= (T **) calloc(this->mx,sizeof(T *));
    for(i=0; i<this->mx; i++){
        this->wx[i] = (T *) calloc((this->kx+1), sizeof(T));
    }
	this->wz= (T **) calloc(this->mz,sizeof(T *));
    for(i=0; i<this->mz; i++){
        this->wz[i] = (T *) calloc((this->kz+1), sizeof(T));
    }
	this->nc = (this->ntz-this->kz-1)*(this->ntx-this->kx-1);
	this->c=  (T *) calloc(this->nc, sizeof(T));
	this->tx= (T *) calloc(this->ntx, sizeof(T));
	this->tz= (T *) calloc(this->ntz, sizeof(T));
	this->lx= (int *) calloc(this->mx, sizeof(int));
	this->lz= (int *) calloc(this->mz, sizeof(int));
	this->x= (T *) calloc(this->mx, sizeof(T));
	this->z= (T *) calloc(this->mz, sizeof(T));
	this->mod= (T *) calloc(this->mx*this->mz, sizeof(T));

	/* Computing the evaluation arrays */
	for (i=0; i<this->mx; i++){
		this->x[i]=i*_dx;
	}

	for (i=0; i<this->mz; i++){
		this->z[i]=i*_dz;
	}

	/* Computing the Knotts array */
	for (i=0; i<this->kx+1; i++){
		this->tx[i]=xmin;
	}

	for (i=kx+1; i<this->ntx-this->kx-1; i++){
		this->tx[i]=xmin+(i-this->kx)*this->dtx;
	}

	for (i=this->ntx-this->kx-1; i<this->ntx; i++){
		this->tx[i]=xmax;
	}

	for (i=0; i<this->kz+1; i++){
		this->tz[i]=zmin;
	}

	for (i=kz+1; i<this->ntz-this->kz-1; i++){
		this->tz[i]=zmin+(i-this->kz)*this->dtz;
	}

	for (i=this->ntz-this->kz-1; i<this->ntz; i++){
		this->tz[i]=zmax;
	}
}
 
template<typename T>
Bspl2D<T>::~Bspl2D()
{
	free(c);
	free(tx);
	free(tz);
	free(x);
	free(z);
	free(lx);
	free(lz);
	free(wx[0]);
	free(wx);
	free(wz[0]);
	free(wz);
	free(mod);
}


template<typename T>
void Bspl2D<T>::bisp()
/*<Evaluates the spline coefficients in c over z.>*/
{

	int i,j;
	int i1, j1, l, l1, l2;
	T arg;
	T h[6];

	int nx=this->ntx;
	int kx=this->kx;
	int mx=this->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	T tb = this->tx[kx1-1];
	T te = this->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = this->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>this->tx[l]){
			while(arg>this->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		this->bspl(h,this->tx,kx,arg,l);
		this->lx[i-1] = l-kx1;
		for(j=1; j<=kx1; j++){
			this->wx[i-1][j-1] = h[j-1];
		}
	}

	int nz=this->ntz;
	int kz=this->kz;
	int mz=this->mz;

	int kz1 = kz+1;
	int nkz1 = nz-kz1;
	tb = this->tz[kz1-1];
	te = this->tz[nkz1];
	l=1;
	for (i=1; i<=mz; i++){
		arg = this->z[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>this->tz[l]){
			while(arg>this->tz[l]){
				l++;
			}
		}
		if(l>nkz1) l=nkz1;
		if(l<kz1) l=kz1;
		this->bspl(h,this->tz,kz,arg,l);
		this->lz[i-1] = l-kz1;
		for (j=1; j<=kz1; j++){
			this->wz[i-1][j-1] = h[j-1];
		}	
	}
	int m = 0;
	T sp = 0.0;
	for (i=1; i<=mz; i++){
		l = this->lz[i-1]*nkx1;
		for (i1=1; i1<=kz1; i1++){
			h[i1-1] = this->wz[i-1][i1-1];
		}
		for (j=1; j<=mx; j++){
			l1 = l+this->lx[j-1];
			sp = 0.0;
			for(i1=1; i1<=kz1; i1++){
				l2 = l1;
				for (j1=1;j1<=kx1; j1++){
					l2 = l2+1;
					sp = sp + this->c[l2-1]*h[i1-1]*this->wx[j-1][j1-1];
				}
				l1 = l1+nkx1;
			}  
			m = m+1;
			mod[m-1] = sp;
		}
	}
}

/* Bspl3D class */
template<typename T>
Bspl3D<T>::Bspl3D(const int _nx, const int _ny, const int _nz, const T _dx, const T _dy,const T _dz, const T _dtx, const T _dty, const T _dtz, const int _kx, const int _ky, const int _kz)
{
	/* Setting size values for the evaluation of the spline */
	this->mx=_nx;
	this->my=_ny;
	this->mz=_nz;
	this->dtx=_dtx;
	this->dty=_dty;
	this->dtz=_dtz;
	this->kx=_kx;
	this->ky=_ky;
	this->kz=_kz;
	
	/* Compute size of knot array */
	double xmax=mx*_dx;
	double xmin=-_dx;
	double ymax=my*_dy;
	double ymin=-_dy;
	double zmax=mz*_dz;
	double zmin=-_dz;

	this->ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(this->mx%2==0){
		if(this->ntx%2 != 0){
			this->ntx++;
			this->dtx=(xmax-xmin)/(this->ntx-2*kx-1);
		}
	}else{
		if(this->ntx%2 == 0){
			this->ntx++;
			this->dtx=(xmax-xmin)/(this->ntx-2*kx-1);
		}
	}
	this->nty=rint((ymax-ymin)/dty) +(2*ky) + 1;
	/* To preserve symmetry, nty and my must be both even or both odd. Thus we need to find the closest dty that ensures this condition. */
	if(this->my%2==0){
		if(this->nty%2 != 0){
			this->nty++;
			this->dty=(ymax-ymin)/(this->nty-2*ky-1);
		}
	}else{
		if(this->nty%2 == 0){
			this->nty++;
			this->dty=(ymax-ymin)/(this->nty-2*ky-1);
		}
	}
	this->ntz=rint((zmax-zmin)/dtz) +(2*kz) + 1;
	/* To preserve symmetry, ntz and mz must be both even or both odd. Thus we need to find the closest dtz that ensures this condition. */
	if(this->mz%2==0){
		if(this->ntz%2 != 0){
			this->ntz++;
			this->dtz=(zmax-zmin)/(this->ntz-2*kz-1);
		}
	}else{
		if(this->ntz%2 == 0){
			this->ntz++;
			this->dtz=(zmax-zmin)/(this->ntz-2*kz-1);
		}
	}

	int i;
	this->wx= (T **) calloc(this->mx,sizeof(T *));
    for(i=0; i<this->mx; i++){
        this->wx[i] = (T *) calloc((this->kx+1), sizeof(T));
    }
	this->wy= (T **) calloc(this->my,sizeof(T *));
    for(i=0; i<this->my; i++){
        this->wy[i] = (T *) calloc((this->ky+1), sizeof(T));
    }
	this->wz= (T **) calloc(this->mz,sizeof(T *));
    for(i=0; i<this->mz; i++){
        this->wz[i] = (T *) calloc((this->kz+1), sizeof(T));
    }
	this->nc = (this->ntz-this->kz-1)*(this->ntx-this->kx-1)*(this->nty-this->ky-1);
	this->c=  (T *) calloc(this->nc, sizeof(T));
	this->tx= (T *) calloc(this->ntx, sizeof(T));
	this->ty= (T *) calloc(this->nty, sizeof(T));
	this->tz= (T *) calloc(this->ntz, sizeof(T));
	this->lx= (int *) calloc(this->mx, sizeof(int));
	this->ly= (int *) calloc(this->my, sizeof(int));
	this->lz= (int *) calloc(this->mz, sizeof(int));
	this->x= (T *) calloc(this->mx, sizeof(T));
	this->y= (T *) calloc(this->my, sizeof(T));
	this->z= (T *) calloc(this->mz, sizeof(T));
	this->mod= (T *) calloc(this->mx*this->my*this->mz, sizeof(T));

    /* Computing the evaluation arrays */
	for (i=0; i<this->mx; i++){
		this->x[i]=i*_dx;
	}

	for (i=0; i<this->my; i++){
		this->y[i]=i*_dy;
	}

	for (i=0; i<this->mz; i++){
		this->z[i]=i*_dz;
	}

	/* Computing the Knotts array */
	for (i=0; i<this->kx+1; i++){
		this->tx[i]=xmin;
	}

	for (i=kx+1; i<this->ntx-this->kx-1; i++){
		this->tx[i]=xmin+(i-this->kx)*this->dtx;
	}

	for (i=this->ntx-this->kx-1; i<this->ntx; i++){
		this->tx[i]=xmax;
	}

	for (i=0; i<this->ky+1; i++){
		this->ty[i]=ymin;
	}

	for (i=ky+1; i<this->nty-this->ky-1; i++){
		this->ty[i]=ymin+(i-this->ky)*this->dty;
	}

	for (i=this->nty-this->ky-1; i<this->nty; i++){
		this->ty[i]=ymax;
	}

	for (i=0; i<this->kz+1; i++){
		this->tz[i]=zmin;
	}

	for (i=kz+1; i<this->ntz-this->kz-1; i++){
		this->tz[i]=zmin+(i-this->kz)*this->dtz;
	}

	for (i=this->ntz-this->kz-1; i<this->ntz; i++){
		this->tz[i]=zmax;
	}
}

template<typename T>
Bspl3D<T>::~Bspl3D()
{
	free(this->wx[0]);
	free(this->wx);
	free(this->wz[0]);
	free(this->wz);
	free(this->wy[0]);
	free(this->wy);
	free(this->c);
	free(this->tx);
	free(this->ty);
	free(this->tz);
	free(this->lx);
	free(this->ly);
	free(this->lz);
	free(this->x);
	free(this->y);
	free(this->z);
	free(this->mod);
}

template<typename T>
void Bspl3D<T>::bisp()
/*<Evaluates the spline coefficients in c over z.>*/
{

	int i,j,k;
	int i1, j1, k1, l, l0, l1, l2, l3;
	T arg;
	T h1[6];
	T h2[6];

	int nx=this->ntx;
	int kx=this->kx;
	int mx=this->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	T tb = this->tx[kx1-1];
	T te = this->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = this->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>this->tx[l]){
			while(arg>this->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		this->bspl(h1,this->tx,kx,arg,l);
		this->lx[i-1] = l-kx1;
		for(j=1; j<=kx1; j++){
			this->wx[i-1][j-1] = h1[j-1];
		}
	}

	int ny=this->nty;
	int ky=this->ky;
	int my=this->my;

	int ky1 = ky+1;
	int nky1 = ny-ky1;

	tb = this->ty[ky1-1];
	te = this->ty[nky1];
	l=1;
	for (i=1; i<=my; i++){
		arg = this->y[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>this->ty[l]){
			while(arg>this->ty[l]){
				l++;
			}
		}
		if(l>nky1) l=nky1;
		if(l<ky1)  l=ky1;
		this->bspl(h1,this->ty,ky,arg,l);
		this->ly[i-1] = l-ky1;
		for(j=1; j<=ky1; j++){
			this->wy[i-1][j-1] = h1[j-1];
		}
	}

	int nz=this->ntz;
	int kz=this->kz;
	int mz=this->mz;

	int kz1 = kz+1;
	int nkz1 = nz-kz1;
	tb = this->tz[kz1-1];
	te = this->tz[nkz1];
	l=1;
	for (i=1; i<=mz; i++){
		arg = this->z[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>this->tz[l]){
			while(arg>this->tz[l]){
				l++;
			}
		}
		if(l>nkz1) l=nkz1;
		if(l<kz1) l=kz1;
		this->bspl(h1,this->tz,kz,arg,l);
		this->lz[i-1] = l-kz1;
		for (j=1; j<=kz1; j++){
			this->wz[i-1][j-1] = h1[j-1];
		}	
	}

	int m = 0;
	T sp = 0.0;
	for (i=1; i<=mz; i++){
		l0=this->lz[i-1]*(nky1*nkx1);
		for (i1=1; i1<=kz1; i1++){
			h1[i1-1] = this->wz[i-1][i1-1];
		}
		for (j=1; j<=my; j++){
			l = l0+this->ly[j-1]*nkx1;
			for (i1=1; i1<=ky1; i1++){
				h2[i1-1] = this->wy[j-1][i1-1];
			}
			for (k=1; k<=mx; k++){
				l1 = l+this->lx[k-1];
				sp = 0.0;
				for(i1=1; i1<=kz1; i1++){
					l2 = l1;
					for (j1=1;j1<=ky1; j1++){
						l3 = l2;
						for (k1=1;k1<=kx1; k1++){
							l3 = l3+1;
							sp = sp + this->c[l3-1]*h1[i1-1]*h2[j1-1]*this->wx[k-1][k1-1];
						}
						l2 = l2+nkx1;
					}
					l1=l1+(nky1*nkx1);
				}
				m = m+1;
				mod[m-1] = sp;
			}
		}
	}
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Bspl<float>;
template class Bspl<double>;
template class Bspl1D<float>;
template class Bspl1D<double>;
template class Bspl2D<float>;
template class Bspl2D<double>;
template class Bspl3D<float>;
template class Bspl3D<double>;

}
