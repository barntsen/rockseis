#include "bsplc.h"

bspl_1dspline bspl_1dspline_init(int nx, float dx, double dtx, int kx) 
/*<Initializing the 1D B-spline struct.>*/
{
	/* Creating structure */
	bspl_1dspline spline;
	spline = (bspl_1dspline) malloc(sizeof(*spline));
	if(spline == NULL) {
		fprintf(stderr, "Error in allocation memory! 1D spline struct!\n");
        exit(0);
	}

	/* Setting size values for the evaluation of the spline */
	spline->mx=nx;
	spline->dtx=dtx;
	spline->kx=kx;
	
	/* Compute size of knot array */
	double xmax=nx*dx;
	double xmin=-dx;

	spline->ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(spline->mx%2==0){
		if(spline->ntx%2 != 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}else{
		if(spline->ntx%2 == 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}

	spline->nc = (spline->ntx-spline->kx-1);
	spline->c = (float *) calloc(spline->nc,sizeof(float));
	spline->tx = (float *) calloc(spline->ntx,sizeof(float));
	spline->x = (float *) calloc(spline->mx,sizeof(float));

	/* Computing the evaluation arrays */
	int i;
	for (i=0; i<spline->mx; i++){
		spline->x[i]=i*dx;
	}

	/* Computing the Knotts array */
	for (i=0; i<spline->kx+1; i++){
		spline->tx[i]=xmin;
	}

	for (i=kx+1; i<spline->ntx-spline->kx-1; i++){
		spline->tx[i]=xmin+(i-spline->kx)*spline->dtx;
	}

	for (i=spline->ntx-spline->kx-1; i<spline->ntx; i++){
		spline->tx[i]=xmax;
	}

        /*Initializing c */
	for(i=0; i<spline->nc; i++)
	{
		spline->c[i] = 0.0;
	}

	return(spline);
}

bspl_2dspline bspl_2dspline_init(int nx,int nz, float dx, float dz, double dtx, double dtz, int kx, int kz) 
/*<Initializing the 2D B-spline struct.>*/
{
	/* Creating structure */
	bspl_2dspline spline;
	spline = (bspl_2dspline) malloc(sizeof(*spline));
	if(spline == NULL) {
		fprintf(stderr, "Error in allocation memory! 2D spline struct!\n");
        exit(0);
	}

	/* Setting size values for the evaluation of the spline */
	spline->mx=nx;
	spline->mz=nz;
	spline->dtx=dtx;
	spline->dtz=dtz;
	spline->kx=kx;
	spline->kz=kz;
	
	/* Compute size of knot array */
	double xmax=nx*dx;
	double xmin=-dx;
	double zmax=nz*dz;
	double zmin=-dz;

	spline->ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(spline->mx%2==0){
		if(spline->ntx%2 != 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}else{
		if(spline->ntx%2 == 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}
	spline->ntz=rint((zmax-zmin)/dtz) +(2*kz) + 1;
	/* To preserve symmetry, ntz and mz must be both even or both odd. Thus we need to find the closest dtz that ensures this condition. */
	if(spline->mz%2==0){
		if(spline->ntz%2 != 0){
			spline->ntz++;
			spline->dtz=(zmax-zmin)/(spline->ntz-2*kz-1);
		}
	}else{
		if(spline->ntz%2 == 0){
			spline->ntz++;
			spline->dtz=(zmax-zmin)/(spline->ntz-2*kz-1);
		}
	}

	int i;
	spline->wx= (float **) calloc(spline->mx,sizeof(float));
    spline->wx[0] = (float *) calloc(spline->mx*(spline->kx+1), sizeof(float));
    for(i=1; i<spline->mx; i++){
        spline->wx[i] = spline->wx[0] + i*(spline->kx+1);
    }
	spline->wz= (float **) calloc(spline->mz,sizeof(float));
    spline->wz[0] = (float *) calloc(spline->mz*(spline->kz+1), sizeof(float));
    for(i=1; i<spline->mz; i++){
        spline->wz[i] = spline->wz[0] + i*(spline->kz+1);
    }
	spline->nc = (spline->ntz-spline->kz-1)*(spline->ntx-spline->kx-1);
	spline->c=  (float *) calloc(spline->nc, sizeof(float));
	spline->tx= (float *) calloc(spline->ntx, sizeof(float));
	spline->tz= (float *) calloc(spline->ntz, sizeof(float));
	spline->lx= (int *) calloc(spline->mx, sizeof(int));
	spline->lz= (int *) calloc(spline->mz, sizeof(int));
	spline->x= (float *) calloc(spline->mx, sizeof(float));
	spline->z= (float *) calloc(spline->mz, sizeof(float));

	/* Computing the evaluation arrays */
	for (i=0; i<spline->mx; i++){
		spline->x[i]=i*dx;
	}

	for (i=0; i<spline->mz; i++){
		spline->z[i]=i*dz;
	}

	/* Computing the Knotts array */
	for (i=0; i<spline->kx+1; i++){
		spline->tx[i]=xmin;
	}

	for (i=kx+1; i<spline->ntx-spline->kx-1; i++){
		spline->tx[i]=xmin+(i-spline->kx)*spline->dtx;
	}

	for (i=spline->ntx-spline->kx-1; i<spline->ntx; i++){
		spline->tx[i]=xmax;
	}

	for (i=0; i<spline->kz+1; i++){
		spline->tz[i]=zmin;
	}

	for (i=kz+1; i<spline->ntz-spline->kz-1; i++){
		spline->tz[i]=zmin+(i-spline->kz)*spline->dtz;
	}

	for (i=spline->ntz-spline->kz-1; i<spline->ntz; i++){
		spline->tz[i]=zmax;
	}
        /*Initializing c */
	for(i=0; i<spline->nc; i++)
	{
		spline->c[i] = 0.0;
	}

        /*Initializing wx and wz */
        int j;
	for(i=0; i<spline->mx; i++){
		for(j=0; j<spline->kx+1; j++){
			spline->wx[i][j]=0.0;
		}
	}
	for(i=0; i<spline->mz; i++){
		for(j=0; j<spline->kz+1; j++){
			spline->wz[i][j]=0.0;
		}
	}

	return(spline);
}

bspl_3dspline bspl_3dspline_init(int nx, int ny, int nz, float dx, float dy, float dz, double dtx, double dty, double dtz, int kx, int ky, int kz) 
/*<Initializing the 2D B-spline struct.>*/
{
	/* Creating structure */
	bspl_3dspline spline;
	spline = (bspl_3dspline) malloc(sizeof(*spline));
	if(spline == NULL) {
		fprintf(stderr, "Error in allocation memory! 3D spline struct!\n");
        exit(0);
	}

	/* Setting size values for the evaluation of the spline */
	spline->mx=nx;
	spline->my=ny;
	spline->mz=nz;
	spline->dtx=dtx;
	spline->dty=dty;
	spline->dtz=dtz;
	spline->kx=kx;
	spline->ky=ky;
	spline->kz=kz;
	
	/* Compute size of knot array */
	double xmax=nx*dx;
	double xmin=-dx;
	double ymax=ny*dy;
	double ymin=-dy;
	double zmax=nz*dz;
	double zmin=-dz;

	spline->ntx=rint((xmax-xmin)/dtx) +(2*kx) + 1;
	/* To preserve symmetry, ntx and mx must be both even or both odd. Thus we need to find the closest dtx that ensures this condition. */
	if(spline->mx%2==0){
		if(spline->ntx%2 != 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}else{
		if(spline->ntx%2 == 0){
			spline->ntx++;
			spline->dtx=(xmax-xmin)/(spline->ntx-2*kx-1);
		}
	}
	spline->nty=rint((ymax-ymin)/dty) +(2*ky) + 1;
	/* To preserve symmetry, nty and my must be both even or both odd. Thus we need to find the closest dty that ensures this condition. */
	if(spline->my%2==0){
		if(spline->nty%2 != 0){
			spline->nty++;
			spline->dty=(ymax-ymin)/(spline->nty-2*ky-1);
		}
	}else{
		if(spline->nty%2 == 0){
			spline->nty++;
			spline->dty=(ymax-ymin)/(spline->nty-2*ky-1);
		}
	}
	spline->ntz=rint((zmax-zmin)/dtz) +(2*kz) + 1;
	/* To preserve symmetry, ntz and mz must be both even or both odd. Thus we need to find the closest dtz that ensures this condition. */
	if(spline->mz%2==0){
		if(spline->ntz%2 != 0){
			spline->ntz++;
			spline->dtz=(zmax-zmin)/(spline->ntz-2*kz-1);
		}
	}else{
		if(spline->ntz%2 == 0){
			spline->ntz++;
			spline->dtz=(zmax-zmin)/(spline->ntz-2*kz-1);
		}
	}

	int i;
	spline->wx= (float **) calloc(spline->mx,sizeof(float));
    spline->wx[0] = (float *) calloc(spline->mx*(spline->kx+1), sizeof(float));
    for(i=1; i<spline->mx; i++){
        spline->wx[i] = spline->wx[0] + i*(spline->kx+1);
    }
	spline->wy= (float **) calloc(spline->my,sizeof(float));
    spline->wy[0] = (float *) calloc(spline->my*(spline->ky+1), sizeof(float));
    for(i=1; i<spline->my; i++){
        spline->wy[i] = spline->wy[0] + i*(spline->ky+1);
    }
	spline->wz= (float **) calloc(spline->mz,sizeof(float));
    spline->wz[0] = (float *) calloc(spline->mz*(spline->kz+1), sizeof(float));
    for(i=1; i<spline->mz; i++){
        spline->wz[i] = spline->wz[0] + i*(spline->kz+1);
    }
	spline->nc = (spline->ntz-spline->kz-1)*(spline->ntx-spline->kx-1)*(spline->nty-spline->ky-1);
	spline->c=  (float *) calloc(spline->nc, sizeof(float));
	spline->tx= (float *) calloc(spline->ntx, sizeof(float));
	spline->ty= (float *) calloc(spline->nty, sizeof(float));
	spline->tz= (float *) calloc(spline->ntz, sizeof(float));
	spline->lx= (int *) calloc(spline->mx, sizeof(int));
	spline->ly= (int *) calloc(spline->my, sizeof(int));
	spline->lz= (int *) calloc(spline->mz, sizeof(int));
	spline->x= (float *) calloc(spline->mx, sizeof(float));
	spline->y= (float *) calloc(spline->my, sizeof(float));
	spline->z= (float *) calloc(spline->mz, sizeof(float));

    /* Computing the evaluation arrays */
	for (i=0; i<spline->mx; i++){
		spline->x[i]=i*dx;
	}

	for (i=0; i<spline->my; i++){
		spline->y[i]=i*dy;
	}

	for (i=0; i<spline->mz; i++){
		spline->z[i]=i*dz;
	}

	/* Computing the Knotts array */
	for (i=0; i<spline->kx+1; i++){
		spline->tx[i]=xmin;
	}

	for (i=kx+1; i<spline->ntx-spline->kx-1; i++){
		spline->tx[i]=xmin+(i-spline->kx)*spline->dtx;
	}

	for (i=spline->ntx-spline->kx-1; i<spline->ntx; i++){
		spline->tx[i]=xmax;
	}

	for (i=0; i<spline->ky+1; i++){
		spline->ty[i]=ymin;
	}

	for (i=ky+1; i<spline->nty-spline->ky-1; i++){
		spline->ty[i]=ymin+(i-spline->ky)*spline->dty;
	}

	for (i=spline->nty-spline->ky-1; i<spline->nty; i++){
		spline->ty[i]=ymax;
	}

	for (i=0; i<spline->kz+1; i++){
		spline->tz[i]=zmin;
	}

	for (i=kz+1; i<spline->ntz-spline->kz-1; i++){
		spline->tz[i]=zmin+(i-spline->kz)*spline->dtz;
	}

	for (i=spline->ntz-spline->kz-1; i<spline->ntz; i++){
		spline->tz[i]=zmax;
	}
        /*Initializing c */
	for(i=0; i<spline->nc; i++)
	{
		spline->c[i] = 0.0;
	}

        /*Initializing wx and wz */
        int j;
	for(i=0; i<spline->mx; i++){
		for(j=0; j<spline->kx+1; j++){
			spline->wx[i][j]=0.0;
		}
	}
	for(i=0; i<spline->my; i++){
		for(j=0; j<spline->ky+1; j++){
			spline->wy[i][j]=0.0;
		}
	}
	for(i=0; i<spline->mz; i++){
		for(j=0; j<spline->kz+1; j++){
			spline->wz[i][j]=0.0;
		}
	}

	return(spline);
}

void bspl_1dspline_free(bspl_1dspline spline)
/*<Free bspl_1dspline structure from memory.>*/
{
	free(spline->c);
	free(spline->tx);
	free(spline->x);
	free(spline);
}

void bspl_2dspline_free(bspl_2dspline spline)
/*<Free bspl_2dspline structure from memory.>*/
{
	free(spline->wx[0]);
	free(spline->wx);
	free(spline->wz[0]);
	free(spline->wz);
	free(spline->c);
	free(spline->tx);
	free(spline->tz);
	free(spline->lx);
	free(spline->lz);
	free(spline->x);
	free(spline->z);
	free(spline);
}

void bspl_3dspline_free(bspl_3dspline spline)
/*<Free bspl_3dspline structure from memory.>*/
{
	free(spline->wx[0]);
	free(spline->wx);
	free(spline->wz[0]);
	free(spline->wz);
	free(spline->wy[0]);
	free(spline->wy);
	free(spline->c);
	free(spline->tx);
	free(spline->ty);
	free(spline->tz);
	free(spline->lx);
	free(spline->ly);
	free(spline->lz);
	free(spline->x);
	free(spline->y);
	free(spline->z);
	free(spline);
}

void bspl(float *h, float *t, int k, float x, int l)
/*<evaluates the (k+1) non-zero b-splines of degree k at t(l) <= x < t(l+1) using the stable recurrence relation of de Boor and Cox.>*/
{
	float hh[6];
	int li, lj;
	float f;
	float one = 1.0; /* 0.1e+01; */ 
	
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


void bspl_bisp1d(float *mod, bspl_1dspline spline)
/*<Evaluates the spline coefficients in spline->c over z.>*/
{

	int i;
	int i1, l, l1;
	float arg;
	float h[6];
    for(i1=0; i1<6; i1++) h[1]=0.0;
	float sp = 0.0;

	int nx=spline->ntx;
	int kx=spline->kx;
	int mx=spline->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	float tb = spline->tx[kx1-1];
	float te = spline->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = spline->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
        if(arg>spline->tx[l]){
			while(arg>spline->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		bspl(h,spline->tx,kx,arg,l);
		l1 = l-kx1;
		sp = 0.0;
		for(i1=1; i1<=kx1; i1++){
			l1 = l1+1;
			sp += h[i1-1]*spline->c[l1-1];
		}  
		mod[i-1] = sp;
	}
}

void bspl_bisp2d(float *mod, bspl_2dspline spline)
	/*<Evaluates the spline coefficients in spline->c over z.>*/
{

	int i,j;
	int i1, j1, l, l1, l2;
	float arg;
	float h[6];

	int nx=spline->ntx;
	int kx=spline->kx;
	int mx=spline->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	float tb = spline->tx[kx1-1];
	float te = spline->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = spline->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>spline->tx[l]){
			while(arg>spline->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		bspl(h,spline->tx,kx,arg,l);
		spline->lx[i-1] = l-kx1;
		for(j=1; j<=kx1; j++){
			spline->wx[i-1][j-1] = h[j-1];
		}
	}

	int nz=spline->ntz;
	int kz=spline->kz;
	int mz=spline->mz;

	int kz1 = kz+1;
	int nkz1 = nz-kz1;
	tb = spline->tz[kz1-1];
	te = spline->tz[nkz1];
	l=1;
	for (i=1; i<=mz; i++){
		arg = spline->z[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>spline->tz[l]){
			while(arg>spline->tz[l]){
				l++;
			}
		}
		if(l>nkz1) l=nkz1;
		if(l<kz1) l=kz1;
		bspl(h,spline->tz,kz,arg,l);
		spline->lz[i-1] = l-kz1;
		for (j=1; j<=kz1; j++){
			spline->wz[i-1][j-1] = h[j-1];
		}	
	}
	int m = 0;
	float sp = 0.0;
	for (i=1; i<=mz; i++){
		l = spline->lz[i-1]*nkx1;
		for (i1=1; i1<=kz1; i1++){
			h[i1-1] = spline->wz[i-1][i1-1];
		}
		for (j=1; j<=mx; j++){
			l1 = l+spline->lx[j-1];
			sp = 0.0;
			for(i1=1; i1<=kz1; i1++){
				l2 = l1;
				for (j1=1;j1<=kx1; j1++){
					l2 = l2+1;
					sp = sp + spline->c[l2-1]*h[i1-1]*spline->wx[j-1][j1-1];
				}
				l1 = l1+nkx1;
			}  
			m = m+1;
			mod[m-1] = sp;
		}
	}
}

void bspl_bisp3d(float *mod, bspl_3dspline spline)
	/*<Evaluates the spline coefficients in spline->c over z.>*/
{

	int i,j,k;
	int i1, j1, k1, l, l0, l1, l2, l3;
	float arg;
	float h1[6];
	float h2[6];

	int nx=spline->ntx;
	int kx=spline->kx;
	int mx=spline->mx;

	int kx1 = kx+1;
	int nkx1 = nx-kx1;

	float tb = spline->tx[kx1-1];
	float te = spline->tx[nkx1];
	l=1;
	for (i=1; i<=mx; i++){
		arg = spline->x[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>spline->tx[l]){
			while(arg>spline->tx[l]){
				l++;
			}
		}
		if(l>nkx1) l=nkx1;
		if(l<kx1)  l=kx1;
		bspl(h1,spline->tx,kx,arg,l);
		spline->lx[i-1] = l-kx1;
		for(j=1; j<=kx1; j++){
			spline->wx[i-1][j-1] = h1[j-1];
		}
	}

	int ny=spline->nty;
	int ky=spline->ky;
	int my=spline->my;

	int ky1 = ky+1;
	int nky1 = ny-ky1;

	tb = spline->ty[ky1-1];
	te = spline->ty[nky1];
	l=1;
	for (i=1; i<=my; i++){
		arg = spline->y[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>spline->ty[l]){
			while(arg>spline->ty[l]){
				l++;
			}
		}
		if(l>nky1) l=nky1;
		if(l<ky1)  l=ky1;
		bspl(h1,spline->ty,ky,arg,l);
		spline->ly[i-1] = l-ky1;
		for(j=1; j<=ky1; j++){
			spline->wy[i-1][j-1] = h1[j-1];
		}
	}

	int nz=spline->ntz;
	int kz=spline->kz;
	int mz=spline->mz;

	int kz1 = kz+1;
	int nkz1 = nz-kz1;
	tb = spline->tz[kz1-1];
	te = spline->tz[nkz1];
	l=1;
	for (i=1; i<=mz; i++){
		arg = spline->z[i-1];
		if(arg<tb) arg = tb;
		if(arg>te) arg = te;
		if(arg>spline->tz[l]){
			while(arg>spline->tz[l]){
				l++;
			}
		}
		if(l>nkz1) l=nkz1;
		if(l<kz1) l=kz1;
		bspl(h1,spline->tz,kz,arg,l);
		spline->lz[i-1] = l-kz1;
		for (j=1; j<=kz1; j++){
			spline->wz[i-1][j-1] = h1[j-1];
		}	
	}

	int m = 0;
	float sp = 0.0;
	for (i=1; i<=mz; i++){
		l0=spline->lz[i-1]*(nky1*nkx1);
		for (i1=1; i1<=kz1; i1++){
			h1[i1-1] = spline->wz[i-1][i1-1];
		}
		for (j=1; j<=my; j++){
			l = l0+spline->ly[j-1]*nkx1;
			for (i1=1; i1<=ky1; i1++){
				h2[i1-1] = spline->wy[j-1][i1-1];
			}
			for (k=1; k<=mx; k++){
				l1 = l+spline->lx[k-1];
				sp = 0.0;
				for(i1=1; i1<=kz1; i1++){
					l2 = l1;
					for (j1=1;j1<=ky1; j1++){
						l3 = l2;
						for (k1=1;k1<=kx1; k1++){
							l3 = l3+1;
							sp = sp + spline->c[l3-1]*h1[i1-1]*h2[j1-1]*spline->wx[k-1][k1-1];
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

