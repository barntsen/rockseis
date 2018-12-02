// Include statements
#include "lsrtm.h"

namespace rockseis {

    // Elastic: Derivate V Source gradients inside scaleGrad
// =============== ABSTRACT LSRTM CLASS =============== //
template<typename T>
Lsrtm<T>::Lsrtm() {
	order = 4;
    snapinc=1;
    snapmethod = FULL;
    ncheck = 0;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 0;
    misfit_type = DIFFERENCE;
    misfit = 0.0;
    noreverse = false;
    filter = false;
    freqs[0]=0;
    freqs[1]=0;
    freqs[2]=0;
    freqs[3]=0;
}

template<typename T>
Lsrtm<T>::Lsrtm(int _order, int _snapinc) {
    if(_order > 1 && _order < 9)
    {
        order = _order;
    }else{
        order = 4;
    }
    if(_snapinc > 0)
    {
        snapinc = _snapinc;
    }else{
        snapinc=1;
    }

    snapmethod = FULL;
    ncheck = 0;
	prog.previous = 0;
	prog.current = 0;
    prog.persec = 1;
    misfit_type = DIFFERENCE;
    misfit = 0.0;
    noreverse = false;
    filter = false;
    freqs[0]=0;
    freqs[1]=0;
    freqs[2]=0;
    freqs[3]=0;
}

template<typename T>
bool Lsrtm<T>::createLog(std::string filename){
	logfile = filename;
	Flog.open(logfile.c_str());
	if(Flog.fail()){
		Flog.close();
		return LSRTM_ERR;
	}else{
		Flog.close();
		return LSRTM_OK;
	}
}

template<typename T>
void Lsrtm<T>::writeLog(std::string text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Lsrtm<T>::writeLog(const char *text){
    if(!logfile.empty()){
        Flog.open(logfile.c_str(), std::ios::app);
        if(!Flog.fail())
            Flog << text << std::endl;
        Flog.close();
    }
}

template<typename T>
void Lsrtm<T>::writeProgressbar(int x, int n, int r, int w){
	if(!logfile.empty()){
		if ( n <= 0 ) n = 1;
		if ( r > n ) r = n;
		if ( w > 48 ) w = 48;
		// Only update r times.
		if (  x % (n/r) != 0 ) return;

		// Calculuate the ratio of complete-to-incomplete.
		float ratio = x/(float)n;
		int   c     = ratio * w;
		prog.current = clock();
		float time_spent  = (float) ( prog.current - prog.previous ) / CLOCKS_PER_SEC;
		prog.previous = prog.current;
		prog.persec = 0.0;
		if(time_spent > 0) prog.persec = (n/r) / time_spent;
		snprintf(prog.speed, 48, "%.2f its/s", prog.persec); 

		// Show the percentage complete.
		sprintf(prog.progress,"%3d%% [", (int)(ratio*100) );

		// Show the load bar.
		for (x=0; x<c; x++)
			strcat(prog.progress,"=");

		for (x=c; x<w; x++)
			strcat(prog.progress," ");

		strcat(prog.progress,"]"); 
		strcat(prog.progress,"  "); 
		strcat(prog.progress, prog.speed); 
		writeLog(prog.progress); // Write to file
	}
}

template<typename T>
void Lsrtm<T>::writeProgress(int x, int n, int r, int w){
	if(!logfile.empty()){
		if ( n <= 0 ) n = 1;
		if ( r > n ) r = n;
		if ( w > 48 ) w = 48;
		// Only update r times.
		if (  x % (n/r) != 0 ) return;

		// Calculuate the ratio of complete-to-incomplete.
		float ratio = x/(float)n;
		prog.current = clock();
		float time_spent  = (float) ( prog.current - prog.previous ) / CLOCKS_PER_SEC;
		prog.previous = prog.current;
		prog.persec = 0.0;
		if(time_spent > 0) prog.persec = (n/r) / time_spent;
		snprintf(prog.speed, 48, "%.2f its/s", prog.persec); 

		// Show the percentage complete.
		sprintf(prog.progress,"%3d%%", (int)(ratio*100) );
		strcat(prog.progress,"  "); 
		strcat(prog.progress, prog.speed); 
		writeLog(prog.progress); // Write to file
	}
}

template<typename T>
void Lsrtm<T>::stoep (int n, T r[], T g[], T f[], T a[])
	/*<Solve a symmetric Toeplitz linear system of equations Rf=g.>*/
{
	int i,j;
	T v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {

		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}

template<typename T>
void Lsrtm<T>::convolve(int lx, int ifx, T *x, int ly, int ify, T *y, int lz, int ifz, T *z)
/*<Compute z = x convolved with y>*/
{
	int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,
		i,j,ilow,ihigh,jlow,jhigh;
	T sa,sb,xa,xb,ya,yb,*t;

	/* if x is longer than y, swap x and y */
	if (lx>ly) {
		i = ifx;  ifx = ify;  ify = i;
		i = ilx;  ilx = ily;  ily = i;
		i = lx;  lx = ly;  ly = i;
		t = x;  x = y;  y = t;
	}
	
	/* adjust pointers for indices of first samples */
	x -= ifx;
	y -= ify;
	z -= ifz;
		
	/* OFF LEFT:  i < ify+ifx */
	
	/* zero output for all i */
	ilow = ifz;
	ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;

	/* ROLLING ON:  ify+ifx <= i < ify+ilx */
	
	/* if necessary, do one i so that number of j in overlap is odd */
	if (i<ify+ilx && i<=ilz) {
		jlow = ifx;
		jhigh = i-ify;
		if ((jhigh-jlow)%2) {
			sa = 0.0;
			for (j=jlow; j<=jhigh; ++j)
				sa += x[j]*y[i-j];
			z[i++] = sa;
		}
	}
	
	/* loop over pairs of i and j */
	ilow = i;
	ihigh = ilx+ify-1;  if (ihigh>ilz) ihigh = ilz;
	jlow = ifx;
	jhigh = ilow-ify;
	for (i=ilow; i<ihigh; i+=2,jhigh+=2) {
		sa = sb = 0.0;
		xb = x[jhigh+1];
		yb = 0.0;
		for (j=jhigh; j>=jlow; j-=2) {
			sa += xb*yb;
			ya = y[i-j];
			sb += xb*ya;
			xa = x[j];
			sa += xa*ya;
			yb = y[i+1-j];
			sb += xa*yb;
			xb = x[j-1];
		}
		z[i] = sa;
		z[i+1] = sb;
	}
	
	/* if number of i is odd */
	if (i==ihigh) {
		jlow = ifx;
		jhigh = i-ify;
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}
	
	/* MIDDLE:  ify+ilx <= i <= ily+ifx */
	
	/* determine limits for i and j */
	ilow = i;
	ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
	jlow = ifx;
	jhigh = ilx;
	
	/* if number of j is even, do j in pairs with no leftover */
	if ((jhigh-jlow)%2) {
		for (i=ilow; i<ihigh; i+=2) {
			sa = sb = 0.0;
			yb = y[i+1-jlow];
			xa = x[jlow];
			for (j=jlow; j<jhigh; j+=2) {
				sb += xa*yb;
				ya = y[i-j];
				sa += xa*ya;
				xb = x[j+1];
				sb += xb*ya;
				yb = y[i-1-j];
				sa += xb*yb;
				xa = x[j+2];
			}
			z[i] = sa;
			z[i+1] = sb;
		}
	
	/* else, number of j is odd, so do j in pairs with leftover */
	} else {
		for (i=ilow; i<ihigh; i+=2) {
			sa = sb = 0.0;
			yb = y[i+1-jlow];
			xa = x[jlow];
			for (j=jlow; j<jhigh; j+=2) {
				sb += xa*yb;
				ya = y[i-j];
				sa += xa*ya;
				xb = x[j+1];
				sb += xb*ya;
				yb = y[i-1-j];
				sa += xb*yb;
				xa = x[j+2];
			}
			z[i] = sa+x[jhigh]*y[i-jhigh];
			z[i+1] = sb+x[jhigh]*y[i+1-jhigh];
		}
	}
	
	/* if number of i is odd */
	if (i==ihigh) {
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}

	/* ROLLING OFF:  ily+ifx < i <= ily+ilx */
	
	/* if necessary, do one i so that number of j in overlap is even */
	if (i<=ily+ilx && i<=ilz) {
		jlow = i-ily;
		jhigh = ilx;
		if (!((jhigh-jlow)%2)) {
			sa = 0.0;
			for (j=jlow; j<=jhigh; ++j)
				sa += x[j]*y[i-j];
			z[i++] = sa;
		}
	}
	
	/* number of j is now even, so loop over both i and j in pairs */
	ilow = i;
	ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
	jlow = ilow-ily;
	jhigh = ilx-2; /* Dave's new patch */
        for (i=ilow; i<ihigh; i+=2,jlow+=2) {
                sa = sb = 0.0;
                xa = x[jlow];
                yb = 0.0;
                for (j=jlow; j<jhigh; j+=2) {
                        sb += xa*yb;
                        ya = y[i-j];
                        sa += xa*ya;
                        xb = x[j+1];
                        sb += xb*ya;
                        yb = y[i-1-j];
                        sa += xb*yb;
                        xa = x[j+2];
                }
                sb += xa*yb;
                ya = y[i-j];
                sa += xa*ya;
                xb = x[j+1];
                sb += xb*ya;
                yb = y[i-1-j];
                sa += xb*yb;
                z[i] = sa;
                z[i+1] = sb;
        }
	
	/* if number of i is odd */
	if (i==ihigh) {
		jlow = i-ily;
		jhigh = ilx;
		sa = 0.0;
		for (j=jlow; j<=jhigh; ++j)
			sa += x[j]*y[i-j];
		z[i++] = sa;
	}
	
	/* OFF RIGHT:  ily+ilx < i */
	
	/* zero output for all i */
	ilow = i;
	ihigh = ilz;
	for (i=ilow; i<=ihigh; ++i)
		z[i] = 0.0;
}

template<typename T>
void Lsrtm<T>::xcor (int lx, int ifx, T *x,int ly, int ify, T *y, int lz, int ifz, T *z)
	/*< Compute z = x cross-correlated with y>*/
{
	int i,j;
	T *xr;

	xr = (T *) calloc(lx, sizeof(T));
	for (i=0,j=lx-1; i<lx; ++i,--j)
		xr[i] = x[j];
	convolve(lx,1-ifx-lx,xr,ly,ify,y,lz,ifz,z);
	free(xr);
}

template<typename T>
T Lsrtm<T>::gauss (int it, T stdev)
{
    if(stdev == 0.0) return 0;
    T val; 
    val = exp(-1.0*(it*it)/(2.0*(stdev*stdev)));
    return val;
}

template<typename T>
T Lsrtm<T>::linear (int it, T stdev)
{
    if(stdev == 0.0) return 0;
    
    if(it < stdev){
        return it;
    }else{
        return stdev;
    }
}


template<typename T>
void Lsrtm<T>::apply_filter (T *data, int nt, T dt)
{
    int i;
    double d_dt = dt;
    float f[4];
    f[0] = freqs[0];
    f[1] = freqs[1];
    f[2] = freqs[2];
    f[3] = freqs[3];
	float *wrk = (float *) calloc(nt, sizeof(float));
	float *flt = (float *) calloc(nt, sizeof(float));
    for(i=0; i< nt; i++) flt[i] = data[i];
    sig_filt(flt, f[0], f[1], f[2], f[3], d_dt, nt); 
    for(i=0; i< nt; i++) data[i] = flt[i];

    free(wrk);
    free(flt);
}

template<typename T>
Lsrtm<T>::~Lsrtm() {
    // Nothing here
}

// =============== ACOUSTIC 2D LSRTM CLASS =============== //

template<typename T>
LsrtmAcoustic2D<T>::LsrtmAcoustic2D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    vpgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
    srcilumset = false;
}

template<typename T>
LsrtmAcoustic2D<T>::LsrtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> _model, std::shared_ptr<Image2D<T>> _pimage, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataP, int order, int snapinc):Lsrtm<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    pimage = _pimage;
    sourceset = true;
    dataPset = true;
    modelset = true;
    vpgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
    srcilumset = false;
}

template<typename T>
void LsrtmAcoustic2D<T>::calcAdjointsource(T *adjsrc_bw, T* wrp, int padr)
{
	if(!pimage->getAllocated()) rs_error("LsrtmAcoustic2D<T>::LsrtmAcoustic2D<T>::calcAdjointsource: pimage is not allocated.");
	int ix, iz, ihx, ihz;
	T *imagedata = pimage->getImagedata();
	int nhx = pimage->getNhx();
	int nhz = pimage->getNhz();
	int nx = pimage->getNx();
	int nxr = nx+2*padr;
	int nz = pimage->getNz();
	int hx, hz;
    //Reset arrays
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            adjsrc_bw[km2D(ix,iz)] = 0.0;
        }
    }
    // Calculate integral
    for (ihx=0; ihx<nhx; ihx++){
        hx= -(nhx-1)/2 + ihx;
        for (ihz=0; ihz<nhz; ihz++){
            hz= -(nhz-1)/2 + ihz;
            for (ix=0; ix<nx; ix++){
                if( ((ix-2*hx) >= 0) && ((ix-2*hx) < nx) && ((ix+2*hx) >= 0) && ((ix+2*hx) < nx))
                {
                    for (iz=0; iz<nz; iz++){
                        if( ((iz-2*hz) >= 0) && ((iz-2*hz) < nz) && ((iz+2*hz) >= 0) && ((iz+2*hz) < nz))
                        {
							adjsrc_bw[km2D(ix,iz)] -= imagedata[kim2D(ix-hx,iz-hz,ihx,ihz)]*wrp[kr2D(ix-2*hx+padr, iz-2*hz+padr)];

                        }

                    }	
                }
            }
        }
    }
}

template<typename T>
void LsrtmAcoustic2D<T>::insertAdjointsource(std::shared_ptr<WavesAcoustic2D<T>> waves_bw, T* adjsrc_bw, T* L)
{
    int ix, iz;
    int nx = pimage->getNx();
    int nz = pimage->getNz();
    T dt = waves_bw->getDt();
    int padw = waves_bw->getLpml();
    int nxw = nx+2*padw;
    T* wr = waves_bw->getP2();
    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            wr[kw2D(ix+padw,iz+padw)] += dt*dt*L[kw2D(ix+padw,iz+padw)]*adjsrc_bw[km2D(ix, iz)];

        }
    }	
}


template<typename T>
void LsrtmAcoustic2D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T* Vp, T* Rho)
{
    if(!vpgradset) rs_error("LsrtmAcoustic2D:crossCorr: No gradient set in lsrtm class");
    if(!vpgrad->getAllocated()) vpgrad->allocateImage();
    if(srcilumset){
        if(!srcilum->getAllocated()) srcilum->allocateImage();
    }

    int ix, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T *srcilumdata;
    if(srcilumset){
        srcilumdata = srcilum->getImagedata();
    }
    T mrxx, mrzz;
    T L;
    int nx = vpgrad->getNx();
    int nz = vpgrad->getNz();
    int dx = vpgrad->getDx();
    int dz = vpgrad->getDz();
    int nxs = nx+2*pads;
    int nxr = nx+2*padr;
	//int nhx = vpgrad->getNhx();
	//int nhz = vpgrad->getNhz();
    for (ix=1; ix<nx-1; ix++){
        {
            for (iz=1; iz<nz-1; iz++){
                L = Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)];
                mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
                mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;
                vpgraddata[ki2D(ix,iz)] += wsp[ks2D(ix+pads, iz+pads)]*L*(mrxx + mrzz);

                if(srcilumset){
                    srcilumdata[ki2D(ix,iz)] += wsp[ks2D(ix+pads, iz+pads)]*wsp[ks2D(ix+pads, iz+pads)];
                }
            }	
        }
    }
}

template<typename T>
void LsrtmAcoustic2D<T>::computeMisfit(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    T res = 0.0;
    T misfit = 0.0;
    T *shaper;
    T *spiker;
    T *autocorr;
    T *crosscorr;
    T *modwei, *recwei;
    T norm; 
    T H;
    T stdev;

    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;

                for(it=0; it<nt; it++){
                    res = (-1.0)*(mod[I(it, itr)]*rec[I(it, itr)]/(norm1*norm2));
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += res;
                }
            }

            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);

                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                if (norm == 0.0) norm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it,stdev);
                    res = H*shaper[it];
                    misfit -= 0.5*res*res/norm;
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);

                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                if (norm == 0.0) norm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it,stdev);
                    res = H*shaper[it];
                    misfit += 0.5*res*res/norm;
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
    }
    // Set the final misfit value
    this->setMisfit(misfit);
}

template<typename T>
void LsrtmAcoustic2D<T>::computeResiduals(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodP->getNt();
    T dt = datamodP->getDt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* res = dataresP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }

    T *shaper;
    T *spiker;
    T *autocorr;
    T *crosscorr;
    T *modwei, *recwei;
    T *fres;
    T pclip; 
    int pos;
    T norm; 
    T misfit;
    T H;
    T stdev;

    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2, norm3;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                norm3 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                   norm3 += mod[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;
                norm3 /= (norm1*norm2);

                for(it=0; it<nt; it++){
                    res[I(it, itr)]=(-1.0)*((rec[I(it, itr)]/(norm1*norm2)) - (mod[I(it, itr)]/(norm1*norm1))*norm3);
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(res[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(res[I(it, itr)]) >= pclip){
                        res[I(it, itr)] = 0.0;
                    }
                }
            }
            free(fres);

            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                // If zero trace we need to add something to avoid division by zero
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);
                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                
                // Ensure norm different than 0.0
                if (norm == 0.0) norm = 1.0;
                misfit = 0.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    misfit += 0.5*(H*shaper[it])*(H*shaper[it])/norm;
                }

                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    res[I(it, itr)] = -1.0*(H*H - 2.0*misfit)*shaper[it]/norm;
                }
                this->stoep(nt, autocorr, &res[I(0, itr)], shaper, spiker);
                this->convolve(nt, 0, shaper, nt, 0, recwei, nt, 0, &res[I(0, itr)]);        
                if(this->getFilter()){
                    this->apply_filter(&res[I(0, itr)], nt, dt);
                }
                if(dataweightset){
                    for (it=0; it < nt; it++)
                    {
                        res[I(it,itr)] *= wei[I(it,itr)];
                    }
                }
            }
            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(res[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(res[I(it, itr)]) >= pclip){
                        res[I(it, itr)] = 0.0;
                    }
                }
            }
            free(fres);
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                // If zero trace we need to add something to avoid division by zero
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);
                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                
                // Ensure norm different than 0.0
                if (norm == 0.0) norm = 1.0;
                misfit = 0.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    misfit += 0.5*(H*shaper[it])*(H*shaper[it])/norm;
                }

                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    res[I(it, itr)] = (H*H - 2.0*misfit)*shaper[it]/norm;
                }
                this->stoep(nt, autocorr, &res[I(0, itr)], shaper, spiker);
                this->convolve(nt, 0, shaper, nt, 0, recwei, nt, 0, &res[I(0, itr)]);        
                if(this->getFilter()){
                    this->apply_filter(&res[I(0, itr)], nt, dt);
                }
                if(dataweightset){
                    for (it=0; it < nt; it++)
                    {
                        res[I(it,itr)] *= wei[I(it,itr)];
                    }
                }
            }
            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(res[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(res[I(it, itr)]) >= pclip){
                        res[I(it, itr)] = 0.0;
                    }
                }
            }
            free(fres);
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int LsrtmAcoustic2D<T>::run(){
     int result = LSRTM_ERR;
     int nt;
     float dt;
	 float ot;
     T *wsp;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();
     if(!this->datamodPset || !this->dataresPset) rs_error("LsrtmAcoustic2D::run: datamodP and dataresP must be set before running the simulation.");

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_adj (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Psnap;
     Psnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

     T* adjsrc;
     adjsrc = (T *) calloc(waves_adj->getNx()*waves_adj->getNz(), sizeof(T));

     this->writeLog("Running 2D Acoustic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);

    	waves_adj->forwardstepAcceleration(model, der);
    	waves_adj->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Inserting adjoint sources
        wsp = waves->getP1();
        this->calcAdjointsource(adjsrc, wsp, waves_adj->getLpml());
        this->insertAdjointsource(waves_adj, adjsrc, model->getL());

        // Recording data 
        if(this->datamodPset){
            waves_adj->recordData(this->datamodP, GMAP, it);
        }

    	//Writting out results to snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field
        Psnap->writeSnap(it);

    	// Roll the pointers P1 and P2
    	waves->roll();
    	waves_adj->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Psnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesAcoustic2D<T>>(model, nt, dt, ot);

    // Create images 
    vpgrad->allocateImage();

    // Compute misfit
    computeMisfit();

    if(!this->getNoreverse()){
        // Compute residuals
        computeResiduals();

        Psnap->openSnap(this->getSnapfile(), 'r');
        Psnap->allocSnap(0);

        this->writeLog("\nDoing reverse-time Loop.");
        // Loop over reverse time
        for(int it=0; it < nt; it++)
        {
            // Time stepping
            waves->forwardstepAcceleration(model, der);
            waves->forwardstepStress(model, der);

            // Inserting residuals
            waves->insertSource(model, dataresP, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Psnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
                T *wrp = waves->getP1();
                T* wrx = waves->getAx(); 
                T* wrz = waves->getAz(); 
                crossCorr(Psnap->getData(0), 0, wrp, wrx, wrz, waves->getLpml(), model->getVp(), model->getR());
            }

            // Roll the pointers P1 and P2
            waves->roll();

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");

    } // End of reverse loop

    //Remove snapshot file
    Psnap->removeSnap();

    // Free arrays
    free(adjsrc);

    result=LSRTM_OK;
    return result;
}

template<typename T>
int LsrtmAcoustic2D<T>::run_optimal(){
     int result = LSRTM_ERR;
     int nt;
     float dt;
	 float ot;
     T *wsp;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic2D<T>> waves_fw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_bw (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic2D<T>> waves_adj (new WavesAcoustic2D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), 1, waves_fw->getNz_pml(), waves_fw->getDx(), 1.0, waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
     vpgrad->allocateImage();

    // Allocate memory for adjoint sources 
    T* adjsrc;
    adjsrc = (T *) calloc(waves_adj->getNx()*waves_adj->getNz(), sizeof(T));

     this->writeLog("Running 2D Acoustic full-waveform inversion gradient with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal->getCapo();
        whatodo = optimal->revolve();
        capo = optimal->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping
                waves_fw->forwardstepAcceleration(model, der);
                waves_fw->forwardstepStress(model, der);



                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);


                // Recording data 
                if(this->datamodPset && !reverse){
                    waves_adj->forwardstepAcceleration(model, der);
                    waves_adj->forwardstepStress(model, der);
                    // Inserting adjoint sources
                    wsp = waves_fw->getP1();
                    this->calcAdjointsource(adjsrc, wsp, waves_adj->getLpml());
                    this->insertAdjointsource(waves_adj, adjsrc, model->getL());
                    waves_adj->recordData(this->datamodP, GMAP, it);
                    waves_adj->roll();
                }

                // Roll the pointers P1 and P2
                waves_fw->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepAcceleration(model, der);
            waves_fw->forwardstepStress(model, der);

            waves_adj->forwardstepAcceleration(model, der);
            waves_adj->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Inserting adjoint sources
            wsp = waves_fw->getP1();
            this->calcAdjointsource(adjsrc, wsp, waves_adj->getLpml());
            this->insertAdjointsource(waves_adj, adjsrc, model->getL());
            waves_adj->recordData(this->datamodP, GMAP, capo);

            // Recording data 
            if(this->datamodPset){
                waves_adj->recordData(this->datamodP, GMAP, capo);
            }

            // Compute misfit
            computeMisfit();

            // Compute residuals
            computeResiduals();

            // Inserting data
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Roll the pointers P1 and P2
            waves_fw->roll();
            waves_bw->roll();
            waves_adj->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);

            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping
            waves_bw->forwardstepAcceleration(model, der);
            waves_bw->forwardstepStress(model, der);

            // Inserting data
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Roll the pointers P1 and P2
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal->writeCheck(waves_fw);
        }
        if (whatodo == restore)
        {
            optimal->readCheck(waves_fw);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));
     this->writeLog("\nGradient computation completed.");


	//Remove snapshot file
	optimal->removeCheck();

    // Free arrays
    free(adjsrc);

    result=LSRTM_OK;
    return result;
}

template<typename T>
LsrtmAcoustic2D<T>::~LsrtmAcoustic2D() {
    // Nothing here
}

// =============== ACOUSTIC 3D LSRTM CLASS =============== //

template<typename T>
LsrtmAcoustic3D<T>::LsrtmAcoustic3D(){
    sourceset = false;
    dataPset = false;
    modelset = false;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodPset = false;
    dataresPset = false;
    dataweightset = false;
}

template<typename T>
LsrtmAcoustic3D<T>::LsrtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataP, int order, int snapinc):Lsrtm<T>(order, snapinc){
    source = _source;
    dataP = _dataP;
    model = _model;
    modelset = true;
    sourceset = true;
    dataPset = true;
    datamodPset = false;
    dataresPset = false;
    vpgradset = false;
    rhogradset = false;
    wavgradset = false;
    dataweightset = false;
}

template<typename T>
void LsrtmAcoustic3D<T>::crossCorr(T *wsp, int pads, T* wrp, T* wrx, T* wry, T*wrz, int padr, T* Vp, T* Rho)
{
    if(!vpgradset && !rhogradset) rs_error("LsrtmAcoustic3D:crossCorr: No gradient set in lsrtm class");
	if(!vpgrad->getAllocated()) vpgrad->allocateImage();
    if(!rhograd->getAllocated()) rhograd->allocateImage();
	int ix, iy, iz;
    T *vpgraddata = vpgrad->getImagedata();
    T *rhograddata = rhograd->getImagedata();
    T mspx, mspy, mspz;
    T mrpx, mrpy, mrpz;
    T mrxx, mryy, mrzz;
    T L;
    T vpscale, rhoscale1, rhoscale2;
	int nx = vpgrad->getNx();
	int ny = vpgrad->getNy();
	int nz = vpgrad->getNz();
    int dx = vpgrad->getDx();
    int dy = vpgrad->getDy();
    int dz = vpgrad->getDz();
	int nxs = nx + 2*pads;
	int nxr = nx + 2*padr;
	int nys = ny + 2*pads;
	int nyr = ny + 2*padr;
    for (ix=1; ix<nx; ix++){
        for (iy=1; iy<ny; iy++){
            for (iz=1; iz<nz; iz++){
                L = Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                vpscale = -2.0/(Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]);
                rhoscale1 = -1.0/(Rho[km3D(ix, iy, iz)]*Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)]);
                rhoscale2 = -1.0/(Rho[km3D(ix, iy, iz)]*Rho[km3D(ix, iy, iz)]);
                mrxx = (wrx[kr3D(ix+padr, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mryy = (wry[kr3D(ix+padr, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrzz = (wrz[kr3D(ix+padr, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                mspx = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads-1, iy+pads, iz+pads)])/dx;
                mspy = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads, iy+pads-1, iz+pads)])/dy;
                mspz = (wsp[ks3D(ix+pads, iy+pads, iz+pads)] - wsp[ks3D(ix+pads, iy+pads, iz+pads-1)])/dz;
                mrpx = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mrpy = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrpz = (wrp[kr3D(ix+padr, iy+padr, iz+padr)] - wrp[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                vpgraddata[ki3D(ix,iy,iz)] -= vpscale*wsp[ks3D(ix+pads, iy+pads, iz+pads)]*L*(mrxx + mryy + mrzz);
                rhograddata[ki3D(ix,iy,iz)] -= rhoscale1*wsp[ks3D(ix+pads, iy+pads, iz+pads)]*L*(mrxx + mryy + mrzz);
                rhograddata[ki3D(ix,iy,iz)] -= rhoscale2*(mspx*mrpx + mspy*mrpy + mspz*mrpz);
            }	
        }
    }
}

template<typename T>
void LsrtmAcoustic3D<T>::computeMisfit(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }
    size_t itr, it;
    T res = 0.0;
    T misfit = 0.0;
    T *shaper;
    T *spiker;
    T *autocorr;
    T *crosscorr;
    T *modwei, *recwei;
    T norm; 
    T H;
    T stdev;

    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;

                for(it=0; it<nt; it++){
                    res = (-1.0)*(mod[I(it, itr)]*rec[I(it, itr)]/(norm1*norm2));
                   if(dataweightset)
                   {
                       res *= wei[I(it, itr)];
                   }
                   misfit += res;
                }
            }

            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);

                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                if (norm == 0.0) norm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it,stdev);
                    res = H*shaper[it];
                    misfit -= 0.5*res*res/norm;
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);

                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                if (norm == 0.0) norm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it,stdev);
                    res = H*shaper[it];
                    misfit += 0.5*res*res/norm;
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res = mod[I(it, itr)] - rec[I(it, itr)];
                    if(dataweightset)
                    {
                       res *= wei[I(it, itr)];
                   }
                   misfit += 0.5*res*res;
                }
            }
            break;
    }
    // Set the final misfit value
    this->setMisfit(misfit);
}

template<typename T>
void LsrtmAcoustic3D<T>::computeResiduals(){
    size_t ntr = datamodP->getNtrace();
    if(dataP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresP->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodP->getNt();
    if(dataP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresP->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* mod = datamodP->getData();
    T* rec = dataP->getData();
    T* res = dataresP->getData();
    T* wei = NULL;
    if(dataweightset) 
    {
        wei = dataweight->getData();
    }

    T *shaper;
    T *spiker;
    T *autocorr;
    T *crosscorr;
    T *modwei, *recwei;
    T norm; 
    T misfit;
    T H;
    T stdev;

    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T norm1, norm2, norm3;
            for(itr=0; itr<ntr; itr++){
                norm1 = 0.0;
                norm2 = 0.0;
                norm3 = 0.0;
                
                for(it=0; it<nt; it++){
                   norm1 += mod[I(it, itr)]*mod[I(it, itr)];
                   norm2 += rec[I(it, itr)]*rec[I(it, itr)];
                   norm3 += mod[I(it, itr)]*rec[I(it, itr)];
                }

                norm1 = sqrt(norm1);
                norm2 = sqrt(norm2);
                if(norm1 ==0 ) norm1= 1.0;
                if(norm2 ==0 ) norm2= 1.0;
                norm3 /= (norm1*norm2);

                for(it=0; it<nt; it++){
                    res[I(it, itr)]=(-1.0)*((rec[I(it, itr)]/(norm1*norm2)) - (mod[I(it, itr)]/(norm1*norm1))*norm3);
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }

            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                // If zero trace we need to add something to avoid division by zero
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);
                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                
                // Ensure norm different than 0.0
                if (norm == 0.0) norm = 1.0;
                misfit = 0.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    misfit += 0.5*(H*shaper[it])*(H*shaper[it])/norm;
                }

                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    res[I(it, itr)] = -1.0*(H*H - 2.0*misfit)*shaper[it]/norm;
                }
                this->stoep(nt, autocorr, &res[I(0, itr)], shaper, spiker);
                this->convolve(nt, 0, shaper, nt, 0, recwei, nt, 0, &res[I(0, itr)]);        
                if(dataweightset){
                    for (it=0; it < nt; it++)
                    {
                        res[I(it,itr)] *= wei[I(it,itr)];
                    }
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            modwei = (T *) calloc(nt, sizeof(T));
            recwei = (T *) calloc(nt, sizeof(T));
            shaper = (T *) calloc(nt, sizeof(T));
            spiker = (T *) calloc(nt, sizeof(T));
            autocorr = (T *) calloc(nt, sizeof(T));
            crosscorr = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightset){
                        modwei[it] = mod[I(it, itr)]*wei[I(it,itr)];
                        recwei[it] = rec[I(it, itr)]*wei[I(it,itr)];
                    }else{
                        modwei[it] = mod[I(it, itr)];
                        recwei[it] = rec[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recwei, nt, 0, recwei, nt, 0, autocorr);  /* for matrix */
                this->xcor(nt, 0, recwei, nt, 0, modwei, nt, 0, crosscorr); /* right hand side */
                // If zero trace we need to add something to avoid division by zero
                if (autocorr[0] == 0.0)  autocorr[0] = 1.0;
                autocorr[0] *= (1.0 + PNOISE_ACOUSTIC);			/* whiten */
                this->stoep(nt, autocorr, crosscorr, shaper, spiker);
                norm = 0.0; 
                for(it=0; it<nt; it++){
                    norm += shaper[it]*shaper[it];
                }
                
                // Ensure norm different than 0.0
                if (norm == 0.0) norm = 1.0;
                misfit = 0.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    misfit += 0.5*(H*shaper[it])*(H*shaper[it])/norm;
                }

                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    res[I(it, itr)] = (H*H - 2.0*misfit)*shaper[it]/norm;
                }
                this->stoep(nt, autocorr, &res[I(0, itr)], shaper, spiker);
                this->convolve(nt, 0, shaper, nt, 0, recwei, nt, 0, &res[I(0, itr)]);        
                if(dataweightset){
                    for (it=0; it < nt; it++)
                    {
                        res[I(it,itr)] *= wei[I(it,itr)];
                    }
                }
            }
            free(modwei);
            free(recwei);
            free(shaper);
            free(spiker);
            free(autocorr);
            free(crosscorr);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    res[I(it, itr)] = mod[I(it, itr)] - rec[I(it, itr)];
                   if(dataweightset)
                   {
                       res[I(it, itr)] *= wei[I(it, itr)]*wei[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int LsrtmAcoustic3D<T>::run(){
     int result = LSRTM_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Psnap;
     Psnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Psnap->openSnap(this->getSnapfile(), 'w'); // Create a new snapshot file
     Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field

     this->writeLog("Running 3D Acoustic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping
    	waves->forwardstepAcceleration(model, der);
    	waves->forwardstepStress(model, der);
    
    	// Inserting source 
    	waves->insertSource(model, source, SMAP, it);

        // Recording data 
        if(this->datamodPset){
            waves->recordData(this->datamodP, GMAP, it);
        }

    	//Writting out results to snapshot file
        Psnap->setData(waves->getP1(), 0); //Set Pressure as snap field
        Psnap->writeSnap(it);

    	// Roll the pointers P1 and P2
    	waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Psnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesAcoustic3D<T>>(model, nt, dt, ot);

    // Create image
    vpgrad->allocateImage();
    rhograd->allocateImage();

    // Compute misfit
    computeMisfit();

    if(!this->getNoreverse()){
        // Compute residuals
        computeResiduals();

        Psnap->openSnap(this->getSnapfile(), 'r');
        Psnap->allocSnap(0);

        this->writeLog("\nDoing reverse-time Loop.");
        // Loop over reverse time
        for(int it=0; it < nt; it++)
        {
            // Time stepping
            waves->forwardstepAcceleration(model, der);
            waves->forwardstepStress(model, der);

            // Inserting source 
            waves->insertSource(model, dataresP, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Psnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Psnap->getEnddiff()) % Psnap->getSnapinc()) == 0){
                T *wrp = waves->getP1();
                T* wrx = waves->getAx(); 
                T* wry = waves->getAy(); 
                T* wrz = waves->getAz(); 
                crossCorr(Psnap->getData(0), 0, wrp, wrx, wry, wrz, waves->getLpml(), model->getVp(), model->getR());
            }

            // Record wavelet gradient
            waves->recordData(this->wavgrad, SMAP, nt-1-it);

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);

            // Roll the pointers P1 and P2
            waves->roll();

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");
    } // End of reverse loop

	//Remove snapshot file
	Psnap->removeSnap();

    result=LSRTM_OK;
    return result;
}

template<typename T>
int LsrtmAcoustic3D<T>::run_optimal(){
     int result = LSRTM_ERR;
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesAcoustic3D<T>> waves_fw (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<WavesAcoustic3D<T>> waves_bw (new WavesAcoustic3D<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), waves_fw->getNy_pml(), waves_fw->getNz_pml(), waves_fw->getDx(), waves_fw->getDy(), waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
     vpgrad->allocateImage();


     this->writeLog("Running 3D Acoustic full-waveform inversion gradient with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal->getCapo();
        whatodo = optimal->revolve();
        capo = optimal->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping
                waves_fw->forwardstepAcceleration(model, der);
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertSource(model, source, SMAP, it);

                // Recording data 
                if(this->datamodPset && !reverse){
                    waves_fw->recordData(this->datamodP, GMAP, it);
                }

                // Roll the pointers P1 and P2
                waves_fw->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping
            waves_fw->forwardstepAcceleration(model, der);
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertSource(model, source, SMAP, capo);

            // Recording data 
            if(this->datamodPset){
                waves_fw->recordData(this->datamodP, GMAP, capo);
            }

            // Compute misfit
            computeMisfit();

            // Compute residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wry = waves_bw->getAy(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wry, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

            // Roll the pointers P1 and P2
            waves_fw->roll();
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);

            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping
            waves_bw->forwardstepAcceleration(model, der);
            waves_bw->forwardstepStress(model, der);

            // Inserting residuals
            waves_bw->insertSource(model, dataresP, GMAP, capo);

            /* Do Crosscorrelation */
            T *wsp = waves_fw->getP1();
            T *wrp = waves_bw->getP1();
            T* wrx = waves_bw->getAx(); 
            T* wry = waves_bw->getAy(); 
            T* wrz = waves_bw->getAz(); 
            crossCorr(wsp, waves_fw->getLpml(), wrp, wrx, wry, wrz, waves_bw->getLpml(), model->getVp(), model->getR());

            // Record wavelet gradient
            waves_bw->recordData(this->wavgrad, SMAP, capo);

            // Roll the pointers P1 and P2
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal->writeCheck(waves_fw);
        }
        if (whatodo == restore)
        {
            optimal->readCheck(waves_fw);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));
     this->writeLog("\nGradient computation completed.");


	//Remove snapshot file
	optimal->removeCheck();

    result=LSRTM_OK;
    return result;
}

template<typename T>
LsrtmAcoustic3D<T>::~LsrtmAcoustic3D() {
    // Nothing here
}

// =============== ELASTIC 2D LSRTM CLASS =============== //

template<typename T>
LsrtmElastic2D<T>::LsrtmElastic2D(){
    sourceset = false;
    dataUxset = false;
    dataUzset = false;
    modelset = false;
    vpgradset = false;
    vsgradset = false;
    datamodUxset = false;
    datamodUzset = false;
    dataresUxset = false;
    dataresUzset = false;
    dataweightxset = false;
    dataweightzset = false;
}

template<typename T>
LsrtmElastic2D<T>::LsrtmElastic2D(std::shared_ptr<ModelElastic2D<T>> _model, std::shared_ptr<Data2D<T>> _source, std::shared_ptr<Data2D<T>> _dataUx, std::shared_ptr<Data2D<T>> _dataUz, int order, int snapinc):Lsrtm<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUz = _dataUz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUzset = true;
    vpgradset = false;
    vsgradset = false;
    rhogradset = false;
    wavgradset = false;
    datamodUxset = false;
    datamodUzset = false;
    dataresUxset = false;
    dataresUzset = false;
    dataweightxset = false;
    dataweightzset = false;
}

template<typename T>
void LsrtmElastic2D<T>::crossCorr(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it)
{
    if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) rs_error("LsrtmElastic2D<T>::crossCorr: No gradient set for computation.");
	int ix, iz;
    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getUx1();
    T* wrz = waves_bw->getUz1();
    T* rsxx = waves_bw->getSxx();
    T* rszz = waves_bw->getSzz();
    T* rsxz = waves_bw->getSxz();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *wavgraddata = NULL;
	T *rhograddata = NULL;
	T msxx=0, mszz=0, msxz=0, mrxx=0, mrzz=0, mrxz=0, uderx=0, uderz=0;
	int nx;
	int nz;
	T dx;
	T dz;

    T* Rx = model->getRx();
    T* Rz = model->getRz();

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
	if(rhogradset){
		if(!rhograd->getAllocated()){
			rhograd->allocateImage();
		}
		rhograddata = rhograd->getImagedata();
	}

    int nt=0;
    int ntrace=0;
    int i=0;
    Point2D<int> *map=NULL;
    
	if(wavgradset){
		wavgraddata = wavgrad->getData();
        nt = wavgrad->getNt();
        ntrace = wavgrad->getNtrace();
        map = (wavgrad->getGeom())->getSmap();
	}

    // Getting sizes
    nx = waves_bw->getNx();
    nz = waves_bw->getNz();
    dx = waves_bw->getDx(); 
    dz = waves_bw->getDz(); 

	int nxs = nx+2*pads;
	int nxr = nx+2*padr;

    for (ix=1; ix<nx-1; ix++){
        for (iz=1; iz<nz-1; iz++){
            msxx = (wsx[ks2D(ix+pads, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads)])/dx;
            mszz = (wsz[ks2D(ix+pads, iz+pads)] - wsz[ks2D(ix+pads, iz+pads-1)])/dz;
            mrxx = (wrx[kr2D(ix+padr, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr)])/dx;
            mrzz = (wrz[kr2D(ix+padr, iz+padr)] - wrz[kr2D(ix+padr, iz+padr-1)])/dz;

            if(vpgradset || rhogradset){
                vpgraddata[ki2D(ix,iz)] -= (msxx + mszz) * (mrxx + mrzz);
            }

            if(vsgradset || rhogradset){
                msxz = 0.5*(wsx[ks2D(ix+pads, iz+pads+1)] - wsx[ks2D(ix+pads, iz+pads)])/dz;
                msxz += 0.5*(wsx[ks2D(ix+pads-1, iz+pads)] - wsx[ks2D(ix+pads-1, iz+pads-1)])/dz;
                msxz += 0.5*(wsz[ks2D(ix+pads+1, iz+pads)] - wsz[ks2D(ix+pads, iz+pads)])/dx;
                msxz += 0.5*(wsz[ks2D(ix+pads, iz+pads-1)] - wsz[ks2D(ix+pads-1, iz+pads-1)])/dx;

                mrxz = 0.5*(wrx[kr2D(ix+padr, iz+padr+1)] - wrx[kr2D(ix+padr, iz+padr)])/dz;
                mrxz += 0.5*(wrx[kr2D(ix+padr-1, iz+padr)] - wrx[kr2D(ix+padr-1, iz+padr-1)])/dz;
                mrxz += 0.5*(wrz[kr2D(ix+padr+1, iz+padr)] - wrz[kr2D(ix+padr, iz+padr)])/dx;
                mrxz += 0.5*(wrz[kr2D(ix+padr, iz+padr-1)] - wrz[kr2D(ix+padr-1, iz+padr-1)])/dx;

                vsgraddata[ki2D(ix,iz)] -= (2.0*msxx*mrxx + 2.0*mszz*mrzz + msxz*mrxz);
            }
            if(wavgradset){
                switch(wavgrad->getField()){
                    case PRESSURE:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = +1.0*(mrxx + mrzz);
                            }
                        }
                        break;
                    case VX:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = 0.5*(wrx[kr2D(ix+padr, iz+padr)] + wrx[kr2D(ix+padr-1, iz+padr)]);
                            }
                        }
                        break;
                    case VZ:
                        for (i=0; i < ntrace; i++) 
                        {
                            if((map[i].x == ix) && (map[i].y == iz))
                            {
                                wavgraddata[kwav(it,i)] = 0.5*(wrz[kr2D(ix+padr, iz+padr)] + wrz[kr2D(ix+padr, iz+padr-1)]);
                            }
                        }
                        break;
                    default:
                        break;
                }
            }

            if(rhogradset){
                uderx = 0.5*wsx[ks2D(ix+pads, iz+pads)]*Rx[kr2D(ix+padr, iz+padr)]*(rsxx[kr2D(ix+padr+1, iz+padr)] - rsxx[kr2D(ix+padr, iz+padr)])/dx;
                uderx += 0.5*wsx[ks2D(ix+pads-1, iz+pads)]*Rx[kr2D(ix+padr-1, iz+padr)]*(rsxx[kr2D(ix+padr, iz+padr)] - rsxx[kr2D(ix+padr-1, iz+padr)])/dx;
                uderx += 0.5*wsx[ks2D(ix+pads, iz+pads)]*Rx[kr2D(ix+padr, iz+padr)]*(rsxz[kr2D(ix+padr, iz+padr)] - rsxz[kr2D(ix+padr, iz+padr-1)])/dz;
                uderx += 0.5*wsx[ks2D(ix+pads-1, iz+pads)]*Rx[kr2D(ix+padr-1, iz+padr)]*(rsxz[kr2D(ix+padr-1, iz+padr)] - rsxz[kr2D(ix+padr-1, iz+padr-1)])/dz;

                uderz = 0.5*wsz[ks2D(ix+pads, iz+pads)]*Rz[kr2D(ix+padr, iz+padr)]*(rsxz[kr2D(ix+padr, iz+padr)] - rsxz[kr2D(ix+padr-1, iz+padr)])/dx;
                uderz += 0.5*wsz[ks2D(ix+pads, iz+pads-1)]*Rz[kr2D(ix+padr, iz+padr-1)]*(rsxz[kr2D(ix+padr, iz+padr-1)] - rsxz[kr2D(ix+padr-1, iz+padr-1)])/dx;
                uderz += 0.5*wsz[ks2D(ix+pads, iz+pads)]*Rz[kr2D(ix+padr, iz+padr)]*(rszz[kr2D(ix+padr, iz+padr+1)] - rszz[kr2D(ix+padr, iz+padr)])/dz;
                uderz += 0.5*wsz[ks2D(ix+pads, iz+pads-1)]*Rz[kr2D(ix+padr, iz+padr-1)]*(rszz[kr2D(ix+padr, iz+padr)] - rszz[kr2D(ix+padr, iz+padr-1)])/dz;
                rhograddata[ki2D(ix,iz)] -= (uderx + uderz);
            }
        }
    }	
}

template<typename T>
void LsrtmElastic2D<T>::scaleGrad(std::shared_ptr<ModelElastic2D<T>> model)
{
	int ix, iz;
    T* Vp = model->getVp();
    T* Vs = model->getVs();
    T* Rho = model->getR();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *rhograddata = NULL;
	T vpscale;
	T vsscale;
    T rhoscale1;
    T rhoscale2;
	int nx;
	int nz;

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
	if(rhogradset){
		if(!rhograd->getAllocated()){
			rhograd->allocateImage();
		}
		rhograddata = rhograd->getImagedata();
	}

        
    // Getting sizes
    nx = model->getNx();
    nz = model->getNz();

    T lambda = 0.0;
    T mu = 0.0;
    T rho = 0.0;

    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            vpscale = 2.0*Rho[km2D(ix, iz)]*Vp[km2D(ix, iz)];
            vsscale = 2.0*Rho[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            rhoscale1 = Vp[km2D(ix, iz)]*Vp[km2D(ix, iz)] -2.0*Vs[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            rhoscale2 = Vs[km2D(ix, iz)]*Vs[km2D(ix, iz)];
            lambda = vpgraddata[ki2D(ix,iz)];
            if(vsgradset || rhogradset){
                mu = vsgraddata[ki2D(ix,iz)];
            }

            if(vpgradset){
                vpgraddata[ki2D(ix,iz)] = vpscale*lambda;
            }

            if(vsgradset){
                vsgraddata[ki2D(ix,iz)] = vsscale*mu -2.0*vsscale*lambda; 
            }

            if(rhogradset){
                rho = rhograddata[ki2D(ix,iz)];
                rhograddata[ki2D(ix,iz)] = rhoscale1*lambda + rhoscale2*mu + rho;
            }

        }
    }	
}

template<typename T>
void LsrtmElastic2D<T>::computeMisfit(){
    size_t ntr = datamodUx->getNtrace();
    if(dataUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodUx->getNt();
    if(dataUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodUx->getData();
    T* recx = dataUx->getData();

    T* modz = datamodUz->getData();
    T* recz = dataUz->getData();
    T resx = 0.0;
    T resz = 0.0;
    T *weix = NULL;
    T *weiz = NULL;
    if(dataweightxset)
    {
        weix = dataweightx->getData();
    }

    if(dataweightzset)
    {
        weiz = dataweightz->getData();
    }

    T *shaperx, *shaperz;
    T *spikerx, *spikerz;
    T *autocorrx, *autocorrz;
    T *crosscorrx, *crosscorrz;
    T *modweix, *modweiz;
    T *recweix, *recweiz;
    T xnorm,znorm; 
    T H;
    T stdev;

    size_t itr, it;
    T misfit = 0.0;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }

                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resz*resz);
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2;
            T znorm1, znorm2;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                
                znorm1 = 0.0;
                znorm2 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;

                for(it=0; it<nt; it++){
                    resx=(-1.0)*(modx[I(it, itr)]*recx[I(it, itr)]/(xnorm1*xnorm2));
                    resz=((-1.0)*(modz[I(it, itr)]*recz[I(it, itr)]/(znorm1*znorm2)));
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += (resx + resz);
                }
            }
            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);

                xnorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it,stdev);
                    resx = H*shaperx[it];
                    resz = H*shaperz[it];
                    misfit -= 0.5*(resx*resx/xnorm + resz*resz/znorm);
                }
            }
            free(modweix);
            free(modweiz);
            free(recweix);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);

                xnorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it,stdev);
                    resx = H*shaperx[it];
                    resz = H*shaperz[it];
                    misfit += 0.5*(resx*resx/xnorm + resz*resz/znorm);
                }
            }
            free(modweix);
            free(modweiz);
            free(recweix);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;
        case ADAPTIVE_SINGLE:
            stdev = HMAX;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(2*nt, sizeof(T));
            recweix = (T *) calloc(2*nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweix[nt+it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweix[nt+it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweix[nt+it] = modz[I(it, itr)];
                        recweix[nt+it] = recz[I(it, itr)];
                    }
                }
                this->xcor(2*nt, 0, recweix, 2*nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(2*nt, 0, recweix, 2*nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);

                xnorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it,stdev);
                    resx = H*shaperx[it];
                    misfit += 0.5*(resx*resx/xnorm);
                }
            }
            free(modweix);
            free(recweix);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resz*resz);
                }
            }
            break;
    }

    // Set the final misfit value
    this->setMisfit(misfit);
}

template<typename T>
void LsrtmElastic2D<T>::computeResiduals(){
    size_t ntr = datamodUx->getNtrace();
    if(dataUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodUx->getNt();
    if(dataUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodUx->getData();
    T* recx = dataUx->getData();
    T* resx = dataresUx->getData();

    T* modz = datamodUz->getData();
    T* recz = dataUz->getData();
    T* resz = dataresUz->getData();
    T *weix = NULL;
    T *weiz = NULL;
    if(dataweightxset)
    {
        weix = dataweightx->getData();
    }
    if(dataweightzset)
    {
        weiz = dataweightz->getData();
    }

    T *shaperx, *shaperz;
    T *spikerx, *spikerz;
    T *autocorrx, *autocorrz;
    T *crosscorrx, *crosscorrz;
    T *modweix, *modweiz;
    T *recweix, *recweiz;
    T *fres;
    T pclip; 
    int pos;
    T xnorm,znorm; 
    T misfitx, misfitz;
    T H;
    T stdev;
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = modx[I(it, itr)] - recx[I(it, itr)];
                   resz[I(it, itr)] = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2, xnorm3;
            T znorm1, znorm2, znorm3;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                xnorm3 = 0.0;
                
                znorm1 = 0.0;
                znorm2 = 0.0;
                znorm3 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];
                    xnorm3 += modx[I(it, itr)]*recx[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                    znorm3 += modz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 == 0 ) xnorm1= 1.0;
                if(xnorm2 == 0 ) xnorm2= 1.0;
                xnorm3 /= (xnorm1*xnorm2);

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);

                if(znorm1 == 0 ) znorm1= 1.0;
                if(znorm2 == 0 ) znorm2= 1.0;
                znorm3 /= (znorm1*znorm2);

                for(it=0; it<nt; it++){
                    resx[I(it, itr)]=((-1.0)*((recx[I(it, itr)]/(xnorm1*xnorm2)) - (modx[I(it, itr)]/(xnorm1*xnorm1))*xnorm3));
                    resz[I(it, itr)]=((-1.0)*((recz[I(it, itr)]/(znorm1*znorm2)) - (modz[I(it, itr)]/(znorm1*znorm1))*znorm3));
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resx[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resz[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            if(fres[pos] > pclip) pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(resx[I(it, itr)]) >= pclip){
                        resx[I(it, itr)] = 0.0;
                        resz[I(it, itr)] = 0.0;
                    }
                }
            }

            free(fres);
            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);
                xnorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                misfitx = 0.0;
                misfitz = 0.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    misfitx += 0.5*(H*shaperx[it])*(H*shaperx[it])/xnorm;
                    misfitz += 0.5*(H*shaperz[it])*(H*shaperz[it])/znorm;
                }

                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    resx[I(it, itr)] = -1.0*(H*H - 2.0*misfitx)*shaperx[it]/xnorm;
                    resz[I(it, itr)] = -1.0*(H*H - 2.0*misfitz)*shaperz[it]/znorm;
                }
                this->stoep(nt, autocorrx, &resx[I(0, itr)], shaperx, spikerx);
                this->stoep(nt, autocorrz, &resz[I(0, itr)], shaperz, spikerz);
                this->convolve(nt, 0, shaperx, nt, 0, recweix, nt, 0, &resx[I(0, itr)]);        
                this->convolve(nt, 0, shaperz, nt, 0, recweiz, nt, 0, &resz[I(0, itr)]);        
                if(dataweightxset){
                    for(it=0; it<nt; it++)
                    {
                        resx[I(it, itr)] *= weix[I(it, itr)];
                    }
                }

                if(dataweightzset){
                    for(it=0; it<nt; it++)
                    {
                        resz[I(it, itr)] *= weiz[I(it, itr)];
                    }
                }
            }

            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resx[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resz[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            if(fres[pos] > pclip) pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(resx[I(it, itr)]) >= pclip){
                        resx[I(it, itr)] = 0.0;
                        resz[I(it, itr)] = 0.0;
                    }
                }
            }

            free(fres);
            free(modweix);
            free(modweiz);
            free(recweix);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;
        case ADAPTIVE_LINEAR:
            stdev = HMAX;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }

                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);
                xnorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                misfitx = 0.0;
                misfitz = 0.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    misfitx += 0.5*(H*shaperx[it])*(H*shaperx[it])/xnorm;
                    misfitz += 0.5*(H*shaperz[it])*(H*shaperz[it])/znorm;
                }

                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    resx[I(it, itr)] = (H*H - 2.0*misfitx)*shaperx[it]/xnorm;
                    resz[I(it, itr)] = (H*H - 2.0*misfitz)*shaperz[it]/znorm;
                }
                this->stoep(nt, autocorrx, &resx[I(0, itr)], shaperx, spikerx);
                this->stoep(nt, autocorrz, &resz[I(0, itr)], shaperz, spikerz);
                this->convolve(nt, 0, shaperx, nt, 0, recweix, nt, 0, &resx[I(0, itr)]);        
                this->convolve(nt, 0, shaperz, nt, 0, recweiz, nt, 0, &resz[I(0, itr)]);        
                if(dataweightxset){
                    for(it=0; it<nt; it++)
                    {
                        resx[I(it, itr)] *= weix[I(it, itr)];
                    }
                }
                if(dataweightzset){
                    for(it=0; it<nt; it++)
                    {
                        resz[I(it, itr)] *= weiz[I(it, itr)];
                    }
                }
            }
            // Thresholding
            fres = (T *) calloc(nt*ntr, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resx[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            pos = (int) (PCLIP*nt*ntr/100);
            pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    fres[I(it, itr)] = ABS(resz[I(it, itr)]);
                }
            }
            std::sort(fres, fres+nt*ntr); 
            if(fres[pos] > pclip) pclip = fres[pos];

            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                    if(ABS(resx[I(it, itr)]) >= pclip){
                        resx[I(it, itr)] = 0.0;
                        resz[I(it, itr)] = 0.0;
                    }
                }
            }

            free(fres);
            free(modweix);
            free(modweiz);
            free(recweix);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;
        case ADAPTIVE_SINGLE:
            stdev = HMAX;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(2*nt, sizeof(T));
            recweix = (T *) calloc(2*nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        recweix[it] = recx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweix[nt+it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweix[nt+it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweix[nt+it] = modz[I(it, itr)];
                        recweix[nt+it] = recz[I(it, itr)];
                    }
                }
                this->xcor(2*nt, 0, recweix, 2*nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(2*nt, 0, recweix, 2*nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                xnorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                misfitx = 0.0;
                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    misfitx += 0.5*(H*shaperx[it])*(H*shaperx[it])/xnorm;
                }

                for(it=0; it<nt; it++){
                    H = this->linear(it, stdev);
                    crosscorrx[it] = (H*H - 2.0*misfitx)*shaperx[it]/xnorm;
                }
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->convolve(nt, 0, shaperx, 2*nt, 0, recweix, 2*nt, 0, modweix);        
                if(dataweightxset){
                    for(it=0; it<nt; it++)
                    {
                        resx[I(it, itr)] = modweix[it]*weix[I(it, itr)];
                    }
                }else{
                    for(it=0; it<nt; it++)
                    {
                        resx[I(it, itr)] = modweix[it];
                    }
                }
                if(dataweightzset){
                    for(it=0; it<nt; it++)
                    {
                        resz[I(it, itr)] = modweix[nt+it]*weiz[I(it, itr)];
                    }
                }else{
                    for(it=0; it<nt; it++)
                    {
                        resz[I(it, itr)] = modweix[nt+it];
                    }
                }
            }
            free(modweix);
            free(recweix);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            break;

        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = modx[I(it, itr)] - recx[I(it, itr)];
                   resz[I(it, itr)] = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            break;
    }
}

template<typename T>
int LsrtmElastic2D<T>::run(){
	int result = LSRTM_ERR;
	if(!vpgradset && !vsgradset && !wavgradset) {
		rs_warning("LsrtmElastic2D::run: No image set");
		return result;
	}
	int nt;
	float dt;
	float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), 1, waves->getNz_pml(), waves->getDx(), 1.0, waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot2D<T>> Uxsnap;
     Uxsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Uxsnap->openSnap(this->getSnapfile() + "-ux", 'w'); // Create a new snapshot file
     Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field

     std::shared_ptr<Snapshot2D<T>> Uzsnap;
     Uzsnap = std::make_shared<Snapshot2D<T>>(waves, this->getSnapinc());
     Uzsnap->openSnap(this->getSnapfile() + "-uz", 'w'); // Create a new snapshot file
     Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field

     this->writeLog("Running 2D Elastic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field
        Uxsnap->writeSnap(it);

        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field
        Uzsnap->writeSnap(it);

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting force source 
    	waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting force source 
    	waves->insertForcesource(model, source, SMAP, it);

        // Recording data (Ux)
        if(this->datamodUxset){
            waves->recordData(model, this->datamodUx, GMAP, it);
        }

        // Recording data (Uz)
        if(this->datamodUzset){
            waves->recordData(model, this->datamodUz, GMAP, it);
        }

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Uxsnap->closeSnap();
    Uzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic2D_DS<T>>(model, nt, dt, ot);

    // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

    Uxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
    Uxsnap->allocSnap(0);

    Uzsnap->openSnap(this->getSnapfile() + "-uz", 'r');
    Uzsnap->allocSnap(0);

    // Compute misfit
    computeMisfit();

    // Compute Residuals
    computeResiduals();

     this->writeLog("\nDoing reverse-time Loop.");
    // Loop over reverse time
    for(int it=0; it < nt; it++)
    {
    	// Time stepping 
    	waves->forwardstepStress(model, der);

    	// Inserting residuals
    	waves->insertPressuresource(model, dataresUx, GMAP, (nt - 1 - it));
    	waves->insertPressuresource(model, dataresUz, GMAP, (nt - 1 - it));

    	// Time stepping 
    	waves->forwardstepDisplacement(model, der);

    	// Inserting residuals
    	waves->insertForcesource(model, dataresUx, GMAP, (nt - 1 - it));
    	waves->insertForcesource(model, dataresUz, GMAP, (nt - 1 - it));

        //Read forward snapshot
        Uxsnap->readSnap(nt - 1 - it);
        Uzsnap->readSnap(nt - 1 - it);

        // Do Crosscorrelation
        crossCorr(Uxsnap->getData(0), Uzsnap->getData(0), 0, waves, model, (nt - 1 - it));

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }
     this->writeLog("\nGradient computation completed.");

    // Scale gradients
    this->scaleGrad(model);
    
	//Remove snapshot file
	Uxsnap->removeSnap();
	Uzsnap->removeSnap();

    result=LSRTM_OK;
    return result;
}

template<typename T>
int LsrtmElastic2D<T>::run_optimal(){
     int result = LSRTM_ERR;
     if(!vpgradset && !vsgradset && !wavgradset) {
         rs_warning("LsrtmElastic2D::run: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic2D_DS<T>> waves_fw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic2D_DS<T>> waves_bw (new WavesElastic2D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), 1, waves_fw->getNz_pml(), waves_fw->getDx(), 1.0, waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

     this->writeLog("Running 2D Elastic full-waveform inversion gradient with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal->getCapo();
        whatodo = optimal->revolve();
        capo = optimal->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping stress
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertPressuresource(model, source, SMAP, it);

                // Time stepping displacement
                waves_fw->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw->insertForcesource(model, source, SMAP, it);

                // Recording data (Ux)
                if(this->datamodUxset && !reverse){
                    waves_fw->recordData(model, this->datamodUx, GMAP, it);
                }

                // Recording data (Uz)
                if(this->datamodUzset && !reverse){
                    waves_fw->recordData(model, this->datamodUz, GMAP, it);
                }

                // Roll pointers
                waves_fw->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }

        }
        if (whatodo == firsturn)
        {
            // Time stepping stess
            waves_fw->forwardstepStress(model, der);

            // Inserting pressure source 
            waves_fw->insertPressuresource(model, source, SMAP, capo);
            
            // Time stepping displacement
            waves_fw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw->insertForcesource(model, source, SMAP, capo);

            // Recording data (Ux)
            if(this->datamodUxset){
                waves_fw->recordData(model, this->datamodUx, GMAP, capo);
            }

            // Recording data (Uz)
            if(this->datamodUzset){
                waves_fw->recordData(model, this->datamodUz, GMAP, capo);
            }

            // Compute misfit
            computeMisfit();

            // Compute Residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertPressuresource(model, dataresUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUz, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getUx1();
            T *wsz = waves_fw->getUz1();
            crossCorr(wsx, wsz, waves_fw->getLpml(), waves_bw, model, capo);

            // Roll pointers
            waves_fw->roll();
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);
      
            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping stress
            waves_bw->forwardstepStress(model, der);

            // Inserting residuals
            waves_bw->insertPressuresource(model, dataresUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUz, GMAP, capo);

            // Time stepping displacement
            waves_bw->forwardstepDisplacement(model, der);

            // Inserting residuals
            waves_bw->insertForcesource(model, dataresUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getUx1();
            T *wsz = waves_fw->getUz1();
            crossCorr(wsx, wsz, waves_fw->getLpml(), waves_bw, model, capo);

            // Roll pointers
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal->writeCheck(waves_fw);
        }
        if (whatodo == restore)
        {
            optimal->readCheck(waves_fw);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));
     this->writeLog("\nGradient computation completed.");
    
    // Scale gradients
    this->scaleGrad(model);

	//Remove snapshot file
	optimal->removeCheck();

    result=LSRTM_OK;
    return result;
}

template<typename T>
LsrtmElastic2D<T>::~LsrtmElastic2D() {
    // Nothing here
}

// =============== ELASTIC 3D LSRTM CLASS =============== //

template<typename T>
LsrtmElastic3D<T>::LsrtmElastic3D(){
    sourceset = false;
    dataUxset = false;
    dataUyset = false;
    dataUzset = false;
    modelset = false;
    rhogradset = false;
    vpgradset = false;
    vsgradset = false;
    wavgradset = false;
    datamodUxset = false;
    datamodUyset = false;
    datamodUzset = false;
    dataresUxset = false;
    dataresUyset = false;
    dataresUzset = false;
    dataweightxset = false;
    dataweightyset = false;
    dataweightzset = false;
}

template<typename T>
LsrtmElastic3D<T>::LsrtmElastic3D(std::shared_ptr<ModelElastic3D<T>> _model, std::shared_ptr<Data3D<T>> _source, std::shared_ptr<Data3D<T>> _dataUx, std::shared_ptr<Data3D<T>> _dataUy, std::shared_ptr<Data3D<T>> _dataUz, int order, int snapinc):Lsrtm<T>(order, snapinc){
    source = _source;
    dataUx = _dataUx;
    dataUy = _dataUy;
    dataUz = _dataUz;
    model = _model;
    sourceset = true;
    modelset = true;
    dataUxset = true;
    dataUyset = true;
    dataUzset = true;
    rhogradset = false;
    vpgradset = false;
    vsgradset = false;
    wavgradset = false;
    datamodUxset = false;
    datamodUyset = false;
    datamodUzset = false;
    dataresUxset = false;
    dataresUyset = false;
    dataresUzset = false;
    dataweightxset = false;
    dataweightyset = false;
    dataweightzset = false;
}

template<typename T>
void LsrtmElastic3D<T>::crossCorr(T *wsx, T*wsy, T *wsz, int pads, std::shared_ptr<WavesElastic3D_DS<T>> waves_bw, std::shared_ptr<ModelElastic3D<T>> model, int it)
{
    if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) rs_error("LsrtmElastic3D<T>::crossCorr: No gradient set for computation.");
    int ix, iy, iz;

    int padr = waves_bw->getLpml();
    T* wrx = waves_bw->getUx1();
    T* wry = waves_bw->getUy1();
    T* wrz = waves_bw->getUz1();
    T* rsxx = waves_bw->getSxx();
    T* rsyy = waves_bw->getSyy();
    T* rszz = waves_bw->getSzz();
    T* rsyz = waves_bw->getSyz();
    T* rsxz = waves_bw->getSxz();
    T* rsxy = waves_bw->getSxy();

    T* Rx = model->getRx();
    T* Ry = model->getRy();
    T* Rz = model->getRz();

    T *vpgraddata = NULL; 
    T *vsgraddata = NULL;
    T *wavgraddata = NULL;
    T *rhograddata = NULL;

    T msxx=0, msyy=0, mszz=0, msyz=0, msxz=0, msxy=0, mrxx=0, mryy=0, mrzz=0, mryz=0, mrxz=0, mrxy=0, uderx=0, udery=0, uderz=0;
    int nx;
    int ny;
    int nz;
    T dx;
    T dy;
    T dz;

    if(vpgradset){
        if(!vpgrad->getAllocated()){
            vpgrad->allocateImage();
        }
        vpgraddata = vpgrad->getImagedata();
    }
    if(vsgradset){
        if(!vsgrad->getAllocated()){
            vsgrad->allocateImage();
        }
        vsgraddata = vsgrad->getImagedata();
    }
    if(rhogradset){
        if(!rhograd->getAllocated()){
            rhograd->allocateImage();
        }
        rhograddata = rhograd->getImagedata();
    }


    int nt=0;
    int ntrace=0;
    int i=0;
    Point3D<int> *map=NULL;

    if(wavgradset){
        wavgraddata = wavgrad->getData();
        nt = wavgrad->getNt();
        ntrace = wavgrad->getNtrace();
        map = (wavgrad->getGeom())->getSmap();
    }

    // Getting sizes
    nx = waves_bw->getNx();
    ny = waves_bw->getNy();
    nz = waves_bw->getNz();
    dx = waves_bw->getDx(); 
    dy = waves_bw->getDy(); 
    dz = waves_bw->getDz(); 

    int nxs = nx + 2*pads;
    int nxr = nx + 2*padr;
    int nys = ny + 2*pads;
    int nyr = ny + 2*padr;

    for (ix=1; ix<nx-1; ix++){
        for (iy=1; iy<ny-1; iy++){
            for (iz=1; iz<nz-1; iz++){
                msxx = (wsx[ks3D(ix+pads, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads)])/dx;
                msyy = (wsy[ks3D(ix+pads, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads)])/dy;
                mszz = (wsz[ks3D(ix+pads, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads-1)])/dz;
                mrxx = (wrx[kr3D(ix+padr, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                mryy = (wry[kr3D(ix+padr, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                mrzz = (wrz[kr3D(ix+padr, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;

                if(vpgradset || rhogradset){
                    vpgraddata[ki3D(ix,iy,iz)] -= (msxx + msyy + mszz) * (mrxx + mryy + mrzz);
                }

                if(vsgradset || rhogradset){
                    msyz = 0.5*(wsz[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msyz += 0.5*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dy;
                    msyz += 0.5*(wsy[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msyz += 0.5*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads, iy+pads-1, iz+pads-1)])/dz;

                    mryz = 0.5*(wrz[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mryz += 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dy;
                    mryz += 0.5*(wry[kr3D(ix+padr, iy+padr, iz+padr+1)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mryz += 0.5*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dz;

                    msxz = 0.5*(wsx[ks3D(ix+pads, iy+pads, iz+pads+1)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dz;
                    msxz += 0.5*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dz;
                    msxz += 0.5*(wsz[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsz[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxz += 0.5*(wsz[ks3D(ix+pads, iy+pads, iz+pads-1)] - wsz[ks3D(ix+pads-1, iy+pads, iz+pads-1)])/dx;

                    mrxz = 0.5*(wrx[kr3D(ix+padr, iy+padr, iz+padr+1)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    mrxz += 0.5*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dz;
                    mrxz += 0.5*(wrz[kr3D(ix+padr+1, iy+padr, iz+padr)] - wrz[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxz += 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr-1)] - wrz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dx;

                    msxy = 0.5*(wsx[ks3D(ix+pads, iy+pads+1, iz+pads)] - wsx[ks3D(ix+pads, iy+pads, iz+pads)])/dy;
                    msxy += 0.5*(wsx[ks3D(ix+pads-1, iy+pads, iz+pads)] - wsx[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dy;
                    msxy += 0.5*(wsy[ks3D(ix+pads+1, iy+pads, iz+pads)] - wsy[ks3D(ix+pads, iy+pads, iz+pads)])/dx;
                    msxy += 0.5*(wsy[ks3D(ix+pads, iy+pads-1, iz+pads)] - wsy[ks3D(ix+pads-1, iy+pads-1, iz+pads)])/dx;

                    mrxy = 0.5*(wrx[kr3D(ix+padr, iy+padr+1, iz+padr)] - wrx[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    mrxy += 0.5*(wrx[kr3D(ix+padr-1, iy+padr, iz+padr)] - wrx[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dy;
                    mrxy += 0.5*(wry[kr3D(ix+padr+1, iy+padr, iz+padr)] - wry[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    mrxy += 0.5*(wry[kr3D(ix+padr, iy+padr-1, iz+padr)] - wry[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dx;
                }
                if(vsgradset){
                    vsgraddata[ki3D(ix,iy,iz)] -= (2.0*msxx*mrxx + 2.0*msyy*mryy + 2.0*mszz*mrzz + msyz*mryz + msxz*mrxz + msxy*mrxy);
                }
                if(wavgradset){
                    switch(wavgrad->getField()){
                        case PRESSURE:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 1.0*(mrxx + mryy + mrzz);
                                }
                            }
                            break;
                        case VX:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wrx[kr3D(ix+padr, iy+padr, iz+padr)] + wrx[kr3D(ix+padr-1, iy+padr, iz+padr)]);
                                }
                            }
                            break;
                        case VY:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wry[kr3D(ix+padr, iy+padr, iz+padr)] + wry[kr3D(ix+padr, iy+padr-1, iz+padr)]);
                                }
                            }
                            break;
                        case VZ:
                            for (i=0; i < ntrace; i++) 
                            {
                                if((map[i].x == ix) && (map[i].y == iy) && (map[i].z == iz))
                                {
                                    wavgraddata[kwav(it,i)] = 0.5*(wrz[kr3D(ix+padr, iy+padr, iz+padr)] + wrz[kr3D(ix+padr, iy+padr, iz+padr-1)]);
                                }
                            }
                            break;
                        default:
                            break;
                    }
                }

                if(rhogradset){
                    uderx = (0.5)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxx[kr3D(ix+padr+1, iy+padr, iz+padr)] - rsxx[kr3D(ix+padr, iy+padr, iz+padr)])/dx;
                    uderx += (0.5)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxx[kr3D(ix+padr, iy+padr, iz+padr)] - rsxx[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                    uderx += (0.5)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                    uderx += (0.5)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr-1, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dy;
                    uderx += (0.5)*wsx[ks3D(ix+pads, iy+pads, iz+pads)]*Rx[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                    uderx += (0.5)*wsx[ks3D(ix+pads-1, iy+pads, iz+pads)]*Rx[kr3D(ix+padr-1, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr-1, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dz;

                    udery = (0.5)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                    udery += (0.5)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsxy[kr3D(ix+padr, iy+padr-1, iz+padr)] - rsxy[kr3D(ix+padr-1, iy+padr-1, iz+padr)])/dx;
                    udery += (0.5)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyy[kr3D(ix+padr, iy+padr+1, iz+padr)] - rsyy[kr3D(ix+padr, iy+padr, iz+padr)])/dy;
                    udery += (0.5)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsyy[kr3D(ix+padr, iy+padr, iz+padr)] - rsyy[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                    udery += (0.5)*wsy[ks3D(ix+pads, iy+pads, iz+pads)]*Ry[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;
                    udery += (0.5)*wsy[ks3D(ix+pads, iy+pads-1, iz+pads)]*Ry[kr3D(ix+padr, iy+padr-1, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr-1, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dz;

                    uderz = (0.5)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr)])/dx;
                    uderz += (0.5)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rsxz[kr3D(ix+padr, iy+padr, iz+padr-1)] - rsxz[kr3D(ix+padr-1, iy+padr, iz+padr-1)])/dx;
                    uderz += (0.5)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr)])/dy;
                    uderz += (0.5)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rsyz[kr3D(ix+padr, iy+padr, iz+padr-1)] - rsyz[kr3D(ix+padr, iy+padr-1, iz+padr-1)])/dy;
                    uderz += (0.5)*wsz[ks3D(ix+pads, iy+pads, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr)]*(rszz[kr3D(ix+padr, iy+padr, iz+padr+1)] - rszz[kr3D(ix+padr, iy+padr, iz+padr)])/dz;
                    uderz += (0.5)*wsz[ks3D(ix+pads, iy+pads-1, iz+pads)]*Rz[kr3D(ix+padr, iy+padr, iz+padr-1)]*(rszz[kr3D(ix+padr, iy+padr, iz+padr)] - rszz[kr3D(ix+padr, iy+padr, iz+padr-1)])/dz;

                    rhograddata[ki3D(ix,iy,iz)] -= (uderx + udery + uderz);
                }

            }
        }	
    }
}

template<typename T>
void LsrtmElastic3D<T>::scaleGrad(std::shared_ptr<ModelElastic3D<T>> model)
{
	int ix, iy, iz;
    T* Vp = model->getVp();
    T* Vs = model->getVs();
    T* Rho = model->getR();

	T *vpgraddata = NULL; 
	T *vsgraddata = NULL;
	T *rhograddata = NULL;
	T vpscale;
	T vsscale;
    T rhoscale1;
    T rhoscale2;
	int nx;
	int ny;
	int nz;

	if(vpgradset){
		if(!vpgrad->getAllocated()){
			vpgrad->allocateImage();
		}
		vpgraddata = vpgrad->getImagedata();
	}
	if(vsgradset){
		if(!vsgrad->getAllocated()){
			vsgrad->allocateImage();
		}
		vsgraddata = vsgrad->getImagedata();
	}
	if(rhogradset){
		if(!rhograd->getAllocated()){
			rhograd->allocateImage();
		}
		rhograddata = rhograd->getImagedata();
	}

        
    // Getting sizes
    nx = model->getNx();
    ny = model->getNy();
    nz = model->getNz();

    T lambda = 0.0;
    T mu = 0.0;
    T rho = 0.0;

    for (ix=0; ix<nx; ix++){
        for (iy=0; iy<ny; iy++){
            for (iz=0; iz<nz; iz++){
                vpscale = 2.0*Rho[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                vsscale = 2.0*Rho[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];
                rhoscale1 = Vp[km3D(ix, iy, iz)]*Vp[km3D(ix, iy, iz)];
                rhoscale2 = Vs[km3D(ix, iy, iz)]*Vs[km3D(ix, iy, iz)];

                lambda = vpgraddata[ki3D(ix,iy,iz)];
                if(vsgradset || rhogradset){
                    mu = vsgraddata[ki3D(ix,iy,iz)];
                }

                if(vpgradset){
                    vpgraddata[ki3D(ix,iy,iz)] = vpscale*lambda;
                }

                if(vsgradset){
                    vsgraddata[ki3D(ix,iy,iz)] = vsscale*mu -2.0*vsscale*lambda; 
                }

                if(rhogradset){
                    rho = rhograddata[ki3D(ix,iy,iz)];
                    rhograddata[ki3D(ix,iy,iz)] = rhoscale1*lambda + rhoscale2*mu + rho;
                }

            }
        }	
    }
}

template<typename T>
void LsrtmElastic3D<T>::computeMisfit(){
    size_t ntr = datamodUx->getNtrace();
    if(dataUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodUx->getNt();
    if(dataUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodUx->getData();
    T* recx = dataUx->getData();

    T* mody = datamodUy->getData();
    T* recy = dataUy->getData();

    T* modz = datamodUz->getData();
    T* recz = dataUz->getData();
    T resx = 0.0;
    T resy = 0.0;
    T resz = 0.0;
    T *weix = NULL;
    T *weiy = NULL;
    T *weiz = NULL;
    if(dataweightxset)
    {
        weix = dataweightx->getData();
    }
    if(dataweightyset)
    {
        weiy = dataweighty->getData();
    }
    if(dataweightzset)
    {
        weiz = dataweightz->getData();
    }

    T *shaperx,*shapery,*shaperz;
    T *spikerx,*spikery,*spikerz;
    T *autocorrx,*autocorry,*autocorrz;
    T *crosscorrx,*crosscorry,*crosscorrz;
    T *modweix, *modweiy, *modweiz;
    T *recweix, *recweiy, *recweiz;
    T xnorm, ynorm, znorm; 
    T H;
    T stdev;

    size_t itr, it;
    T misfit = 0.0;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resy = mody[I(it, itr)] - recy[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy *= weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resy*resy + resz*resz);
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2;
            T ynorm1, ynorm2;
            T znorm1, znorm2;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;

                ynorm1 = 0.0;
                ynorm2 = 0.0;
                
                znorm1 = 0.0;
                znorm2 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];

                    ynorm1 += mody[I(it, itr)]*mody[I(it, itr)];
                    ynorm2 += recy[I(it, itr)]*recy[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;

                ynorm1 = sqrt(ynorm1);
                ynorm2 = sqrt(ynorm2);
                if(ynorm1 ==0 ) ynorm1= 1.0;
                if(ynorm2 ==0 ) ynorm2= 1.0;

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;

                for(it=0; it<nt; it++){
                    resx=(-1.0)*(modx[I(it, itr)]*recx[I(it, itr)]/(xnorm1*xnorm2));
                    resy=(-1.0)*(mody[I(it, itr)]*recy[I(it, itr)]/(ynorm1*ynorm2));
                    resz=((-1.0)*(modz[I(it, itr)]*recz[I(it, itr)]/(znorm1*znorm2)));
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy *= weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += (resx + resy + resz);
                }
            }
            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shapery = (T *) calloc(nt, sizeof(T));
            spikery = (T *) calloc(nt, sizeof(T));
            autocorry = (T *) calloc(nt, sizeof(T));
            crosscorry = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiy = (T *) calloc(nt, sizeof(T));
            recweiy = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightyset){
                        recweiy[it] = recy[I(it, itr)]*weiy[I(it,itr)];
                        recweiy[it] = recy[I(it, itr)]*weiy[I(it,itr)];
                    }else{
                        modweiy[it] = mody[I(it, itr)];
                        recweiy[it] = recy[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiy, nt, 0, recweiy, nt, 0, autocorry);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiy, nt, 0, modweiy, nt, 0, crosscorry); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorry[0] == 0.0)  autocorry[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorry[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorry, crosscorry, shapery, spikery);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);

                xnorm = 0.0; 
                ynorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    ynorm += shapery[it]*shapery[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (ynorm == 0.0) ynorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it,stdev);
                    resx = H*shaperx[it];
                    resy = H*shapery[it];
                    resz = H*shaperz[it];
                    misfit += 0.5*(resx*resx/xnorm + resy*resy/ynorm + resz*resz/znorm);
                }
            }
            free(modweix);
            free(modweiy);
            free(modweiz);
            free(recweix);
            free(recweiy);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shapery);
            free(spikery);
            free(autocorry);
            free(crosscorry);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;
        default:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx = modx[I(it, itr)] - recx[I(it, itr)];
                   resy = mody[I(it, itr)] - recy[I(it, itr)];
                   resz = modz[I(it, itr)] - recz[I(it, itr)];
                   if(dataweightxset)
                   {
                       resx *= weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy *= weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz *= weiz[I(it, itr)];
                   }
                   misfit += 0.5*(resx*resx + resy*resy + resz*resz);
                }
            }
            break;
    }

    // Set the final misfit value
    this->setMisfit(misfit);
}


template<typename T>
void LsrtmElastic3D<T>::computeResiduals(){
    size_t ntr = datamodUx->getNtrace();
    if(dataUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUx->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    size_t nt = datamodUx->getNt();
    if(dataUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUx->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUy->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUy->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    if(dataUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and recorded data.");
    if(dataresUz->getNtrace() != ntr) rs_error("Mismatch between number of traces in the modelled and residual data.");
    if(dataUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and recorded data.");
    if(dataresUz->getNt() != nt) rs_error("Mismatch between number of time samples in the modelled and residual data.");

    T* modx = datamodUx->getData();
    T* recx = dataUx->getData();
    T* resx = dataresUx->getData();

    T* mody = datamodUy->getData();
    T* recy = dataUy->getData();
    T* resy = dataresUy->getData();

    T* modz = datamodUz->getData();
    T* recz = dataUz->getData();
    T* resz = dataresUz->getData();
    T *weix = NULL;
    T *weiy = NULL;
    T *weiz = NULL;
    if(dataweightxset)
    {
        weix = dataweightx->getData();
    }
    if(dataweightyset)
    {
        weiy = dataweighty->getData();
    }
    if(dataweightzset)
    {
        weiz = dataweightz->getData();
    }

    T *shaperx,*shapery,*shaperz;
    T *spikerx,*spikery,*spikerz;
    T *autocorrx,*autocorry,*autocorrz;
    T *crosscorrx,*crosscorry,*crosscorrz;
    T *modweix, *modweiy, *modweiz;
    T *recweix, *recweiy, *recweiz;
    T xnorm, ynorm, znorm; 
    T misfitx, misfity, misfitz;
    T H;
    T stdev;
    size_t itr, it;
    Index I(nt, ntr);
    switch(this->getMisfit_type()){
        case DIFFERENCE:
            for(itr=0; itr<ntr; itr++){
                for(it=0; it<nt; it++){
                   resx[I(it, itr)] = (modx[I(it, itr)] - recx[I(it, itr)]);
                   resy[I(it, itr)] = (mody[I(it, itr)] - recy[I(it, itr)]);
                   resz[I(it, itr)] = (modz[I(it, itr)] - recz[I(it, itr)]);
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy[I(it, itr)] *= weiy[I(it, itr)]*weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            break;
        case CORRELATION:
            T xnorm1, xnorm2, xnorm3;
            T ynorm1, ynorm2, ynorm3;
            T znorm1, znorm2, znorm3;
            for(itr=0; itr<ntr; itr++){
                xnorm1 = 0.0;
                xnorm2 = 0.0;
                xnorm3 = 0.0;

                ynorm1 = 0.0;
                ynorm2 = 0.0;
                ynorm3 = 0.0;

                znorm1 = 0.0;
                znorm2 = 0.0;
                znorm3 = 0.0;
                for(it=0; it<nt; it++){
                    xnorm1 += modx[I(it, itr)]*modx[I(it, itr)];
                    xnorm2 += recx[I(it, itr)]*recx[I(it, itr)];
                    xnorm3 += modx[I(it, itr)]*recx[I(it, itr)];

                    ynorm1 += mody[I(it, itr)]*mody[I(it, itr)];
                    ynorm2 += recy[I(it, itr)]*recy[I(it, itr)];
                    ynorm3 += mody[I(it, itr)]*recy[I(it, itr)];

                    znorm1 += modz[I(it, itr)]*modz[I(it, itr)];
                    znorm2 += recz[I(it, itr)]*recz[I(it, itr)];
                    znorm3 += modz[I(it, itr)]*recz[I(it, itr)];
                }

                xnorm1 = sqrt(xnorm1);
                xnorm2 = sqrt(xnorm2);
                if(xnorm1 ==0 ) xnorm1= 1.0;
                if(xnorm2 ==0 ) xnorm2= 1.0;
                xnorm3 /= (xnorm1*xnorm2);

                ynorm1 = sqrt(ynorm1);
                ynorm2 = sqrt(ynorm2);
                if(ynorm1 ==0 ) ynorm1= 1.0;
                if(ynorm2 ==0 ) ynorm2= 1.0;
                ynorm3 /= (ynorm1*ynorm2);

                znorm1 = sqrt(znorm1);
                znorm2 = sqrt(znorm2);
                if(znorm1 ==0 ) znorm1= 1.0;
                if(znorm2 ==0 ) znorm2= 1.0;
                znorm3 /= (znorm1*znorm2);

                for(it=1; it<nt; it++){
                    resx[I(it, itr)]=((-1.0)*((recx[I(it, itr)]/(xnorm1*xnorm2)) - (modx[I(it, itr)]/(xnorm1*xnorm1))*xnorm3));
                    resy[I(it, itr)]=((-1.0)*((recy[I(it, itr)]/(ynorm1*ynorm2)) - (mody[I(it, itr)]/(ynorm1*ynorm1))*ynorm3));
                    resz[I(it, itr)]=((-1.0)*((recz[I(it, itr)]/(znorm1*znorm2)) - (modz[I(it, itr)]/(znorm1*znorm1))*znorm3));
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy[I(it, itr)] *= weiy[I(it, itr)]*weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            break;
        case ADAPTIVE_GAUSS:
            stdev = STDEV*nt;
            shaperx = (T *) calloc(nt, sizeof(T));
            spikerx = (T *) calloc(nt, sizeof(T));
            autocorrx = (T *) calloc(nt, sizeof(T));
            crosscorrx = (T *) calloc(nt, sizeof(T));
            shapery = (T *) calloc(nt, sizeof(T));
            spikery = (T *) calloc(nt, sizeof(T));
            autocorry = (T *) calloc(nt, sizeof(T));
            crosscorry = (T *) calloc(nt, sizeof(T));
            shaperz = (T *) calloc(nt, sizeof(T));
            spikerz = (T *) calloc(nt, sizeof(T));
            autocorrz = (T *) calloc(nt, sizeof(T));
            crosscorrz = (T *) calloc(nt, sizeof(T));
            modweix = (T *) calloc(nt, sizeof(T));
            recweix = (T *) calloc(nt, sizeof(T));
            modweiy = (T *) calloc(nt, sizeof(T));
            recweiy = (T *) calloc(nt, sizeof(T));
            modweiz = (T *) calloc(nt, sizeof(T));
            recweiz = (T *) calloc(nt, sizeof(T));
            for(itr=0; itr<ntr; itr++){
                for (it=0; it < nt; it++)
                {
                    if(dataweightxset){
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                        modweix[it] = modx[I(it, itr)]*weix[I(it,itr)];
                    }else{
                        modweix[it] = modx[I(it, itr)];
                        recweix[it] = recx[I(it, itr)];
                    }
                    if(dataweightyset){
                        recweiy[it] = recy[I(it, itr)]*weiy[I(it,itr)];
                        recweiy[it] = recy[I(it, itr)]*weiy[I(it,itr)];
                    }else{
                        modweiy[it] = mody[I(it, itr)];
                        recweiy[it] = recy[I(it, itr)];
                    }
                    if(dataweightzset){
                        modweiz[it] = modz[I(it, itr)]*weiz[I(it,itr)];
                        recweiz[it] = recz[I(it, itr)]*weiz[I(it,itr)];
                    }else{
                        modweiz[it] = modz[I(it, itr)];
                        recweiz[it] = recz[I(it, itr)];
                    }
                }
                this->xcor(nt, 0, recweix, nt, 0, recweix, nt, 0, autocorrx);  /* for matrix */
                this->xcor(nt, 0, recweiy, nt, 0, recweiy, nt, 0, autocorry);  /* for matrix */
                this->xcor(nt, 0, recweiz, nt, 0, recweiz, nt, 0, autocorrz);  /* for matrix */
                this->xcor(nt, 0, recweix, nt, 0, modweix, nt, 0, crosscorrx); /* right hand side */
                this->xcor(nt, 0, recweiy, nt, 0, modweiy, nt, 0, crosscorry); /* right hand side */
                this->xcor(nt, 0, recweiz, nt, 0, modweiz, nt, 0, crosscorrz); /* right hand side */
                if (autocorrx[0] == 0.0)  autocorrx[0] = 1.0;
                if (autocorry[0] == 0.0)  autocorry[0] = 1.0;
                if (autocorrz[0] == 0.0)  autocorrz[0] = 1.0;
                autocorrx[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorry[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                autocorrz[0] *= (1.0 + PNOISE_ELASTIC);			/* whiten */
                this->stoep(nt, autocorrx, crosscorrx, shaperx, spikerx);
                this->stoep(nt, autocorry, crosscorry, shapery, spikery);
                this->stoep(nt, autocorrz, crosscorrz, shaperz, spikerz);
                xnorm = 0.0; 
                ynorm = 0.0; 
                znorm = 0.0; 
                for(it=0; it<nt; it++){
                    xnorm += shaperx[it]*shaperx[it];
                    ynorm += shapery[it]*shapery[it];
                    znorm += shaperz[it]*shaperz[it];
                }
                if (xnorm == 0.0) xnorm = 1.0;
                if (ynorm == 0.0) ynorm = 1.0;
                if (znorm == 0.0) znorm = 1.0;
                misfitx = 0.0;
                misfity = 0.0;
                misfitz = 0.0;
                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    misfitx += 0.5*(H*shaperx[it])*(H*shaperx[it])/xnorm;
                    misfity += 0.5*(H*shapery[it])*(H*shapery[it])/ynorm;
                    misfitz += 0.5*(H*shaperz[it])*(H*shaperz[it])/znorm;
                }

                for(it=0; it<nt; it++){
                    H = this->gauss(it, stdev);
                    resx[I(it, itr)] = -1.0*(H*H - 2.0*misfitx)*shaperx[it]/xnorm;
                    resy[I(it, itr)] = -1.0*(H*H - 2.0*misfity)*shapery[it]/ynorm;
                    resz[I(it, itr)] = -1.0*(H*H - 2.0*misfitz)*shaperz[it]/znorm;
                }
                this->stoep(nt, autocorrx, &resx[I(0, itr)], shaperx, spikerx);
                this->stoep(nt, autocorry, &resy[I(0, itr)], shapery, spikery);
                this->stoep(nt, autocorrz, &resz[I(0, itr)], shaperz, spikerz);
                this->convolve(nt, 0, shaperx, nt, 0, &recx[I(0, itr)], nt, 0, &resx[I(0, itr)]);        
                this->convolve(nt, 0, shapery, nt, 0, &recy[I(0, itr)], nt, 0, &resy[I(0, itr)]);        
                this->convolve(nt, 0, shaperz, nt, 0, &recz[I(0, itr)], nt, 0, &resz[I(0, itr)]);        
                if(dataweightxset)
                {
                    resx[I(it, itr)] *= weix[I(it, itr)];
                }
                if(dataweightyset)
                {
                    resy[I(it, itr)] *= weiy[I(it, itr)];
                }
                if(dataweightzset)
                {
                    resz[I(it, itr)] *= weiz[I(it, itr)];
                }
            }
            free(modweix);
            free(modweiy);
            free(modweiz);
            free(recweix);
            free(recweiy);
            free(recweiz);
            free(shaperx);
            free(spikerx);
            free(autocorrx);
            free(crosscorrx);
            free(shapery);
            free(spikery);
            free(autocorry);
            free(crosscorry);
            free(shaperz);
            free(spikerz);
            free(autocorrz);
            free(crosscorrz);
            break;

        default:
            for(itr=0; itr<ntr; itr++){
                for(it=1; it<nt; it++){
                   resx[I(it, itr)] = (modx[I(it, itr)] - recx[I(it, itr)]);
                   resy[I(it, itr)] = (mody[I(it, itr)] - recy[I(it, itr)]);
                   resz[I(it, itr)] = (modz[I(it, itr)] - recz[I(it, itr)]);
                   if(dataweightxset)
                   {
                       resx[I(it, itr)] *= weix[I(it, itr)]*weix[I(it, itr)];
                   }
                   if(dataweightyset)
                   {
                       resy[I(it, itr)] *= weiy[I(it, itr)]*weiy[I(it, itr)];
                   }
                   if(dataweightzset)
                   {
                       resz[I(it, itr)] *= weiz[I(it, itr)]*weiz[I(it, itr)];
                   }
                }
            }
            break;
    }
}


template<typename T>
int LsrtmElastic3D<T>::run(){
     int result = LSRTM_ERR;
     if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) {
         rs_warning("LsrtmElastic3D::run: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D_DS<T>> waves (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves->getNx_pml(), waves->getNy_pml(), waves->getNz_pml(), waves->getDx(), waves->getDy(), waves->getDz(), this->getOrder()));

     // Create snapshots
     std::shared_ptr<Snapshot3D<T>> Uxsnap;
     Uxsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uxsnap->openSnap(this->getSnapfile() + "-ux", 'w'); // Create a new snapshot file
     Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field

     std::shared_ptr<Snapshot3D<T>> Uysnap;
     Uysnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uysnap->openSnap(this->getSnapfile() + "-uy", 'w'); // Create a new snapshot file
     Uysnap->setData(waves->getUy1(), 0); //Set Uy as snap field

     std::shared_ptr<Snapshot3D<T>> Uzsnap;
     Uzsnap = std::make_shared<Snapshot3D<T>>(waves, this->getSnapinc());
     Uzsnap->openSnap(this->getSnapfile() + "-uz", 'w'); // Create a new snapshot file
     Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field

     this->writeLog("Running 3D Elastic full-waveform inversion gradient with full checkpointing.");
     this->writeLog("Doing forward Loop.");
    // Loop over forward time
    for(int it=0; it < nt; it++)
    {
    	//Writting out results to snapshot files
        Uxsnap->setData(waves->getUx1(), 0); //Set Ux as snap field
        Uxsnap->writeSnap(it);

        Uysnap->setData(waves->getUy1(), 0); //Set Uy as snap field
        Uysnap->writeSnap(it);

        Uzsnap->setData(waves->getUz1(), 0); //Set Uz as snap field
        Uzsnap->writeSnap(it);

    	// Time stepping stress
    	waves->forwardstepStress(model, der);

    	// Inserting pressure source 
    	waves->insertPressuresource(model, source, SMAP, it);

    	// Time stepping displacement
    	waves->forwardstepDisplacement(model, der);
    
    	// Inserting pressure source 
    	waves->insertForcesource(model, source, SMAP, it);

        // Recording data (Ux)
        if(this->datamodUxset){
            waves->recordData(model, this->datamodUx, GMAP, it);
        }

        // Recording data (Uy)
        if(this->datamodUyset){
            waves->recordData(model, this->datamodUy, GMAP, it);
        }

        // Recording data (Uz)
        if(this->datamodUzset){
            waves->recordData(model, this->datamodUz, GMAP, it);
        }

        // Roll pointers
        waves->roll();

        // Output progress to logfile
        this->writeProgress(it, nt-1, 20, 48);
    }//End of forward loop
    
    
    //Close snapshot file
    Uxsnap->closeSnap();
    Uysnap->closeSnap();
    Uzsnap->closeSnap();

    // Reset waves
    waves.reset();
    waves  = std::make_shared<WavesElastic3D_DS<T>>(model, nt, dt, ot);

    // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

    Uxsnap->openSnap(this->getSnapfile() + "-ux", 'r');
    Uxsnap->allocSnap(0);

    Uysnap->openSnap(this->getSnapfile() + "-uy", 'r');
    Uysnap->allocSnap(0);

    Uzsnap->openSnap(this->getSnapfile() + "-uz", 'r');
    Uzsnap->allocSnap(0);

    if(!this->getNoreverse()){
        // Compute misfit
        computeMisfit();

        // Compute Residuals
        computeResiduals();

        this->writeLog("\nDoing reverse-time Loop.");
        // Loop over reverse time
        for(int it=0; it < nt; it++)
        {
            // Time stepping stress
            waves->forwardstepStress(model, der);

            // Inserting residuals
            waves->insertPressuresource(model, dataresUx, GMAP, (nt - 1 - it));
            waves->insertPressuresource(model, dataresUy, GMAP, (nt - 1 - it));
            waves->insertPressuresource(model, dataresUz, GMAP, (nt - 1 - it));

            // Time stepping displacement
            waves->forwardstepDisplacement(model, der);

            // Inserting residuals
            waves->insertForcesource(model, dataresUx, GMAP, (nt - 1 - it));
            waves->insertForcesource(model, dataresUy, GMAP, (nt - 1 - it));
            waves->insertForcesource(model, dataresUz, GMAP, (nt - 1 - it));

            //Read forward snapshot
            Uxsnap->readSnap(nt - 1 - it);
            Uysnap->readSnap(nt - 1 - it);
            Uzsnap->readSnap(nt - 1 - it);

            // Do Crosscorrelation
            if((((nt - 1 - it)-Uxsnap->getEnddiff()) % Uxsnap->getSnapinc()) == 0){
                crossCorr(Uxsnap->getData(0), Uysnap->getData(0), Uzsnap->getData(0), 0, waves, model, (nt - 1 - it));
            }

            // Roll pointers
            waves->roll();

            // Output progress to logfile
            this->writeProgress(it, nt-1, 20, 48);
        }
        this->writeLog("\nGradient computation completed.");
    } // End of reverse loop

    // Scale gradients
    this->scaleGrad(model);
    
	//Remove snapshot file
	Uxsnap->removeSnap();
	Uysnap->removeSnap();
	Uzsnap->removeSnap();

    result=LSRTM_OK;
    return result;
}

template<typename T>
int LsrtmElastic3D<T>::run_optimal(){
     int result = LSRTM_ERR;
     if(!vpgradset && !vsgradset && !rhogradset && !wavgradset) {
         rs_warning("LsrtmElastic3D::run: No image set");
         return result;
     }
     int nt;
     float dt;
	 float ot;

     nt = source->getNt();
     dt = source->getDt();
     ot = source->getOt();

	// Create log file
     this->createLog(this->getLogfile());

     // Create the classes 
     std::shared_ptr<WavesElastic3D_DS<T>> waves_fw (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<WavesElastic3D_DS<T>> waves_bw (new WavesElastic3D_DS<T>(model, nt, dt, ot));
     std::shared_ptr<Der<T>> der (new Der<T>(waves_fw->getNx_pml(), waves_fw->getNy_pml(), waves_fw->getNz_pml(), waves_fw->getDx(), waves_fw->getDy(), waves_fw->getDz(), this->getOrder()));
     std::shared_ptr<Revolve<T>> optimal (new Revolve<T>(nt, this->getNcheck(), this->getIncore()));
     revolve_action whatodo;
     int oldcapo,capo;
     capo = 0;

     // Create checkpoint file
     optimal->openCheck(this->getSnapfile(), waves_fw, 'w');

     // Create image
    if(this->vpgradset) vpgrad->allocateImage();
    if(this->vsgradset) vsgrad->allocateImage();
    if(this->rhogradset) rhograd->allocateImage();

     this->writeLog("Running 3D Elastic full-waveform inversion gradient with optimal checkpointing.");
     this->writeLog("Doing forward Loop.");
     bool reverse = false;
    // Loop over forward time
    do
    {
        oldcapo=optimal->getCapo();
        whatodo = optimal->revolve();
        capo = optimal->getCapo();
        if (whatodo == advance)
        {
            for(int it=oldcapo; it < capo; it++)
            {
                // Time stepping stress
                waves_fw->forwardstepStress(model, der);

                // Inserting source 
                waves_fw->insertPressuresource(model, source, SMAP, it);

                // Time stepping displacement
                waves_fw->forwardstepDisplacement(model, der);

                // Inserting source 
                waves_fw->insertForcesource(model, source, SMAP, it);

                // Recording data (Ux)
                if(this->datamodUxset && !reverse){
                    waves_fw->recordData(model, this->datamodUx, GMAP, it);
                }

                // Recording data (Uy)
                if(this->datamodUyset && !reverse){
                    waves_fw->recordData(model, this->datamodUy, GMAP, it);
                }

                // Recording data (Uz)
                if(this->datamodUzset && !reverse){
                    waves_fw->recordData(model, this->datamodUz, GMAP, it);
                }

                // Roll pointers
                waves_fw->roll();

                if(!reverse){
                    // Output progress to logfile
                    this->writeProgress(it, nt-1, 20, 48);
                }
            }
        }
        if (whatodo == firsturn)
        {
            // Time stepping stress
            waves_fw->forwardstepStress(model, der);

            // Inserting source 
            waves_fw->insertPressuresource(model, source, SMAP, capo);

            // Time stepping displacement
            waves_fw->forwardstepDisplacement(model, der);

            // Inserting source 
            waves_fw->insertForcesource(model, source, SMAP, capo);

            // Recording data (Ux)
            if(this->datamodUxset){
                waves_fw->recordData(model, this->datamodUx, GMAP, capo);
            }

            // Recording data (Uy)
            if(this->datamodUyset){
                waves_fw->recordData(model, this->datamodUy, GMAP, capo);
            }

            // Recording data (Uz)
            if(this->datamodUzset){
                waves_fw->recordData(model, this->datamodUz, GMAP, capo);
            }

            // Compute misfit
            computeMisfit();

            // Compute Residuals
            computeResiduals();

            // Inserting residuals
            waves_bw->insertPressuresource(model, dataresUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUy, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUz, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUy, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUz, GMAP, capo);

            // Do Crosscorrelation 
            T *wsx = waves_fw->getUx1();
            T *wsy = waves_fw->getUy1();
            T *wsz = waves_fw->getUz1();
            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), waves_bw, model, capo);

            // Roll pointers
            waves_fw->roll();
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(capo, nt-1, 20, 48);
      
            //Close checkpoint file for w and reopen for rw
            optimal->closeCheck();
            optimal->openCheck(this->getSnapfile(), waves_fw, 'a');
            reverse = true;
            // Output progress to logfile
            this->writeLog("\nDoing reverse-time Loop.");
            this->writeProgress(0, nt-1, 20, 48);
        }
        if (whatodo == youturn)
        {
            // Time stepping stress
            waves_bw->forwardstepStress(model, der);

            // Inserting residuals
            waves_bw->insertPressuresource(model, dataresUx, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUy, GMAP, capo);
            waves_bw->insertPressuresource(model, dataresUz, GMAP, capo);

            // Time stepping displacement
            waves_bw->forwardstepDisplacement(model, der);

            // Inserting residuals
            waves_bw->insertForcesource(model, dataresUx, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUy, GMAP, capo);
            waves_bw->insertForcesource(model, dataresUz, GMAP, capo);

            // Do Crosscorrelation
            T *wsx = waves_fw->getUx1();
            T *wsy = waves_fw->getUy1();
            T *wsz = waves_fw->getUz1();
            crossCorr(wsx, wsy, wsz, waves_fw->getLpml(), waves_bw, model, capo);
            
            // Roll pointers
            waves_bw->roll();

            // Output progress to logfile
            this->writeProgress(nt-1-capo, nt-1, 20, 48);
        }
        if (whatodo == takeshot)
        {
            optimal->writeCheck(waves_fw);
        }
        if (whatodo == restore)
        {
            optimal->readCheck(waves_fw);
        }

        if(whatodo == error){
            std::cerr << "Error!" << std::endl;
        }

    } while((whatodo != terminate) && (whatodo != error));
     this->writeLog("\nGradient computation completed.");

    // Scale gradients
    this->scaleGrad(model);

	//Remove snapshot file
	optimal->removeCheck();

    result=LSRTM_OK;
    return result;
}

template<typename T>
LsrtmElastic3D<T>::~LsrtmElastic3D() {
    // Nothing here
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Lsrtm<float>;
template class Lsrtm<double>;
template class LsrtmAcoustic2D<float>;
template class LsrtmAcoustic2D<double>;
template class LsrtmAcoustic3D<float>;
template class LsrtmAcoustic3D<double>;

template class LsrtmElastic2D<float>;
template class LsrtmElastic2D<double>;
template class LsrtmElastic3D<float>;
template class LsrtmElastic3D<double>;

}
