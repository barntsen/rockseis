// Include statements
#include "rays.h"

#define I2D(i,j) ((i)*nx + (j)) 
#define I3D(i,j,k) ((k)*nx*ny+(j)*nx + (i)) 
#define S(i) (1.0/Vp[i]) 

namespace rockseis {

// =============== ABSTRACT RAYS CLASS =============== //
template<typename T>
Rays<T>::Rays() {
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, 1);
    geometry->setN(2, 1);
    geometry->setN(3, 1);
    geometry->setD(1, 1.);
    geometry->setD(2, 1.);
    geometry->setD(3, 1.);
    geometry->setO(1, 0.);
    geometry->setO(2, 0.);
    geometry->setO(3, 0.);
    dim=0;
    lpml=10;
}

template<typename T>
Rays<T>::~Rays() {
    // Nothing here
}

template<typename T>
Rays<T>::Rays(const int _dim, const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz) {
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, _nx);
    geometry->setN(2, _ny);
    geometry->setN(3, _nz);
    geometry->setD(1, _dx);
    geometry->setD(2, _dy);
    geometry->setD(3, _dz);
    geometry->setO(1, _ox);
    geometry->setO(2, _oy);
    geometry->setO(3, _oz);
    dim = _dim;
    lpml = _lpml;

}

// =============== 2D ACOUSTIC RAYS CLASS =============== //
template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(){
    int nx, nz;
    nx = this->getNx_pml();
    nz = this->getNz_pml();

    /* Allocate memory variables */
    TT = (T *) calloc(nx*nz, sizeof(T));
    lam = (T *) calloc(nx*nz, sizeof(T));
    recmask = (bool *) calloc(nx*nz, sizeof(bool));

    /* Initialize arrays */
    for (int i=0; i < nx*nz; i++){
        TT[i] = 2.0*TMAX;
    }

    Index Ilam(nx,nz);
    for (int i=1; i < nx-1; i++){
        for (int j=1; j < nz-1; j++){
            lam[Ilam(i,j)] = 10.0*TMAX;
        }
    }
}

template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz): Rays<T>(2, _nx, 1, _nz, _lpml, _dx, 1.0, _dz, _ox, 0.0, _oz) {

    int nx, nz;
    nx = _nx + 2*_lpml;
    nz = _nz + 2*_lpml;
    /* Allocate memory variables */
    TT = (T *) calloc(nx*nz, sizeof(T));
    lam = (T *) calloc(nx*nz, sizeof(T));
    recmask = (bool *) calloc(nx*nz, sizeof(bool));

    /* Initialize TT */
    for (int i=0; i < nx*nz; i++){
        TT[i] = 2.0*TMAX;
    }
    Index Ilam(nx,nz);
    for (int i=1; i < nx-1; i++){
        for (int j=1; j < nz-1; j++){
            lam[Ilam(i,j)] = 10.0*TMAX;
        }
    }
}


template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(std::shared_ptr<rockseis::ModelEikonal2D<T>> _model): Rays<T>(){

    int _nx, _ny, _nz, _lpml;
    int nx_pml, nz_pml;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim;

    /* Get necessary parameters from model class */
    _nx=_model->getNx();
    _ny=_model->getNy();
    _nz=_model->getNz();
    _dx=_model->getDx();
    _dy=_model->getDy();
    _dz=_model->getDz();
    _ox=_model->getOx();
    _oy=_model->getOy();
    _oz=_model->getOz();
    _dim = _model->getDim();
    _lpml = _model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setDim(_dim);
    this->setLpml(_lpml);

    // Setting model pointer
    model = _model;

    nx_pml= _nx + 2*_lpml;
    nz_pml= _nz + 2*_lpml;
    /* Allocate memory variables */
    TT = (T *) calloc(nx_pml*nz_pml, sizeof(T));
    lam = (T *) calloc(nx_pml*nz_pml, sizeof(T));
    recmask = (bool *) calloc(nx_pml*nz_pml, sizeof(bool));

    /* Initialize TT */
    for (int i=0; i < nx_pml*nz_pml; i++){
        TT[i] = 2.0*TMAX;
    }
    Index Ilam(nx_pml,nz_pml);
    for (int i=1; i < nx_pml-1; i++){
        for (int j=1; j < nz_pml-1; j++){
            lam[Ilam(i,j)] = 10.0*TMAX;
        }
    }

}

template<typename T>
RaysAcoustic2D<T>::~RaysAcoustic2D() {
    /* Free allocated variables */
    free(TT);
    free(lam);
    free(recmask);
}

template<typename T>
void RaysAcoustic2D<T>::sweep(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy)
{ 
/*------------------------------------------------------------------------
 *  Solve eikonal equation by fast sweeping method according to Zhao (2004)
 *
 *  D. Koehn
 *  Kiel, 09/12/2015
 *  ----------------------------------------------------------------------*/

    int nx = model->getNx_pml();
    int ny = model->getNz_pml();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *Vp = model->getL();

    /* local variables */
    int i, j, h, k;
    T a, b, Tt;


    /* sweep over FD-grid */
    h = nx1;
    for (i=0;i<nx;i++){
        k = ny1;
        for (j=0;j<ny;j++){	

            /* model interior */
            if((h>0)&&(h<nx-1)&&(k>0)&&(k<ny-1)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* model borders */
            /* left */
            if((h==0)&&(k>0)&&(k<ny-1)){
                a = MIN(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* right */
            if((h==nx-1)&&(k>0)&&(k<ny-1)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* top */
            if((k==0)&&(h>0)&&(h<nx-1)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }

            /* bottom */
            if((k==ny-1)&&(h>0)&&(h<nx-1)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* model corners */
            /* upper-left */
            if((h==0)&&(k==0)){
                a = MIN(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }

            /* lower-left */
            if((h==0)&&(k==ny-1)){
                a = MIN(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* lower-right */
            if((h==nx-1)&&(k==ny-1)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = MIN(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* upper-right */
            if((h==nx-1)&&(k==0)){
                a = MIN(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = MIN(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }


            /* calculate solution */
            if(ABS(a-b)>=(S(I2D(k,h))*dh)){
                Tt = MIN(a,b) + S(I2D(k,h))*dh;
            }else{
                Tt = (a + b + sqrt((2.0*pow(S(I2D(k,h)),2.0)*pow(dh,2.0)) - pow((a-b),2.0)))/2.0;
            }

            TT[I2D(k,h)] = MIN(TT[I2D(k,h)],Tt);

            k += ndy;
        }
        h += ndx;
    }
}

template<typename T>
void RaysAcoustic2D<T>::sweep_adj(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy)
{
    int nx = model->getNx_pml();
    int ny = model->getNz_pml();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *lam = this->getLam();
    bool *recmask = this->getRecmask();

    /* local variables */
    int i, j, h, k;
    T app, amp, apm, amm;
    T bpp, bmp, bpm, bmm;
    T ap, am , bp, bm;
    T lhs, rhs, lamt;

    /* sweep over FD-grid */
    h = nx1;
    for (i=1;i<nx-1;i++){
        k = ny1;
        for (j=1;j<ny-1;j++){
            if(recmask[I2D(k,h)]==0){
                /* assemble equation (3.6) in Leung & Qian (2006) */
                ap = -(TT[I2D(k,h+1)]-TT[I2D(k,h)])/dh;
                am = -(TT[I2D(k,h)]-TT[I2D(k,h-1)])/dh;

                bp = -(TT[I2D(k+1,h)]-TT[I2D(k,h)])/dh;
                bm = -(TT[I2D(k,h)]-TT[I2D(k-1,h)])/dh;

                app = (ap + ABS(ap))/2.0;
                apm = (ap - ABS(ap))/2.0;

                amp = (am + ABS(am))/2.0;
                amm = (am - ABS(am))/2.0;

                bpp = (bp + ABS(bp))/2.0;
                bpm = (bp - ABS(bp))/2.0;

                bmp = (bm + ABS(bm))/2.0;
                bmm = (bm - ABS(bm))/2.0;


                /* Leung & Qian (2006) */
                lhs = (app-amm)/dh + (bpp-bmm)/dh;
                rhs = (amp*lam[I2D(k,h-1)]-apm*lam[I2D(k,h+1)])/dh + (bmp*lam[I2D(k-1,h)]-bpm*lam[I2D(k+1,h)])/dh;
                //
                lamt = rhs*lhs/(SQ(lhs)+EPS_ADJ);
                lam[I2D(k,h)] = MIN(lam[I2D(k,h)],lamt);

            }

            k += ndy;
        }
        h += ndx;
    }
}

template<typename T>
void RaysAcoustic2D<T>::clearTT() {
    /* Re-initialize TT array */

    int nx = model->getNx_pml();
    int ny = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny; i++){
        TT[i] = 2.0*TMAX;
    }
}

template<typename T>
void RaysAcoustic2D<T>::clearLam() {
    /* Re-initialize lam array */

    int nx = model->getNx_pml();
    int ny = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny; i++){
        lam[i] = 10.0*TMAX;
    }
}

template<typename T>
int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
    Point2D<int> *map;
    int ntrace = source->getNtrace();
    int nx, nz;
    int lpml;
    int nr = 0;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    int i;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0)
        { 
            TT[I(lpml+map[i].x, lpml+map[i].y)] = 0.0;
            nr++;
        }
    }
    return nr;
}

template<typename T>
int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int traceno){
    Point2D<int> *map;
    int ntrace = source->getNtrace();
    if(traceno < 0 || traceno > ntrace-1){
        rs_error("RaysAcoustic2D<T>::insertSource: traceno out of bounds.");
    }
    int nx, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();
    int nr = 0;

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    int i = traceno;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    if(map[i].x >= 0 && map[i].y >=0)
    { 
        TT[I(lpml+map[i].x, lpml+map[i].y)] = 0.0;
        nr++;
    }
    return nr;
}

template<typename T>
void RaysAcoustic2D<T>::recordData(std::shared_ptr<rockseis::Data2D<T>> data, bool maptype){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *dataarray; 
    T *Fielddata;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    dataarray = data->getData();
    int i;
    Index I(nx, nz);
    Index Idat(nt, ntrace);
    Fielddata = this->getTT();
    for (i=0; i < ntrace; i++) 
    { 
        if(map[i].x >= 0 && map[i].y >=0)
        {
            dataarray[Idat(0,i)] = Fielddata[I(lpml+map[i].x, lpml+map[i].y)];
        }
    }

}

template<typename T>
T RaysAcoustic2D<T>::norm1(T *TT, T *TTold)
{
    /* local variables */
    int NX = this->getNx_pml();
    int NY = this->getNz_pml();

    int i;
    T sum, norml1;

    /* estimate L1 norm */
    sum = 0.0;
    for (i=0;i<NX*NY;i++){
        sum += (TT[i]-TTold[i]);
    }

    norml1 = sum/(NX*NY);
    norml1 = ABS(norml1);

    return norml1;
}

template<typename T>
void RaysAcoustic2D<T>::createRecmask(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
    Point2D<int> *map;
    int ntrace = source->getNtrace();
    int nx, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    int i;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    Index Idat(1, ntrace); // Data indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0)
        { 
            recmask[I(lpml+map[i].x, lpml+map[i].y)] = true;
        }
    }
}

template<typename T>
void RaysAcoustic2D<T>::insertResiduals(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
    Point2D<int> *map;
    int ntrace = source->getNtrace();
    int nx, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    T *res;
    res = source->getData();

    int i;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    Index Idat(1, ntrace); // Data indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0)
        { 
            lam[I(lpml+map[i].x, lpml+map[i].y)] = res[Idat(0,i)];
        }
    }
}

template<typename T>
void RaysAcoustic2D<T>::solve()
{
    int nx = this->getNx_pml();
    int nz = this->getNz_pml();
    T tmax = TMAX;
    T *TTold = (T *) calloc(nx*nz, sizeof(T));
    int i;
    int iter=0;

    /* initialize l1 norm of TT - TTold */
    T lnorm1 = 10.0 * tmax;

    /* apply fast sweeping method to solve the eikonal equation */
    while(lnorm1>=TTNORM){

        /* save old TT values */
        for (i=0; i<nx*nz; i++){
            TTold[i] = TT[i];
        }

        /* sweep with order according to Zhao (2004) */
        sweep(0,nx-1,1,0,nz-1,1);
        sweep(nx-1,0,-1,0,nz-1,1);
        sweep(nx-1,0,-1,nz-1,0,-1);
        sweep(0,nx-1,1,nz-1,0,-1);

        /* calculate l1 norm of TT - TTold */
        lnorm1 = norm1(TT,TTold);
        iter++;
        if(iter > MAXITER) rs_error("RaysAcoustic2D<T>::solve(): Number of iterations exceeded limit. Eikonal solver is not converging. Try increasing the smoothening of the model."); 
    }
    free(TTold);
}

template<typename T>
void RaysAcoustic2D<T>::solve_adj()
{
    int nx = this->getNx_pml();
    int nz = this->getNz_pml();
    T tmax = TMAX;
    T *lamold = (T *) calloc(nx*nz, sizeof(T));
    int i;
    int iter=0;

    /* initialize l1 norm of lam - lamold */
    T lnorm1 = 10.0 * tmax;

    /* apply fast sweeping method to solve the eikonal equation */
    while(lnorm1>=TTNORM){

        /* save old lam values */
        for (i=0; i<nx*nz; i++){
            lamold[i] = lam[i];
        }

        /* sweep with order according to Zhao (2004) */
		sweep_adj(1,nx-2,1,1,nz-2,1);
		sweep_adj(1,nx-2,1,nz-2,1,-1);
		sweep_adj(nx-2,1,-1,1,nz-2,1);
		sweep_adj(nx-2,1,-1,nz-2,1,-1);

		/* calculate l1 norm of lam - lamold */
        lnorm1 = norm1(lam,lamold);

        iter++;
        if(iter > MAXITER) rs_error("RaysAcoustic2D<T>::solve_adj(): Number of iterations exceeded limit. Eikonal solver is not converging. Try increasing the smoothening of the model."); 
    }
    free(lamold);
}

// =============== 3D ACOUSTIC RAYS CLASS =============== //
template<typename T>
RaysAcoustic3D<T>::RaysAcoustic3D(){
    size_t nx, ny, nz;
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    /* Allocate memory variables */
    TT = (T *) calloc(nx*ny*nz, sizeof(T));
    lam = (T *) calloc(nx*ny*nz, sizeof(T));
    recmask = (bool *) calloc(nx*ny*nz, sizeof(bool));

    /* Initialize arrays */
    for (size_t i=0; i < nx*ny*nz; i++){
        TT[i] = 2.0*TMAX;
    }

    Index Ilam(nx,ny,nz);
    for (size_t i=1; i < nx-1; i++){
        for (size_t j=1; j < ny-1; j++){
            for (size_t k=1; k < nz-1; k++){
                lam[Ilam(i,j,k)] = 10.0*TMAX;
            }
        }
    }
}

template<typename T>
RaysAcoustic3D<T>::RaysAcoustic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz): Rays<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz) {

    size_t nx, ny, nz;
    nx = _nx + 2*_lpml;
    ny = _ny + 2*_lpml;
    nz = _nz + 2*_lpml;
    /* Allocate memory variables */
    TT = (T *) calloc(nx*ny*nz, sizeof(T));
    lam = (T *) calloc(nx*ny*nz, sizeof(T));
    recmask = (bool *) calloc(nx*ny*nz, sizeof(bool));

    /* Initialize TT */
    for (size_t i=0; i < nx*ny*nz; i++){
        TT[i] = 2.0*TMAX;
    }
    Index Ilam(nx,ny,nz);
    for (size_t i=1; i < nx-1; i++){
        for (size_t j=1; j < ny-1; j++){
            for (size_t k=1; k < nz-1; k++){
                lam[Ilam(i,j,k)] = 10.0*TMAX;
            }
        }
    }
}


template<typename T>
RaysAcoustic3D<T>::RaysAcoustic3D(std::shared_ptr<rockseis::ModelEikonal3D<T>> _model): Rays<T>(){

    size_t _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim;
    int _lpml;

    /* Get necessary parameters from model class */
    _nx=_model->getNx();
    _ny=_model->getNy();
    _nz=_model->getNz();
    _dx=_model->getDx();
    _dy=_model->getDy();
    _dz=_model->getDz();
    _ox=_model->getOx();
    _oy=_model->getOy();
    _oz=_model->getOz();
    _dim = _model->getDim();
    _lpml = _model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setDim(_dim);
    this->setLpml(_lpml);

    // Setting model pointer
    model = _model;

    size_t nx, ny, nz;
    nx = _nx + 2*_lpml;
    ny = _ny + 2*_lpml;
    nz = _nz + 2*_lpml;
    /* Allocate memory variables */
    TT = (T *) calloc(nx*ny*nz, sizeof(T));
    lam = (T *) calloc(nx*ny*nz, sizeof(T));
    recmask = (bool *) calloc(nx*ny*nz, sizeof(bool));

    /* Initialize TT */
    for (size_t i=0; i < nx*ny*nz; i++){
        TT[i] = 2.0*TMAX;
    }
    Index Ilam(nx,ny,nz);
    for (size_t i=1; i < nx-1; i++){
        for (size_t j=1; j < ny-1; j++){
            for (size_t k=1; k < nz-1; k++){
                lam[Ilam(i,j,k)] = 10.0*TMAX;
            }
        }
    }
}

template<typename T>
RaysAcoustic3D<T>::~RaysAcoustic3D() {
    /* Free allocated variables */
    free(TT);
    free(lam);
    free(recmask);
}

template<typename T>
void RaysAcoustic3D<T>::sweep(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int nz1, int nz2, int ndz)
{ 
    size_t nx = model->getNx_pml();
    size_t ny = model->getNy_pml();
    size_t nz = model->getNz_pml();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *Vp = model->getL();

    /* local variables */
    size_t i, j, k, m, n, l;
    T a[3];
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
    T Tt;

    T slo;

    /* sweep over FD-grid */
    m = nx1;
    for (i=0;i<nx;i++){
        n = ny1;
        for (j=0;j<ny;j++){	
            l = nz1;
            for (k=0;k<nz;k++){	

                /* model interior */
                if((m>0)&&(m<nx-1)&&(n>0)&&(n<ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

             
                /* model borders */
                /* left */
                if((m==0)&&(n>0)&&(n<ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* right */
                if((m==nx-1)&&(n>0)&&(n<ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* top */
                if((m>0)&&(m<nx-1)&&(n>0)&&(n<ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* bottom */
                if((m>0)&&(m<nx-1)&&(n>0)&&(n<ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }


                /* front */
                if((m>0)&&(m<nx-1)&&(n==0)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* back */
                if((m>0)&&(m<nx-1)&&(n==ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }


                /* model edges */
                /* upper-front */
                if((m>0)&&(m<nx-1)&&(n==0)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* upper-back */
                if((m>0)&&(m<nx-1)&&(n==ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* lower-front */
                if((m>0)&&(m<nx-1)&&(n==0)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* lower-back */
                if((m>0)&&(m<nx-1)&&(n==ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* left-front */
                if((m==0)&&(n==0)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* left-back */
                if((m==0)&&(n==ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* right-front */
                if((m==nx-1)&&(n==0)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* right-back */
                if((m==nx-1)&&(n==ny-1)&&(l>0)&&(l<nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);
                }

                /* upper-left */
                if((m==0)&&(n>0)&&(n<ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                 /* upper-right */
                if((m==nx-1)&&(n>0)&&(n<ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* lower-left */
                if((m==0)&&(n>0)&&(n<ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                 /* lower-right */
                if((m==nx-1)&&(n>0)&&(n<ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* model corners */
                /* upper-left-front */
                if((m==0)&&(n==0)&&(l==0)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* upper-left-back */
                if((m==0)&&(n==ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* upper-right-front */
                if((m==nx-1)&&(n==0)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* upper-right-back */
                if((m==nx-1)&&(n==ny-1)&&(l==0)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n,l+1)]);
                }

                /* lower-left-front */
                if((m==0)&&(n==0)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* lower-left-back */
                if((m==0)&&(n==ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m,n,l)],TT[I3D(m+1,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* lower-right-front */
                if((m==nx-1)&&(n==0)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n,l)],TT[I3D(m,n+1,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }

                /* lower-right-back */
                if((m==nx-1)&&(n==ny-1)&&(l==nz-1)){
                    a[0] = MIN(TT[I3D(m-1,n,l)],TT[I3D(m,n,l)]);
                    a[1] = MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n,l)]);
                    a[2] = MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l)]);
                }


                /* calculate solution */
                std::sort(a,a+3);
                while(1){
                    slo = S(I3D(m,n,l))*dh;
                    Tt = a[0] + slo;
                    if(Tt <= a[1]) break;
                    Tt = 0.5*(a[0]+a[1]+sqrt(2.*slo*slo - (a[0]-a[1])*(a[0]-a[1])));
                    if(Tt <= a[2]) break;
                    Tt = 1./3. * ((a[0] + a[1] + a[2]) + sqrt(-2.*a[0]*a[0] + 2.*a[0]*a[1] - 2.*a[1]*a[1] + 2.*a[0]*a[2] + 2.*a[1]*a[2] - 2.*a[2]*a[2] + 3.*slo*slo));
                    break;
                }


                TT[I3D(m,n,l)] = MIN(TT[I3D(m,n,l)],Tt);
                l += ndz;
            }
            n += ndy;
        }
        m += ndx;
    }
}

template<typename T>
void RaysAcoustic3D<T>::sweep_adj(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int nz1, int nz2, int ndz)
{
    int nx = model->getNx_pml();
    int ny = model->getNy_pml();
    int nz = model->getNz_pml();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *lam = this->getLam();
    bool *recmask = this->getRecmask();

    /* local variables */
    int i, j, k, m, n, l;
    T app, amp, apm, amm;
    T bpp, bmp, bpm, bmm;
    T cpp, cmp, cpm, cmm;
    T ap, am , bp, bm, cp, cm;
    T lhs, rhs, lamt;

    /* sweep over FD-grid */
    m = nx1;
    for (i=1;i<nx-1;i++){
        n = ny1;
        for (j=1;j<ny-1;j++){
            l = nz1;
            for (k=1;k<nz-1;k++){
                if(recmask[I3D(m,n,l)]==0){
                    /* assemble equation (3.6) in Leung & Qian (2006) */
                    ap = -(TT[I3D(m+1,n,l)]-TT[I3D(m,n,l)])/dh;
                    am = -(TT[I3D(m,n,l)]-TT[I3D(m-1,n,l)])/dh;

                    bp = -(TT[I3D(m,n+1,l)]-TT[I3D(m,n,l)])/dh;
                    bm = -(TT[I3D(m,n,l)]-TT[I3D(m,n-1,l)])/dh;

                    cp = -(TT[I3D(m,n,l+1)]-TT[I3D(m,n,l)])/dh;
                    cm = -(TT[I3D(m,n,l)]-TT[I3D(m,n,l-1)])/dh;

                    app = (ap + ABS(ap))/2.0;
                    apm = (ap - ABS(ap))/2.0;

                    amp = (am + ABS(am))/2.0;
                    amm = (am - ABS(am))/2.0;

                    bpp = (bp + ABS(bp))/2.0;
                    bpm = (bp - ABS(bp))/2.0;

                    bmp = (bm + ABS(bm))/2.0;
                    bmm = (bm - ABS(bm))/2.0;

                    cpp = (cp + ABS(cp))/2.0;
                    cpm = (cp - ABS(cp))/2.0;

                    cmp = (cm + ABS(cm))/2.0;
                    cmm = (cm - ABS(cm))/2.0;


                    /* Leung & Qian (2006) */
                    lhs = (app-amm)/dh + (bpp-bmm)/dh + (cpp-cmm)/dh;
                    rhs = (amp*lam[I3D(m-1,n,l)]-apm*lam[I3D(m+1,n,l)])/dh + (bmp*lam[I3D(m,n-1,l)]-bpm*lam[I3D(m,n+1,l)])/dh + (cmp*lam[I3D(m,n,l-1)]-cpm*lam[I3D(m,n,l+1)])/dh;

                    lamt = rhs*lhs/(SQ(lhs)+EPS_ADJ);
                    lam[I3D(m,n,l)] = MIN(lam[I3D(m,n,l)],lamt);
                }
                l += ndz;
            }

            n += ndy;
        }
        m += ndx;
    }
}


template<typename T>
void RaysAcoustic3D<T>::solve()
{
    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();

    T tmax = TMAX;
    T *TTold = (T *) calloc(nx*ny*nz, sizeof(T));
    size_t i;
    int iter=0;

    /* initialize l1 norm of TT - TTold */
    T lnorm1 = 10.0 * tmax;

    /* apply fast sweeping method to solve the eikonal equation */
    while(lnorm1>=TTNORM){

        /* save old TT values */
        for (i=0; i<nx*ny*nz; i++){
            TTold[i] = TT[i];
        }

        sweep(0, nx-1, 1,  0, ny-1, 1,  0, nz-1, 1);
        sweep(nx-1, 0, -1,  0, ny-1, 1,  0, nz-1, 1);
        sweep(0, nx-1, 1,  ny-1, 0, -1,  0, nz-1, 1);
        sweep(nx-1, 0, -1,  ny-1, 0, -1,  0, nz-1, 1);
        sweep(0, nx-1, 1,  0, ny-1, 1,  nz-1, 0, -1);
        sweep(nx-1, 0, -1,  0, ny-1, 1,  nz-1, 0, -1);
        sweep(0, nx-1, 1,  ny-1, 0, -1,  nz-1, 0, -1);
        sweep(nx-1, 0, -1,  ny-1, 0, -1,  nz-1, 0, -1);

        /* calculate l1 norm of TT - TTold */
        lnorm1 = norm1(TT,TTold);
        iter++;
        if(iter > MAXITER) rs_error("RaysAcoustic3D<T>::solve(): Number of iterations exceeded limit. Eikonal solver is not converging. Try increasing the smoothening of the model."); 
    }
    free(TTold);
}

template<typename T>
void RaysAcoustic3D<T>::solve_adj()
{
    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    T tmax = TMAX;
    T *lamold = (T *) calloc(nx*ny*nz, sizeof(T));
    size_t i;
    int iter=0;

    /* initialize l1 norm of lam - lamold */
    T lnorm1 = 10.0 * tmax;

    /* apply fast sweeping method to solve the eikonal equation */
    while(lnorm1>=TTNORM){

        /* save old lam values */
        for (i=0; i<nx*ny*nz; i++){
            lamold[i] = lam[i];
        }

        /* sweep with order according to Zhao (2004) */
        sweep_adj(1, nx-2, 1,  1, ny-2, 1,  1, nz-2, 1);
        sweep_adj(1, nx-2, 1,  1, ny-2, 1,  nz-2, 1, -1);
        sweep_adj(1, nx-2, 1,  ny-2, 1, -1,  1, nz-2, 1);
        sweep_adj(1, nx-2, 1,  ny-2, 1, -1,  nz-2, 1, -1);
        sweep_adj(nx-2, 1, -1,  1, ny-2, 1,  1, nz-2, 1);
        sweep_adj(nx-2, 1, -1,  1, ny-2, 1,  nz-2, 1, -1);
        sweep_adj(nx-2, 1, -1,  ny-2, 1, -1,  1, nz-2, 1);
        sweep_adj(nx-2, 1, -1,  ny-2, 1, -1,  nz-2, 1, -1);

		/* calculate l1 norm of lam - lamold */
        lnorm1 = norm1(lam,lamold);

        iter++;
        if(iter > MAXITER) rs_error("RaysAcoustic3D<T>::solve_adj(): Number of iterations exceeded limit. Eikonal solver is not converging. Try increasing the smoothening of the model."); 
    }
    free(lamold);
}

template<typename T>
T RaysAcoustic3D<T>::norm1(T *TT, T *TTold)
{
    /* local variables */
    size_t NX = this->getNx_pml();
    size_t NY = this->getNy_pml();
    size_t NZ = this->getNz_pml();

    size_t i;
    T sum, norml1;

    /* estimate L1 norm */
    sum = 0.0;
    for (i=0;i<NX*NY*NZ;i++){
        sum += (TT[i]-TTold[i]);
    }

    norml1 = sum/(NX*NY*NZ);
    norml1 = ABS(norml1);

    return (norml1);
}

template<typename T>
int RaysAcoustic3D<T>::insertSource(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype){
    Point3D<int> *map;
    size_t ntrace = source->getNtrace();
    size_t nx, ny, nz;
    int nr = 0;
    int lpml;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    size_t i;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            TT[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)] = 0.0;
            nr++;
        }
    }
    return nr;
}

template<typename T>
int RaysAcoustic3D<T>::insertSource(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int traceno){
    Point3D<int> *map;
    size_t ntrace = source->getNtrace();

    if(traceno < 0 || traceno > ntrace-1){
        rs_error("RaysAcoustic3D<T>::insertSource: traceno out of bounds.");
    }

    size_t nx, ny, nz;
    int lpml;
    int nr = 0;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    size_t i = traceno;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
    { 
        TT[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)] = 0.0;
        nr++;
    }
    return nr;
}

template<typename T>
void RaysAcoustic3D<T>::recordData(std::shared_ptr<rockseis::Data3D<T>> data, bool maptype){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *dataarray; 
    T *Fielddata;
    size_t ntrace = data->getNtrace();
    size_t nt = data->getNt();
    size_t nx, ny, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    dataarray = data->getData();
    size_t i;
    Index I(nx, ny, nz);
    Index Idat(nt, ntrace);
    Fielddata = this->getTT();
    for (i=0; i < ntrace; i++) 
    { 
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        {
            dataarray[Idat(0,i)] = Fielddata[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)];
        }
    }

}

template<typename T>
void RaysAcoustic3D<T>::createRecmask(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype){
    Point3D<int> *map;
    size_t ntrace = source->getNtrace();
    size_t nx, ny, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    size_t i;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    Index Idat(1, ntrace); // Data indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            recmask[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)] = true;
        }
    }
}

template<typename T>
void RaysAcoustic3D<T>::clearTT() {
    /* Re-initialize TT array */

    int nx = model->getNx_pml();
    int ny = model->getNy_pml();
    int nz = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny*nz; i++){
        TT[i] = 2.0*TMAX;
    }
}

template<typename T>
void RaysAcoustic3D<T>::clearLam() {
    /* Re-initialize lam array */

    int nx = model->getNx_pml();
    int ny = model->getNy_pml();
    int nz = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny*nz; i++){
        lam[i] = 10.0*TMAX;
    }
}

template<typename T>
void RaysAcoustic3D<T>::insertResiduals(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype){
    Point3D<int> *map;
    int ntrace = source->getNtrace();
    size_t nx, ny, nz;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    T *res;
    res = source->getData();

    size_t i;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    Index Idat(1, ntrace); // Data indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            lam[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)] = res[Idat(0,i)];
        }
    }
}




// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class RaysAcoustic2D<float>;
template class RaysAcoustic2D<double>;
template class RaysAcoustic3D<float>;
template class RaysAcoustic3D<double>;

}

