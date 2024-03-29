// Include statements
#include "rays.h"

#define I2D(i,j) ((j)*nx + (i)) 
#define I3D(i,j,k) ((k)*nx*ny+(j)*nx + (i)) 
#define S(i) (1.0/Vp[i]) 
#define SWAPF(a,b,t) t = b; b = a; a = t;

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

template<typename T>
void Rays<T>::bubble_sort(T *a, T *b, int n) {
   int i = 0, j = 0;
   T tmp;
   for (i = 0; i < n; i++) {   // loop n times - 1 per element
      for (j = 0; j < n - i - 1; j++) { // last i elements are sorted already
         if (a[j] > a[j + 1]) {  // swop if order is broken
            tmp = a[j];
            a[j] = a[j + 1];
            a[j + 1] = tmp;

            tmp = b[j];
            b[j] = b[j + 1];
            b[j + 1] = tmp;
         }
      }
   }
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
    adjsource = (T *) calloc(nx*nz, sizeof(T));

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
    adjsource = (T *) calloc(nx*nz, sizeof(T));

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
    adjsource = (T *) calloc(nx_pml*nz_pml, sizeof(T));

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
    free(adjsource);
}

template<typename T>
void RaysAcoustic2D<T>::sweep(int nx1, int nx2, int ndx, int nz1, int nz2, int ndz)
{ 
/*------------------------------------------------------------------------
 *  Solve eikonal equation by fast sweeping method according to Zhao (2004)
 *
 *  ----------------------------------------------------------------------*/

    int nx = model->getNx_pml();
    int nz = model->getNz_pml();
    T dx = model->getDx();
    T dz = model->getDz();
    T *TT = this->getTT();
    T *Vp = model->getL();

    /* local variables */
    int i, j, h, k;
    T a, b, Tta, Ttb, Tt;

    /* sweep over FD-grid */
    h = nz1;
    for (j=0;j<nz;j++){	
       k = nx1;
       for (i=0;i<nx;i++){
          b = 0 == k ? TT[I2D(1,h)] : ((nx - 1) == k) ? TT[I2D(nx-2,h)] : MIN(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
          a = 0 == h ? TT[I2D(k,1)] : ((nz - 1) == h) ? TT[I2D(k,nz-2)] : MIN(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);

          /* calculate solution */
          Tta = a + S(I2D(k,h))*dz;
          Ttb = b + S(I2D(k,h))*dx;
          if(Tta <= b){
             Tt = Tta;
          }else if (Ttb <= a){
             Tt = Ttb;
          }else{
             Tt = (a/SQ(dz) + b/SQ(dx) + sqrt((SQ(dx)+SQ(dz))*SQ(S(I2D(k,h))) - SQ(a-b))/(dx*dz))/(SQ(1./dx)+SQ(1./dz));
          }
          TT[I2D(k,h)] = MIN(TT[I2D(k,h)],Tt);

          k += ndx;
       }
       h += ndz;
    }
}

template<typename T>
void RaysAcoustic2D<T>::sweep_adj(int nx1, int nx2, int ndx, int nz1, int nz2, int ndz)
{
    int nx = model->getNx_pml();
    int nz = model->getNz_pml();
    T dx = model->getDx();
    T dz = model->getDz();
    T *TT = this->getTT();
    T *lam = this->getLam();
    T *adjsource = this->getAdjsource();

    /* local variables */
    int i, j, h, k;
    T app, amp, apm, amm;
    T bpp, bmp, bpm, bmm;
    T ap, am , bp, bm;
    T lhs, rhs, lamt;

    /* sweep over FD-grid */
    k = nx1;
    for (i=1;i<nx-1;i++){
       h = nz1;
       for (j=1;j<nz-1;j++){
          /* assemble equation (3.6) in Leung & Qian (2006) */
          ap = -(TT[I2D(k,h+1)]-TT[I2D(k,h)])/dz;
          am = -(TT[I2D(k,h)]-TT[I2D(k,h-1)])/dz;

          bp = -(TT[I2D(k+1,h)]-TT[I2D(k,h)])/dx;
          bm = -(TT[I2D(k,h)]-TT[I2D(k-1,h)])/dx;

          app = (ap + ABS(ap))/2.0;
          apm = (ap - ABS(ap))/2.0;

          amp = (am + ABS(am))/2.0;
          amm = (am - ABS(am))/2.0;

          bpp = (bp + ABS(bp))/2.0;
          bpm = (bp - ABS(bp))/2.0;

          bmp = (bm + ABS(bm))/2.0;
          bmm = (bm - ABS(bm))/2.0;


          /* Leung & Qian (2006) */
          lhs = (app-amm)/dz + (bpp-bmm)/dx;
          rhs = (amp*lam[I2D(k,h-1)]-apm*lam[I2D(k,h+1)])/dz + (bmp*lam[I2D(k-1,h)]-bpm*lam[I2D(k+1,h)])/dx - adjsource[I2D(k,h)];
          //
          lamt = rhs*lhs/(SQ(lhs)+EPS_ADJ);
          lam[I2D(k,h)] = MIN(lam[I2D(k,h)],lamt);
          h += ndz;
       }
       k += ndx;
    }
}
template<typename T>
void RaysAcoustic2D<T>::copyTT(T *TTin) {
    /* Copy traveltimes from other array */

    int nx = model->getNx_pml();
    int nz = model->getNz_pml();
    int lpml = model->getLpml();
    Index Irays(nx, nz);

    /* Copy TT */
    for(size_t ix=lpml; ix < nx-lpml; ix++){
        for(size_t iz=lpml; iz < nz-lpml; iz++){
            TT[Irays(ix, iz)] = TTin[Irays(ix, iz)]; 
        }
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
void RaysAcoustic2D<T>::clearLam(T val) {
    /* Re-initialize lam array */

    int nx = model->getNx_pml();
    int ny = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny; i++){
        lam[i] = val;
    }
}

//template<typename T>
//int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
//    Point2D<int> *map;
//    int ntrace = source->getNtrace();
//    int nx, nz;
//    int lpml;
//    int nr = 0;
//    lpml = this->getLpml();
//    nx = this->getNx_pml();
//    nz = this->getNz_pml();
//
//    // Get correct map (source or receiver mapping)
//    if(maptype == SMAP) {
//        map = (source->getGeom())->getSmap();
//    }else{
//        map = (source->getGeom())->getGmap();
//    }
//
//    int i;
//    //Indexes 
//    Index I(nx, nz); //Model and Field indexes
//    for (i=0; i < ntrace; i++) 
//    {
//        if(map[i].x >= 0 && map[i].y >=0)
//        { 
//            TT[I(lpml+map[i].x, lpml+map[i].y)] = 0.0;
//            nr++;
//        }
//    }
//    return nr;
//}
//
//template<typename T>
//int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int traceno){
//    Point2D<int> *map;
//    int ntrace = source->getNtrace();
//    if(traceno < 0 || traceno > ntrace-1){
//        rs_error("RaysAcoustic2D<T>::insertSource: traceno out of bounds.");
//    }
//    int nx, nz;
//    int lpml;
//    lpml = this->getLpml();
//    nx = this->getNx_pml();
//    nz = this->getNz_pml();
//    int nr = 0;
//
//    // Get correct map (source or receiver mapping)
//    if(maptype == SMAP) {
//        map = (source->getGeom())->getSmap();
//    }else{
//        map = (source->getGeom())->getGmap();
//    }
//
//    int i = traceno;
//    //Indexes 
//    Index I(nx, nz); //Model and Field indexes
//    if(map[i].x >= 0 && map[i].y >=0)
//    { 
//        TT[I(lpml+map[i].x, lpml+map[i].y)] = 0.0;
//        nr++;
//    }
//    return nr;
//}

template<typename T>
int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
    Point2D<int> *map;
    Point2D<T> *shift;
    int ntrace = source->getNtrace();
    int nx, nz;
    T dx, dz;
    int lpml;
    int nr = 0;
    int isx, isz;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();
    dx = this->getDx();
    dz = this->getDz();
    T x,z;
    T *V = model->getL();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    int i;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0)
        { 
            for(isx =-1; isx < 2 ; isx++)
            {
                if(((map[i].x + isx) >= 0) && ((map[i].x + isx) < nx)){
                    x = shift[i].x -  isx;
                    for(isz = -1; isz < 2 ; isz++){
                        z = shift[i].y -  isz;
                        if(((map[i].y + isz) >= 0) && ((map[i].y + isz) < nz)){
                            TT[I(lpml+map[i].x+isx, lpml+map[i].y+isz)] = sqrt(SQ(x*dx) + SQ(z*dz))/V[I(lpml+map[i].x, lpml+map[i].y)];
                        }
                    }
                }
            }
            nr++;
        }
    }
    return nr;
}

template<typename T>
int RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int traceno){
    Point2D<int> *map;
    Point2D<T> *shift;
    int nx, nz;
    T dx, dz;
    int lpml;
    int nr = 0;
    int isx, isz;
    lpml = this->getLpml();
    nx = this->getNx_pml();
    nz = this->getNz_pml();
    dx = this->getDx();
    dz = this->getDz();
    T x,z;
    T *V = model->getL();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    int i = traceno;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
        if(map[i].x >= 0 && map[i].y >=0)
        { 
            for(isx = -1; isx < 2 ; isx++)
            {
                if(((map[i].x + isx) >= 0) && ((map[i].x + isx) < nx)){
                    x = shift[i].x -  isx;
                    for(isz = -1; isz < 2 ; isz++){
                        z = shift[i].y -  isz;
                        if(((map[i].y + isz) >= 0) && ((map[i].y + isz) < nz)){
                            TT[I(lpml+map[i].x+isx, lpml+map[i].y+isz)] = sqrt(SQ(x*dx) + SQ(z*dz))/V[I(lpml+map[i].x, lpml+map[i].y)];
                        }
                    }
                }
            }
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
    int i, i1, i2;
    Index I(nx, nz);
    Index Idat(nt, ntrace);
    Fielddata = this->getTT();
    for (i=0; i < ntrace; i++) 
    { 
        if(map[i].x >= 0 && map[i].y >=0)
        {
            for(i1=0; i1<2*LANC_SIZE; i1++){
                for(i2=0; i2<2*LANC_SIZE; i2++){
                    dataarray[Idat(0,i)] += Fielddata[I(lpml + map[i].x - (LANC_SIZE) + i2, lpml + map[i].y - (LANC_SIZE) + i1)]*LANC(shift[i].x + LANC_SIZE - i2, LANC_SIZE)*LANC(shift[i].y + LANC_SIZE - i1 ,LANC_SIZE);
                }
            }
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
            adjsource[I(lpml+map[i].x, lpml+map[i].y)] =  -res[Idat(0,i)];
        }
    }
}

template<typename T>
void RaysAcoustic2D<T>::insertImageresiduals(T *res){
    int nx, nz;
    int nx_pml, nz_pml;
    int lpml;
    lpml = this->getLpml();
    nx = this->getNx();
    nz = this->getNz();
    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();

    int ix,iz;
    //Indexes 
    Index Iray(nx_pml, nz_pml); //Model and Field indexes
    Index Ires(nx,nz); // Data indexes
    for (iz=0; iz < nz; iz++) {
        for (ix=0; ix < nx; ix++) {
            adjsource[Iray(lpml+ix, lpml+iz)] =  +res[Ires(ix,iz)];
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
    adjsource = (T *) calloc(nx*ny*nz, sizeof(T));

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
    adjsource = (T *) calloc(nx*ny*nz, sizeof(T));

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
    adjsource = (T *) calloc(nx*ny*nz, sizeof(T));

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
    free(adjsource);
}

template<typename T>
void RaysAcoustic3D<T>::sweep(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy, int nz1, int nz2, int ndz)
{ 
    size_t nx = model->getNx_pml();
    size_t ny = model->getNy_pml();
    size_t nz = model->getNz_pml();
    T dx = model->getDx();
    T dy = model->getDy();
    T dz = model->getDz();
    T *TT = this->getTT();
    T *Vp = model->getL();

    /* local variables */
    size_t i, j, k, m, n, l;
    T a[3], d[3];
    T Tt;

    T slo;
    T swap;

    /* sweep over FD-grid */
    m = nx1;
    for (i=0;i<nx;i++){
        n = ny1;
        for (j=0;j<ny;j++){	
            l = nz1;
            for (k=0;k<nz;k++){	
                a[0] = (0 == m) ? TT[I3D(1,n,l)] : ((nx - 1) == m) ? TT[I3D(nx-2,n,l)] :
                    MIN(TT[I3D(m-1,n,l)],TT[I3D(m+1,n,l)]);

                a[1] = (0 == n) ? TT[I3D(m,1,l)] : ((ny - 1) == n) ? TT[I3D(m,ny-2,l)] :
                    MIN(TT[I3D(m,n-1,l)],TT[I3D(m,n+1,l)]);

                a[2] = (0 == l) ? TT[I3D(m,n,1)] : ((nz - 1) == l) ? TT[I3D(m,n,nz-2)] :
                    MIN(TT[I3D(m,n,l-1)],TT[I3D(m,n,l+1)]);


                /* calculate solution */
                slo = S(I3D(m,n,l));
                d[0] = dx;
                d[1] = dy;
                d[2] = dz;
                if (a[0] > a[2]) {
                    SWAPF (a[0], a[2], swap)
                    SWAPF (d[0], d[2], swap)
                }
                if (a[0] > a[1]) {
                    SWAPF (a[0], a[1], swap)
                    SWAPF (d[0], d[1], swap)
                }
                if (a[1] > a[2]) {
                    SWAPF (a[1], a[2], swap)
                    SWAPF (d[1], d[2], swap)
                }

                Tt = a[0] + d[0]*slo;
                if(Tt > a[1]){
                    Tt = (d[0] * d[1] * sqrt ((SQ(d[0]) + SQ(d[1])) * slo * slo - (a[0] - a[1]) * (a[0] - a[1])) +
                            a[0] * SQ(d[1]) + a[1] * SQ(d[0])) / (SQ(d[0]) + SQ(d[1]));
                    if(Tt > a[2]){
                        Tt = (d[0] * d[1] * d[2] * sqrt ((SQ(d[1]) * SQ(d[2]) + SQ(d[0]) * SQ(d[2]) + SQ(d[0]) * SQ(d[1])) * slo * slo
                                    + 2.0 * SQ(d[2]) * a[0] * a[1]
                                    + 2.0 * SQ(d[1]) * a[0] * a[2]
                                    + 2.0 * SQ(d[0]) * a[1] * a[2]
                                    - (SQ(d[2]) + SQ(d[1])) * a[0] * a[0]
                                    - (SQ(d[2]) + SQ(d[0])) * a[1] * a[1]
                                    - (SQ(d[1]) + SQ(d[0])) * a[2] * a[2])
                                + SQ(d[0]) * SQ(d[1]) * a[2] + SQ(d[0]) * SQ(d[2]) * a[1] + SQ(d[1]) * SQ(d[2]) * a[0])
                            / (SQ(d[1]) * SQ(d[2]) + SQ(d[0]) * SQ(d[2]) + SQ(d[0]) * SQ(d[1]));
                    }
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
    T dx = model->getDx();
    T dy = model->getDy();
    T dz = model->getDz();
    T *TT = this->getTT();
    T *lam = this->getLam();
    T *adjsource = this->getAdjsource();

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
                    /* assemble equation (3.6) in Leung & Qian (2006) */
                    ap = -(TT[I3D(m+1,n,l)]-TT[I3D(m,n,l)])/dx;
                    am = -(TT[I3D(m,n,l)]-TT[I3D(m-1,n,l)])/dx;

                    bp = -(TT[I3D(m,n+1,l)]-TT[I3D(m,n,l)])/dy;
                    bm = -(TT[I3D(m,n,l)]-TT[I3D(m,n-1,l)])/dy;

                    cp = -(TT[I3D(m,n,l+1)]-TT[I3D(m,n,l)])/dz;
                    cm = -(TT[I3D(m,n,l)]-TT[I3D(m,n,l-1)])/dz;

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
                    lhs = (app-amm)/dx + (bpp-bmm)/dy + (cpp-cmm)/dz;
                    rhs = (amp*lam[I3D(m-1,n,l)]-apm*lam[I3D(m+1,n,l)])/dx + (bmp*lam[I3D(m,n-1,l)]-bpm*lam[I3D(m,n+1,l)])/dy + (cmp*lam[I3D(m,n,l-1)]-cpm*lam[I3D(m,n,l+1)])/dz - adjsource[I3D(m,n,l)];

                    lamt = rhs*lhs/(SQ(lhs)+EPS_ADJ);
                    lam[I3D(m,n,l)] = MIN(lam[I3D(m,n,l)],lamt);
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
    while((lnorm1>=TTNORM)){

        /* save old TT values */
        for (i=0; i<nx*ny*nz; i++){
            TTold[i] = TT[i];
        }

        sweep(0, nx-1, 1,  0, ny-1, 1,  0, nz-1, 1);
        sweep(0, nx-1, 1,  0, ny-1, 1,  nz-1, 0, -1);
        sweep(0, nx-1, 1,  ny-1, 0, -1,  0, nz-1, 1);
        sweep(0, nx-1, 1,  ny-1, 0, -1,  nz-1, 0, -1);
        sweep(nx-1, 0, -1,  0, ny-1, 1,  0, nz-1, 1);
        sweep(nx-1, 0, -1,  0, ny-1, 1,  nz-1, 0, -1);
        sweep(nx-1, 0, -1,  ny-1, 0, -1,  0, nz-1, 1);
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
    Point3D<T> *shift;
    size_t ntrace = source->getNtrace();
    size_t nx, ny, nz;
    T dx, dy, dz;
    int nr = 0;
    int lpml;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    dx = this->getDx();
    dy = this->getDy();
    dz = this->getDz();

    T x,y,z;
    T *V = model->getL();


    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    size_t i;
    int isx, isy, isz;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            for(isx = -1; isx < 2 ; isx++)
            {
                if(((map[i].x + isx) >= 0) && ((map[i].x + isx) < nx)){
                    x = shift[i].x -  isx;
                    for(isy = -1; isy < 2 ; isy++){
                        if(((map[i].y + isy) >= 0) && ((map[i].y + isy) < ny)){
                            y = shift[i].y -  isy;
                            for(isz = -1; isz < 2 ; isz++){
                                if(((map[i].z + isz) >= 0) && ((map[i].z + isz) < nz)){
                                    z = shift[i].z -  isz;
                                    TT[I(lpml+map[i].x + isx, lpml+map[i].y + isy, lpml+map[i].z + isz)] = sqrt(SQ(x*dx) + SQ(y*dy) + SQ(z*dz))/V[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)];
                                    nr++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return nr;
}

template<typename T>
int RaysAcoustic3D<T>::insertSource(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, size_t traceno){
    Point3D<int> *map;
    Point3D<T> *shift;
    size_t nx, ny, nz;
    T dx, dy, dz;
    int nr = 0;
    int lpml;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    dx = this->getDx();
    dy = this->getDy();
    dz = this->getDz();

    T x,y,z;
    T *V = model->getL();


    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    size_t i = traceno;
    int isx, isy, isz;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            for(isx = -1; isx < 2 ; isx++)
            {
                if(((map[i].x + isx) >= 0) && ((map[i].x + isx) < nx)){
                    x = shift[i].x -  isx;
                    for(isy = -1; isy < 2 ; isy++){
                        if(((map[i].y + isy) >= 0) && ((map[i].y + isy) < ny)){
                            y = shift[i].y -  isy;
                            for(isz = -1; isz < 2 ; isz++){
                                if(((map[i].z + isz) >= 0) && ((map[i].z + isz) < nz)){
                                    z = shift[i].z -  isz;
                                    TT[I(lpml+map[i].x + isx, lpml+map[i].y + isy, lpml+map[i].z + isz)] = sqrt(SQ(x*dx) + SQ(y*dy) + SQ(z*dz))/V[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)];
                                }
                            }
                        }
                    }
                }
            }
            nr++;
        }
    return nr;
}

template<typename T>
int RaysAcoustic3D<T>::solveHomogen(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype){
    Point3D<int> *map;
    Point3D<T> *shift;
    size_t ntrace = source->getNtrace();
    size_t nx, ny, nz;
    T dx, dy, dz;
    int nr = 0;
    int lpml;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    dx = this->getDx();
    dy = this->getDy();
    dz = this->getDz();

    T x,y,z;
    T *V = model->getL();


    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    size_t i;
    int isx, isy, isz;
    int nxi,nyi,nzi;
    nxi = (int) nx;
    nyi = (int) ny;
    nzi = (int) nz;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    for (i=0; i < ntrace; i++) 
    {
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        { 
            for(isx = -map[i].x; isx < (nxi-map[i].x) ; isx++)
                {
                if(((map[i].x + isx + lpml) >= 0) && ((map[i].x + isx + lpml) < nx)){
                    x = shift[i].x -  isx;
                    for(isy = -map[i].y; isy < (nyi-map[i].y) ; isy++){
                        if(((map[i].y + isy + lpml) >= 0) && ((map[i].y + isy + lpml) < ny)){
                            y = shift[i].y -  isy;
                            for(isz = -map[i].z; isz < (nzi-map[i].z) ; isz++){
                                if(((map[i].z + isz + lpml) >= 0) && ((map[i].z + isz + lpml) < nz)){
                                    z = shift[i].z -  isz;
                                    TT[I(lpml+map[i].x + isx, lpml+map[i].y + isy, lpml+map[i].z + isz)] = sqrt(SQ(x*dx) + SQ(y*dy) + SQ(z*dz))/V[I(lpml, lpml, lpml)];
                                    nr++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return nr;
}

template<typename T>
int RaysAcoustic3D<T>::solveHomogen(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, size_t traceno){
    Point3D<int> *map;
    Point3D<T> *shift;
    size_t nx, ny, nz;
    T dx, dy, dz;
    int nr = 0;
    int lpml;

    lpml = this->getLpml();
    nx = this->getNx_pml();
    ny = this->getNy_pml();
    nz = this->getNz_pml();

    dx = this->getDx();
    dy = this->getDy();
    dz = this->getDz();

    T x,y,z;
    T *V = model->getL();


    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    size_t i = traceno;
    int isx, isy, isz;
    int nxi,nyi,nzi;
    nxi = (int) nx;
    nyi = (int) ny;
    nzi = (int) nz;
    //Indexes 
    Index I(nx, ny, nz); //Model and Field indexes
    if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
    { 
        for(isx = -map[i].x; isx < (nxi-map[i].x) ; isx++)
        {
            if(((map[i].x + isx + lpml) >= 0) && ((map[i].x + isx + lpml) < nx)){
                x = shift[i].x -  isx;
                for(isy = -map[i].y; isy < (nyi-map[i].y) ; isy++){
                    if(((map[i].y + isy + lpml) >= 0) && ((map[i].y + isy + lpml) < ny)){
                        y = shift[i].y -  isy;
                        for(isz = -map[i].z; isz < (nzi-map[i].z) ; isz++){
                            if(((map[i].z + isz + lpml) >= 0) && ((map[i].z + isz + lpml) < nz)){
                                z = shift[i].z -  isz;
                                TT[I(lpml+map[i].x + isx, lpml+map[i].y + isy, lpml+map[i].z + isz)] = sqrt(SQ(x*dx) + SQ(y*dy) + SQ(z*dz))/V[I(lpml, lpml, lpml)];
                                nr++;
                            }
                        }
                    }
                }
            }
        }
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
    size_t i,i1,i2,i3;
    Index I(nx, ny, nz);
    Index Idat(nt, ntrace);
    Fielddata = this->getTT();
    for (i=0; i < ntrace; i++) 
    { 
        if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
        {
            for(i1=0; i1<2*LANC_SIZE; i1++){
                for(i2=0; i2<2*LANC_SIZE; i2++){
                    for(i3=0; i3<2*LANC_SIZE; i3++){
                        dataarray[Idat(0,i)] += Fielddata[I(lpml + map[i].x - (LANC_SIZE) + i3, lpml + map[i].y - (LANC_SIZE) + i2, lpml + map[i].z - (LANC_SIZE) + i1)]*LANC(shift[i].x + LANC_SIZE - i3, LANC_SIZE)*LANC(shift[i].y + LANC_SIZE - i2 ,LANC_SIZE)*LANC(shift[i].z + LANC_SIZE - i1 ,LANC_SIZE);
                    }
                }
            }
        }
    }

}

template<typename T>
void RaysAcoustic3D<T>::copyTT(T *TTin) {
    /* Copy traveltimes from other array */

    int nx = model->getNx_pml();
    int ny = model->getNy_pml();
    int nz = model->getNz_pml();
    int lpml = model->getLpml();
    Index Irays(nx, ny, nz);

    /* Copy TT */
    for(size_t ix=lpml; ix < nx-lpml; ix++){
        for(size_t iy=lpml; iy < ny-lpml; iy++){
            for(size_t iz=lpml; iz < nz-lpml; iz++){
                TT[Irays(ix, iy, iz)] = TTin[Irays(ix, iy, iz)]; 
            }
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
void RaysAcoustic3D<T>::clearLam(T val) {
    /* Re-initialize lam array */

    int nx = model->getNx_pml();
    int ny = model->getNy_pml();
    int nz = model->getNz_pml();

    /* Initialize TT */
    for (int i=0; i < nx*ny*nz; i++){
        lam[i] = val;
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
            adjsource[I(lpml+map[i].x, lpml+map[i].y, lpml+map[i].z)] = -res[Idat(0,i)];
        }
    }
}

template<typename T>
void RaysAcoustic3D<T>::insertImageresiduals(T *res){
    int nx, ny, nz;
    int nx_pml, ny_pml, nz_pml;
    int lpml;
    lpml = this->getLpml();

    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    nx_pml = this->getNx_pml();
    ny_pml = this->getNy_pml();
    nz_pml = this->getNz_pml();

    int ix,iy,iz;
    //Indexes 
    Index Iray(nx_pml, ny_pml, nz_pml); //Model and Field indexes
    Index Ires(nx, ny, nz); // Data indexes
    for (iz=0; iz < nz; iz++) {
        for (iy=0; iy < ny; iy++) {
            for (ix=0; ix < nx; ix++) {
                adjsource[Iray(lpml+ix, lpml+iy, lpml+iz)] =  -res[Ires(ix,iy,iz)];
            }
        }
    }
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class RaysAcoustic2D<float>;
template class RaysAcoustic2D<double>;
template class RaysAcoustic3D<float>;
template class RaysAcoustic3D<double>;

}

