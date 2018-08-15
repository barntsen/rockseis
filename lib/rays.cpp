// Include statements
#include "rays.h"

#define I2D(i,j) ((i)*nx + (j)) 
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
    tmax = 0.;
    lpml=3;
}

template<typename T>
Rays<T>::~Rays() {
    // Nothing here
}

template<typename T>
Rays<T>::Rays(const int _dim, const int _nx, const int _ny, const int _nz, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const T _tmax) {
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
    tmax = _tmax;
    lpml = 3;

}

// =============== 2D ACOUSTIC RAYS CLASS =============== //
template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(){
    int nx, nz;
    nx = this->getNx();
    nz = this->getNz();

    /* Allocate memory variables */
    TT = (T *) calloc(nx*nz, sizeof(T));
    lam = (T *) calloc(nx*nz, sizeof(T));
    recmask = (bool *) calloc(nx*nz, sizeof(bool));

    /* Initialize arrays */
    for (int i=0; i < nx*nz; i++){
        TT[i] = 2.0*this->getTmax();
        lam[i] = 10.0*this->getTmax();
    }
}

template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(const int _nx, const int _nz, const T _dx, const T _dz, const T _ox, const T _oz, const T _tmax): Rays<T>(2, _nx, 1, _nz, _dx, 1.0, _dz, _ox, 0.0, _oz, _tmax) {

    /* Allocate memory variables */
    TT = (T *) calloc(_nx*_nz, sizeof(T));
    lam = (T *) calloc(_nx*_nz, sizeof(T));
    recmask = (bool *) calloc(_nx*_nz, sizeof(bool));

    /* Initialize TT */
    for (int i=0; i < _nx*_nz; i++){
        TT[i] = 2.0*this->getTmax();
        lam[i] = 10.0*this->getTmax();
    }
}


template<typename T>
RaysAcoustic2D<T>::RaysAcoustic2D(std::shared_ptr<rockseis::ModelAcoustic2D<T>> _model, T _tmax): Rays<T>(){

    int _nx, _ny, _nz;
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
    this->setTmax(_tmax);

    // Setting model pointer
    model = _model;

    /* Allocate memory variables */
    TT = (T *) calloc(_nx*_nz, sizeof(T));
    lam = (T *) calloc(_nx*_nz, sizeof(T));
    recmask = (bool *) calloc(_nx*_nz, sizeof(bool));

    /* Initialize TT */
    for (int i=0; i < _nx*_nz; i++){
        TT[i] = 2.0*_tmax;
        lam[i] = 10.0*_tmax;
    }
}

template<typename T>
RaysAcoustic2D<T>::~RaysAcoustic2D() {
    /* Free allocated variables */
    free(TT);
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

    int nx = model->getNx();
    int ny = model->getNz();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *Vp = model->getVp();

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
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* model borders */
            /* left */
            if((h==0)&&(k>0)&&(k<ny-1)){
                a = fminf(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* right */
            if((h==nx-1)&&(k>0)&&(k<ny-1)){
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k+1,h)]);
            }

            /* top */
            if((k==0)&&(h>0)&&(h<nx-1)){
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }

            /* bottom */
            if((k==ny-1)&&(h>0)&&(h<nx-1)){
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* model corners */
            /* upper-left */
            if((h==0)&&(k==0)){
                a = fminf(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }

            /* lower-left */
            if((h==0)&&(k==ny-1)){
                a = fminf(TT[I2D(k,h)],TT[I2D(k,h+1)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* lower-right */
            if((h==nx-1)&&(k==ny-1)){
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = fminf(TT[I2D(k-1,h)],TT[I2D(k,h)]);
            }

            /* upper-right */
            if((h==nx-1)&&(k==0)){
                a = fminf(TT[I2D(k,h-1)],TT[I2D(k,h)]);
                b = fminf(TT[I2D(k,h)],TT[I2D(k+1,h)]);
            }


            /* calculate solution */
            if(fabs(a-b)>=(S(I2D(k,h))*dh)){
                Tt = fminf(a,b) + S(I2D(k,h))*dh;
            }else{
                Tt = (a + b + sqrt((2.0*pow(S(I2D(k,h)),2.0)*pow(dh,2.0)) - pow((a-b),2.0)))/2.0;
            }

            TT[I2D(k,h)] = fminf(TT[I2D(k,h)],Tt);

            k += ndy;
        }
        h += ndx;
    }
}

template<typename T>
void RaysAcoustic2D<T>::sweep_adj(int nx1, int nx2, int ndx, int ny1, int ny2, int ndy)
{
    int nx = model->getNx();
    int ny = model->getNz();
    T dh = model->getDx();
    T *TT = this->getTT();
    T *lam = this->getLam();
    bool *recmask = this->getRecmask();

    /* local variables */
    int i, j, h, k;
    float app, amp, apm, amm;
    float bpp, bmp, bpm, bmm;
    float ap, am , bp, bm;
    float lhs, rhs, lamt;

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

                app = (ap + fabs(ap))/2.0;
                apm = (ap - fabs(ap))/2.0;

                amp = (am + fabs(am))/2.0;
                amm = (am - fabs(am))/2.0;

                bpp = (bp + fabs(bp))/2.0;
                bpm = (bp - fabs(bp))/2.0;

                bmp = (bm + fabs(bm))/2.0;
                bmm = (bm - fabs(bm))/2.0;


                /* Leung & Qian (2006) */
                lhs = (app-amm)/dh + (bpp-bmm)/dh;
                rhs = (amp*lam[I2D(k,h-1)]-apm*lam[I2D(k,h+1)])/dh + (bmp*lam[I2D(k-1,h)]-bpm*lam[I2D(k+1,h)])/dh;

                /* Taillandier et al. (2009) */
                /* lhs = (apm-amp)/dh + (bpm-bmp)/dh;
                   rhs = (amm*lam[I2D(k,h-1)]-app*lam[I2D(k,h+1)])/dh + (bmm*lam[I2D(k-1,h)]-bpp*lam[I2D(k+1,h)])/dh; */

                lamt = rhs/(lhs+EPS_ADJ);
                lam[I2D(k,h)] = fminf(lam[I2D(k,h)],lamt);
            }

            k += ndy;
        }
        h += ndx;
    }
}


template<typename T>
void RaysAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype){
    Point2D<int> *map;
    int ntrace = source->getNtrace();
    int nx, nz;
    nx = this->getNx();
    nz = this->getNz();

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
            TT[I(map[i].x, map[i].y)] = 0.0;
        }
    }
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
    nx = this->getNx();
    nz = this->getNz();

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
            dataarray[Idat(0,i)] = Fielddata[I(map[i].x, map[i].y)];
        }
    }

}

template<typename T>
T RaysAcoustic2D<T>::norm1(T *TT, T *TTold)
{
    /* local variables */
    int NX = this->getNx();
    int NY = this->getNz();

    int i;
    T sum, norml1;

    /* estimate L1 norm */
    sum = 0.0;
    for (i=0;i<NX*NY;i++){
        sum += (TT[i]-TTold[i]);
    }

    norml1 = sum/(NX*NY);

    return norml1;
}


template<typename T>
void RaysAcoustic2D<T>::solve()
{
    int nx = this->getNx();
    int nz = this->getNz();
    T tmax = this->getTmax();
    T *TTold = (T *) calloc(nx*nz, sizeof(T));
    int i;

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
    }
    free(TTold);
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class RaysAcoustic2D<float>;
template class RaysAcoustic2D<double>;

}


