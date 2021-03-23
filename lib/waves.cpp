// Include statements
#include "waves.h"

#define I2D(i,j) ((j)*nx + (i)) 
#define I2D_lr(i,j) ((j)*lpml + (i)) 
#define I2D_tb(i,j) ((j)*nx + (i)) 

#define I3D(i,j,k) ((k)*nx*ny + (j)*nx + (i)) 
#define I3D_lr(i,j,k) ((k)*lpml*ny + (j)*lpml + (i)) 
#define I3D_tb(i,j,k) ((k)*nx*ny + (j)*nx + (i)) 
#define I3D_fb(i,j,k) ((k)*nx*lpml + (j)*nx + (i)) 

#define Idat(i,j) ((j)*nt + (i))

namespace rockseis {

// =============== ABSTRACT WAVES CLASS =============== //
template<typename T>
Waves<T>::Waves() {
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, 1);
    geometry->setN(2, 1);
    geometry->setN(3, 1);
    geometry->setN(4, 1);
    geometry->setD(1, 1.);
    geometry->setD(2, 1.);
    geometry->setD(3, 1.);
    geometry->setD(4, 1.);
    geometry->setO(1, 0.);
    geometry->setO(2, 0.);
    geometry->setO(3, 0.);
    geometry->setO(4, 0.);
    dim=0;
    lpml = 10;
    domdec = false;
}

template<typename T>
Waves<T>::~Waves() {
    // Nothing here
}

template<typename T>
Waves<T>::Waves(const int _dim, const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt,  const T _ox, const T _oy, const T _oz, const T _ot) {
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, _nx);
    geometry->setN(2, _ny);
    geometry->setN(3, _nz);
    geometry->setN(4, _nt);
    geometry->setD(1, _dx);
    geometry->setD(2, _dy);
    geometry->setD(3, _dz);
    geometry->setD(4, _dt);
    geometry->setO(1, _ox);
    geometry->setO(2, _oy);
    geometry->setO(3, _oz);
    geometry->setO(4, _ot);
    dim = _dim;
    lpml = _lpml;
    domdec = false;
}


template<typename T>
int Waves<T>::getNx_pml() {
    if(this->getDomdec()){
        return geometry->getN(1);
    }else{
        return geometry->getN(1) + 2*lpml;
    }
}

template<typename T>
int Waves<T>::getNy_pml() {
    if(this->getDomdec()){
        return geometry->getN(2);
    }else{
        return geometry->getN(2) + 2*lpml;
    }
}

template<typename T>
int Waves<T>::getNz_pml() {
    if(this->getDomdec()){
        return geometry->getN(3);
    }else{
        return geometry->getN(3) + 2*lpml;
    }
}

// =============== 1D ACOUSTIC MODEL CLASS =============== //
template<typename T>
WavesAcoustic1D<T>::WavesAcoustic1D(){
    int nz, lpml, nz_pml;
    T dt;
    nz = this->getNz();
    lpml = this->getLpml();
    dt = this->getDt();
    this->setDim(1);
    nz_pml = nz + 2*lpml;
    Pml = std::make_shared<PmlAcoustic1D<T>>(nz, lpml, dt);
    P1 = (T *) calloc(nz_pml,sizeof(T));
    P2 = (T *) calloc(nz_pml,sizeof(T));
    Az = (T *) calloc(nz_pml,sizeof(T));
}

template<typename T>
WavesAcoustic1D<T>::WavesAcoustic1D(const int _nz, const int _nt, const int _lpml, const T _dz, const T _dt, const T _oz, const T _ot): Waves<T>(1, 1, 1, _nz, _nt, _lpml, 1.0, 1.0, _dz, _dt, 0.0, 0.0, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlAcoustic1D<T>>(_nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nz_pml;
    this->setDim(1);
    nz_pml = _nz + 2*_lpml;
    P1 = (T *) calloc(nz_pml,sizeof(T));
    P2 = (T *) calloc(nz_pml,sizeof(T));
    Az = (T *) calloc(nz_pml,sizeof(T));
}


template<typename T>
WavesAcoustic1D<T>::WavesAcoustic1D(std::shared_ptr<rockseis::ModelAcoustic1D<T>> model, int _nt, T _dt, T _ot): Waves<T>(){

    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    /* Create associated PML class */
    Pml = std::make_shared<PmlAcoustic1D<T>>(_nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nz_pml;
    this->setDim(1);
    nz_pml = _nz + 2*_lpml;
    P1 = (T *) calloc(nz_pml,sizeof(T));
    P2 = (T *) calloc(nz_pml,sizeof(T));
    Az = (T *) calloc(nz_pml,sizeof(T));
}

template<typename T>
WavesAcoustic1D<T>::~WavesAcoustic1D() {
    /* Free allocated variables */
    free(P1);
    free(P2);
    free(Az);
}

template<typename T>
void WavesAcoustic1D<T>::forwardstepAcceleration(std::shared_ptr<rockseis::ModelAcoustic1D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, iz, nz, lpml;
    lpml = model->getLpml();
    nz = model->getNz() + 2*lpml;
    T *Rz, *df;
    Rz = model->getRz();
    df = der->getDf();
    
    // Derivate P forward with respect to z
    der->ddz_fw(P1);
    // Compute Az
    for(iz=0; iz < nz; iz++){
            Az[iz] = Rz[iz]*df[iz];
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
            // Top
            Pml->P_top[iz] = Pml->B_ltf_stag[iz]*Pml->P_top[iz] + Pml->A_ltf_stag[iz]*df[iz];
            
            Az[iz] -= Rz[iz]*(Pml->P_top[iz] + Pml->C_ltf_stag[iz]*df[iz]);
            i = iz + nz - lpml;
            //Bottom
            Pml->P_bottom[iz] = Pml->B_rbb_stag[iz]*Pml->P_bottom[iz] + Pml->A_rbb_stag[iz]*df[i];
            Az[i] -= Rz[i]*(Pml->P_bottom[iz] + Pml->C_rbb_stag[iz]*df[i]);
        }
    
} // End of forwardstepAcceleration

template<typename T>
void WavesAcoustic1D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelAcoustic1D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, iz, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *L, *df;
    L = model->getL();
    df = der->getDf();
    
    // Derivate Az backward with respect to z
    der->ddz_bw(Az);
    // Compute P
    for(iz=0; iz < nz; iz++){
            P2[iz] = 2.0 * P1[iz] - P2[iz] +  dt*dt*L[iz]*df[iz];
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
            // Top
            Pml->Azz_top[iz] = Pml->B_ltf[iz]*Pml->Azz_top[iz] + Pml->A_ltf[iz]*df[iz];
            
            P2[iz] -= dt*dt*L[iz]*(Pml->Azz_top[iz] + Pml->C_ltf[iz]*df[iz]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Azz_bottom[iz] = Pml->B_rbb[iz]*Pml->Azz_bottom[iz] + Pml->A_rbb[iz]*df[i];
            P2[i] -= dt*dt*L[i]*(Pml->Azz_bottom[iz] + Pml->C_rbb[iz]*df[i]);
    }
    
}

template<typename T>
void WavesAcoustic1D<T>::insertSource(std::shared_ptr<rockseis::ModelAcoustic1D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    T *wav; 
    int ntrace = source->getNtrace();
    int lpml;
    lpml = this->getLpml();
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i;
    //Indexes 
    switch(sourcetype)
    {
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].y >=0)
                { 
                    Az[lpml + map[i].y] += Mod[lpml + map[i].y]*wav[Idat(it,i)]; 
                }
            }
            break;
        case PRESSURE:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].y >=0)
                { 
                    P2[lpml + map[i].y] += dt*dt*Mod[lpml + map[i].y]*wav[Idat(it,i)]; 
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesAcoustic1D<T>::recordData(std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it){
    Point2D<int> *map;
    T *dataarray; 
    T *Fielddata;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int lpml;
    lpml = this->getLpml();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
    }else{
        map = (data->getGeom())->getGmap();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i;
    switch(field)
    {
        case PRESSURE:
            Fielddata = this->getP1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = Fielddata[lpml + map[i].y];
                }
            }
	    break;
        case VZ:
            Fielddata = this->getAz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = Fielddata[lpml + map[i].y];
                }
            }
            break;
        default:
            break;
    }

}


// =============== 2D ACOUSTIC WAVES CLASS =============== //
template<typename T>
WavesAcoustic2D<T>::WavesAcoustic2D(){
    int nx, nz, lpml, nx_pml, nz_pml;
    T dt;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    dt = this->getDt();
    this->setDim(2);
    nx_pml = nx + 2*lpml;
    nz_pml = nz + 2*lpml;
    Pml = std::make_shared<PmlAcoustic2D<T>>(nx, nz, lpml, dt);
    P1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    P2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ax = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Az = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesAcoustic2D<T>::WavesAcoustic2D(const int _nx, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot): Waves<T>(2, _nx, 1, _nz, _nt, _lpml, _dx, 1.0, _dz, _dt, _ox, 0.0, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlAcoustic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    P1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    P2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ax = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Az = (T *) calloc(nx_pml*nz_pml,sizeof(T));

}


template<typename T>
WavesAcoustic2D<T>::WavesAcoustic2D(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, int _nt, T _dt, T _ot): Waves<T>(){

    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    if((model->getDomain())->getStatus()){
       // Domain decomposition mode
       this->setDomain(model->getDomain());
       int ix0, nxo, iz0, nzo;
       bool low[2],high[2];

       for (int i=0; i<2; i++){
          low[i] = false;
          high[i] = false;
        }
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
       if(ix0 < _lpml) low[0]=true;
       if(ix0 + _nx >= nxo-_lpml) high[0]=true;
       if(iz0 < _lpml) low[1]=true;
       if(iz0 + _nz >= nzo-_lpml) high[1]=true;

       /* Create associated PML class */
       Pml = std::make_shared<PmlAcoustic2D<T>>(_nx, _nz, _lpml, _dt, &low[0], &high[0]);

       int nx_pml, nz_pml;
       nx_pml = _nx;
       nz_pml = _nz;

       P1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       P2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Ax = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Az = (T *) calloc(nx_pml*nz_pml,sizeof(T));

    }else{
       /* Create associated PML class */
       Pml = std::make_shared<PmlAcoustic2D<T>>(_nx, _nz, _lpml, _dt);

       /* Allocate memory variables */
       int nx_pml, nz_pml;
       this->setDim(2);
       nx_pml = _nx + 2*_lpml;
       nz_pml = _nz + 2*_lpml;
       P1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       P2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Ax = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Az = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    }

}

template<typename T>
WavesAcoustic2D<T>::~WavesAcoustic2D() {
    /* Free allocated variables */
    free(P1);
    free(P2);
    free(Ax);
    free(Az);
}

template<typename T>
void WavesAcoustic2D<T>::forwardstepAcceleration(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iz0, iz, nx, nz, lpml, nxo, nzo;
    lpml = model->getLpml();

    T *Rx, *Rz, *df;
    Rx = model->getRx();
    Rz = model->getRz();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
       // Domain decomposition 
       nx = model->getNx();
       nz = model->getNz();
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
    }else{
       nx = model->getNx() + 2*lpml;
       nz = model->getNz() + 2*lpml;
       ix0 = 0;
       iz0 = 0;
       nxo = nx;
       nzo = nz;
    }
    
    // Derivate P forward with respect to x
    der->ddx_fw(P1);
    // Compute Ax
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Ax[I2D(ix,iz)] = Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(ix=0; ix < lpml; ix++){
             if(Pml->getApplypml(0)){
                 if(ix >= ix0 && ix < (ix0 + nx)){
                     // Left
                     Pml->P_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->P_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix-ix0,iz)];
                     Ax[I2D(ix-ix0,iz)] -= Rx[I2D(ix-ix0,iz)]*(Pml->P_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix-ix0,iz)]);
                 }
             }
             if(Pml->getApplypml(1)){
                 i = ix + nxo - lpml;
                 if(i >= ix0 && i < (ix0 + nx)){
                     // Right
                     Pml->P_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->P_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i-ix0,iz)];
                     Ax[I2D(i-ix0,iz)] -= Rx[I2D(i-ix0,iz)]*(Pml->P_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i-ix0,iz)]);
                 }
             }
          }
       }
    }

    
    // Derivate P forward with respect to z
    der->ddz_fw(P1);
    // Compute Az
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Az[I2D(ix,iz)] = Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(ix=0; ix < nx; ix++){
              if(Pml->getApplypml(4)){
                  if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->P_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->P_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz-iz0)];
                      Az[I2D(ix,iz-iz0)] -= Rz[I2D(ix,iz-iz0)]*(Pml->P_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz-iz0)]);
                  }
              }
              if(Pml->getApplypml(5)){
                  i = iz + nzo - lpml;
                  if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->P_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->P_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i-iz0)];
                      Az[I2D(ix,i-iz0)] -= Rz[I2D(ix,i-iz0)]*(Pml->P_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i-iz0)]);
                  }
              }
          }
       }
    }
    
    
} // End of forwardstepAcceleration

template<typename T>
void WavesAcoustic2D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iz0, iz, nx, nz, lpml, nxo, nzo;
    T dt;
    lpml = model->getLpml();
    dt = this->getDt();
    T *L, *df;
    L = model->getL();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
       // Domain decomposition 
       nx = model->getNx();
       nz = model->getNz();
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
    }else{
       nx = model->getNx() + 2*lpml;
       nz = model->getNz() + 2*lpml;
       ix0 = 0;
       iz0 = 0;
       nxo = nx;
       nzo = nz;
    }
 
    
    // Derivate Ax backward with respect to x
    der->ddx_bw(Ax);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            P2[I2D(ix,iz)] = 2.0 * P1[I2D(ix,iz)] - P2[I2D(ix,iz)] + dt*dt*L[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
        for(iz=0; iz < nz; iz++){
            for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                    if(ix >= ix0 && ix < (ix0 + nx)){
                        // Left
                        Pml->Axx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Axx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix-ix0,iz)];
                        P2[I2D(ix-ix0,iz)] -= dt*dt*L[I2D(ix-ix0,iz)]*(Pml->Axx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix-ix0,iz)]);
                    }
                }
                if(Pml->getApplypml(1)){
                    i = ix + nxo - lpml;
                    if(i >= ix0 && i < (ix0 + nx)){
                        // Right
                        Pml->Axx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Axx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i-ix0,iz)];
                        P2[I2D(i-ix0,iz)] -= dt*dt*L[I2D(i-ix0,iz)]*(Pml->Axx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i-ix0,iz)]);
                    }
                }
            }
        }
    } 
    
    // Derivate Az backward with respect to z
    der->ddz_bw(Az);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            P2[I2D(ix,iz)] +=  dt*dt*L[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
        for(iz=0; iz < lpml; iz++){
            for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                    if(iz >= iz0 && iz < (iz0 + nz)){
                        // Top
                        Pml->Azz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Azz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz-iz0)];
                        P2[I2D(ix,iz-iz0)] -= dt*dt*L[I2D(ix,iz-iz0)]*(Pml->Azz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz-iz0)]);
                    }
                }
                if(Pml->getApplypml(5)){
                    i = iz + nzo - lpml;
                    if(i >= iz0 && i < (iz0 + nz)){
                        //Bottom
                        Pml->Azz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Azz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i-iz0)];
                        P2[I2D(ix,i-iz0)] -= dt*dt*L[I2D(ix,i-iz0)]*(Pml->Azz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i-iz0)]);
                    }
                }
            }
        }
    }
    
}

template<typename T>
void WavesAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, lpml;
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
    }
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2;
    switch(sourcetype)
    {
        case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    Ax[I2D(lpml + map[i].x, lpml + map[i].y)] += Mod[I2D(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    Az[I2D(lpml + map[i].x, lpml + map[i].y)] += Mod[I2D(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case PRESSURE:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    // Interpolation using sinc lanczos
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            P2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesAcoustic2D<T>::recordData(std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *dataarray; 
    T *Fielddata;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, lpml;

    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
    }

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i, i1, i2;
    switch(field)
    {
        case PRESSURE:
            Fielddata = this->getP1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
	    break;
        case VX:
            Fielddata = this->getAx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VZ:
            Fielddata = this->getAz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        default:
            break;
    }

}


// Roll the pressure pointers
template<typename T>
void WavesAcoustic2D<T>::roll()
{
    T *tmp;
    tmp = P2;
    P2 = P1;
    P1 = tmp;
}

// =============== 3D ACOUSTIC WAVES CLASS =============== //
template<typename T>
WavesAcoustic3D<T>::WavesAcoustic3D(){
    // Nothing here
}

template<typename T>
WavesAcoustic3D<T>::WavesAcoustic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot): Waves<T>(3, _nx, _ny, _nz, _nt, _lpml, _dx, _dy, _dz, _dt, _ox, _oy, _oz, _ot) {
    
    /* Create associated PML class */
    Pml =std::make_shared<PmlAcoustic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    P1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    P2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ax = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ay = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Az = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}


template<typename T>
WavesAcoustic3D<T>::WavesAcoustic3D(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;
    int nx_pml, ny_pml, nz_pml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    if((model->getDomain())->getStatus()){
        // Domain decomposition mode
        this->setDomain(model->getDomain());
        int ix0, nxo, iy0, nyo, iz0, nzo;
        bool low[3],high[3];
        for (int i=0; i<3; i++){
           low[i] = false;
           high[i] = false;
        }
        
        ix0 = (model->getDomain())->getIx0();
        iy0 = (model->getDomain())->getIy0();
        iz0 = (model->getDomain())->getIz0();
        nxo = (model->getDomain())->getNx_orig();
        nyo = (model->getDomain())->getNy_orig();
        nzo = (model->getDomain())->getNz_orig();
        if(ix0 < _lpml) low[0]=true;
        if(ix0 + _nx >= nxo-_lpml) high[0]=true;
        if(iy0 < _lpml) low[1]=true;
        if(iy0 + _ny >= nyo-_lpml) high[1]=true;
        if(iz0 < _lpml) low[2]=true;
        if(iz0 + _nz >= nzo-_lpml) high[2]=true;

        /* Create associated PML class */
        Pml = std::make_shared<PmlAcoustic3D<T>>(_nx, _ny, _nz, _lpml, _dt, &low[0], &high[0]);

        nx_pml = _nx;
        ny_pml = _ny;
        nz_pml = _nz;

        P1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        P2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Ax = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Ay = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Az = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    }else{
        /* Create associated PML class */
        Pml = std::make_shared<PmlAcoustic3D<T>>(_nx, _ny, _nz, _lpml, _dt);

        /* Allocate memory variables */
        this->setDim(3);
        nx_pml = _nx + 2*_lpml;
        ny_pml = _ny + 2*_lpml;
        nz_pml = _nz + 2*_lpml;
        P1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        P2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Ax = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Ay = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
        Az = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    }
}

template<typename T>
void WavesAcoustic3D<T>::forwardstepAcceleration(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iy0, iy, iz0, iz, nx, ny, nz, lpml, nxo, nyo, nzo;
    lpml = model->getLpml();
    T *Rx, *Ry, *Rz, *df;
    Rx = model->getRx();
    Ry = model->getRx();
    Rz = model->getRz();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
        // Domain decomposition 
        nx = model->getNx();
        ny = model->getNy();
        nz = model->getNz();
        ix0 = (model->getDomain())->getIx0();
        iy0 = (model->getDomain())->getIy0();
        iz0 = (model->getDomain())->getIz0();
        nxo = (model->getDomain())->getNx_orig();
        nyo = (model->getDomain())->getNy_orig();
        nzo = (model->getDomain())->getNz_orig();
    }else{
        nx = model->getNx() + 2*lpml;
        ny = model->getNy() + 2*lpml;
        nz = model->getNz() + 2*lpml;
        ix0 = 0;
        iy0 = 0;
        iz0 = 0;
        nxo = nx;
        nyo = ny;
        nzo = nz;
    }

    // Derivate P forward with respect to x
    der->ddx_fw(P1);
    // Compute Ax
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ax[I3D(ix,iy,iz)] = Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
        for(iz=0; iz < nz; iz++){
            for(iy=0; iy < ny; iy++){
                for(ix=0; ix < lpml; ix++){
                    if(Pml->getApplypml(0)){
                        if(ix >= ix0 && ix < (ix0 + nx)){
                            // Left
                            Pml->P_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->P_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)];

                            Ax[I3D(ix-ix0,iy,iz)] -= Rx[I3D(ix-ix0,iy,iz)]*(Pml->P_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)]);
                        }
                    }

                    if(Pml->getApplypml(1)){
                        i = ix + nxo - lpml;
                        if(i >= ix0 && i < (ix0 + nx)){
                        // Right
                            Pml->P_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->P_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)];
                            Ax[I3D(i-ix0,iy,iz)] -= Rx[I3D(i-ix0,iy,iz)]*(Pml->P_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)]);
                        }
                    }
                }
            }
        }
    }
    
    // Derivate P forward with respect to y
    der->ddy_fw(P1);
    // Compute Ay
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ay[I3D(ix,iy,iz)] = Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
        for(iz=0; iz < nz; iz++){
            for(iy=0; iy < lpml; iy++){
                for(ix=0; ix < nx; ix++){
                    if(Pml->getApplypml(2)){
                        if(iy >= iy0 && iy < (iy0 + ny)){
                            // Front
                            Pml->P_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->P_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)];
                            Ay[I3D(ix,iy-iy0,iz)] -= Ry[I3D(ix,iy-iy0,iz)]*(Pml->P_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)]);

                        }
                    }
                    if(Pml->getApplypml(3)){
                        i = iy + nyo - lpml;
                        if(i >= iy0 && i < (iy0 + ny)){
                            //Back
                            Pml->P_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->P_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)];
                            Ay[I3D(ix,i-iy0,iz)] -= Ry[I3D(ix,i-iy0,iz)]*(Pml->P_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)]);
                        }
                    }
                }
            }
        }
    }
    
    // Derivate P forward with respect to z
    der->ddz_fw(P1);
    // Compute Az
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Az[I3D(ix,iy,iz)] = Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate bottom and top using staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
        for(iz=0; iz < lpml; iz++){
            for(iy=0; iy < ny; iy++){
                for(ix=0; ix < nx; ix++){
                    if(Pml->getApplypml(4)){
                        if(iz >= iz0 && iz < (iz0 + nz)){
                            // Top
                            Pml->P_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->P_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)];
                            Az[I3D(ix,iy,iz-iz0)] -= Rz[I3D(ix,iy,iz-iz0)]*(Pml->P_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)]);
                        }
                    }

                    if(Pml->getApplypml(5)){
                        i = iz + nzo - lpml;
                        if(i >= iz0 && i < (iz0 + nz)){
                            //Bottom
                            Pml->P_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->P_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)];
                            Az[I3D(ix,iy,i-iz0)] -= Rz[I3D(ix,iy,i-iz0)]*(Pml->P_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)]);
                        }
                    }
                }
            }
        }
    }
} // End of forwardstepAcceleration

template<typename T>
void WavesAcoustic3D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iy0, iy, iz0, iz, nx, ny, nz, lpml, nxo, nyo, nzo;
    T dt;
    lpml = model->getLpml();
    dt = this->getDt();
    T *L, *df;
    L = model->getL();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
        // Domain decomposition 
        nx = model->getNx();
        ny = model->getNy();
        nz = model->getNz();
        ix0 = (model->getDomain())->getIx0();
        iy0 = (model->getDomain())->getIy0();
        iz0 = (model->getDomain())->getIz0();
        nxo = (model->getDomain())->getNx_orig();
        nyo = (model->getDomain())->getNy_orig();
        nzo = (model->getDomain())->getNz_orig();
    }else{
        nx = model->getNx() + 2*lpml;
        ny = model->getNy() + 2*lpml;
        nz = model->getNz() + 2*lpml;
        ix0 = 0;
        iy0 = 0;
        iz0 = 0;
        nxo = nx;
        nyo = ny;
        nzo = nz;
    }

    // Derivate Ax backward with respect to x
    der->ddx_bw(Ax);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I3D(ix,iy,iz)] = 2.0 * P1[I3D(ix,iy,iz)] - P2[I3D(ix,iy,iz)] + dt*dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
        for(iz=0; iz < nz; iz++){
            for(iy=0; iy < ny; iy++){
                for(ix=0; ix < lpml; ix++){
                    if(Pml->getApplypml(0)){
                        if(ix >= ix0 && ix < (ix0 + nx)){
                            // Left
                            Pml->Axx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Axx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix-ix0,iy,iz)];

                            P2[I3D(ix-ix0,iy,iz)] -= dt*dt*L[I3D(ix-ix0,iy,iz)]*(Pml->Axx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                        }
                    }
                    if(Pml->getApplypml(1)){
                        i = ix + nxo - lpml;
                        if(i >= ix0 && i < (ix0 + nx)){
                            // Right
                            Pml->Axx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Axx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i-ix0,iy,iz)];
                            P2[I3D(i-ix0,iy,iz)] -= dt*dt*L[I3D(i-ix0,iy,iz)]*(Pml->Axx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                        }
                    }
                }
            }
        }
    }
    
    // Derivate Ay backward with respect to y
    der->ddy_bw(Ay);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I3D(ix,iy,iz)] +=  dt*dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using non-staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
        for(iz=0; iz < nz; iz++){
            for(iy=0; iy < lpml; iy++){
                for(ix=0; ix < nx; ix++){
                    if(Pml->getApplypml(2)){
                        if(iy >= iy0 && iy < (iy0 + ny)){
                            // Front
                            Pml->Ayy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Ayy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy-iy0,iz)];
                            P2[I3D(ix,iy-iy0,iz)] -= dt*dt*L[I3D(ix,iy-iy0,iz)]*(Pml->Ayy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                        }
                    }

                    if(Pml->getApplypml(3)){
                        i = iy + nyo - lpml;
                        if(i >= iy0 && i < (iy0 + ny)){
                            //Back
                            Pml->Ayy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Ayy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i-iy0,iz)];
                            P2[I3D(ix,i-iy0,iz)] -= dt*dt*L[I3D(ix,i-iy0,iz)]*(Pml->Ayy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                        }
                    }
                }
            }
        }
    }
    
    
    // Derivate Az backward with respect to z
    der->ddz_bw(Az);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I3D(ix,iy,iz)] +=  dt*dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
        for(iz=0; iz < lpml; iz++){
            for(iy=0; iy < ny; iy++){
                for(ix=0; ix < nx; ix++){
                    if(Pml->getApplypml(4)){
                        if(iz >= iz0 && iz < (iz0 + nz)){
                            // Top
                            Pml->Azz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Azz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz-iz0)];

                            P2[I3D(ix,iy,iz-iz0)] -= dt*dt*L[I3D(ix,iy,iz-iz0)]*(Pml->Azz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                        }
                    }
                    if(Pml->getApplypml(5)){
                        i = iz + nzo - lpml;
                        if(i >= iz0 && i < (iz0 + nz)){
                            //Bottom
                            Pml->Azz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Azz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i-iz0)];
                            P2[I3D(ix,iy,i-iz0)] -= dt*dt*L[I3D(ix,iy,i-iz0)]*(Pml->Azz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
void WavesAcoustic3D<T>::insertSource(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, lpml;
    lpml = this->getLpml();
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
        ny = this->getNy();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
        ny = this->getNy() + 2*lpml;
    }
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP ) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }
    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2,i3;
    switch(sourcetype)
    {
        case PRESSURE:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    // Interpolation using sinc lanczos
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                P2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
	    case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Ax[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case VY:
            Mod = model->getRy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Ay[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Az[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I3D(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesAcoustic3D<T>::recordData(std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *dataarray; 
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, ny, lpml;
    lpml = this->getLpml();

    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
        ny = this->getNy();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
        ny = this->getNy() + 2*lpml;
    }


    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i,i1,i2,i3;
    switch(field)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += P1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        case VX:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Ax[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        case VY:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Ay[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;

        case VZ:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Az[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}


template<typename T>
WavesAcoustic3D<T>::~WavesAcoustic3D() {
    // Freein variables
    free(P1);
    free(P2);
    free(Ax);
    free(Ay);
    free(Az);
}

// Roll the pressure pointers
template<typename T>
void WavesAcoustic3D<T>::roll()
{
    T *tmp;
    tmp = P2;
    P2 = P1;
    P1 = tmp;
}

// =============== 2D ELASTIC VELOCITY-STRESS WAVES CLASS =============== //
/** The 2D Elastic WAVES model class
 *
 */

template<typename T>
WavesElastic2D<T>::WavesElastic2D()	///< Constructor
{
    // Do nothing
}


template<typename T>
WavesElastic2D<T>::WavesElastic2D(const int _nx, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot): Waves<T>(2, _nx, 1, _nz, _nt, _lpml, _dx, 1.0, _dz, _dt, _ox, 0.0, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesElastic2D<T>::WavesElastic2D(std::shared_ptr<rockseis::ModelElastic2D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);
   
    if((model->getDomain())->getStatus()){
       // Domain decomposition mode
       this->setDomain(model->getDomain());
       int ix0, nxo, iz0, nzo;
       bool low[2],high[2];

       for (int i=0; i<2; i++){
          low[i] = false;
          high[i] = false;
       }
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
       if(ix0 < _lpml) low[0]=true;
       if(ix0 + _nx >= nxo-_lpml) high[0]=true;
       if(iz0 < _lpml) low[1]=true;
       if(iz0 + _nz >= nzo-_lpml) high[1]=true;

       /* Create associated PML class */
       Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt, &low[0], &high[0]);

       int nx_pml, nz_pml;
       nx_pml = _nx;
       nz_pml = _nz;

       Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Vx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Vz = (T *) calloc(nx_pml*nz_pml,sizeof(T));

    }else{

       /* Create associated PML class */
       Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);

       /* Allocate memory variables */
       int nx_pml, nz_pml;
       this->setDim(2);
       nx_pml = _nx + 2*_lpml;
       nz_pml = _nz + 2*_lpml;
       Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Vx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
       Vz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    }
}


template<typename T>
void WavesElastic2D<T>::forwardstepVelocity(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iz0, iz, nx, nz, lpml, nxo, nzo;
    lpml = model->getLpml();
    T dt;
    dt = this->getDt();
    T *Rx, *Rz, *df;
    Rx = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    
    if((model->getDomain())->getStatus()){
       // Domain decomposition 
       nx = model->getNx();
       nz = model->getNz();
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
    }else{
       nx = model->getNx() + 2*lpml;
       nz = model->getNz() + 2*lpml;
       ix0 = 0;
       iz0 = 0;
       nxo = nx;
       nzo = nz;
    }

    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I2D(ix,iz)] += dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
       for(ix=0; ix < lpml; ix++){
          if(Pml->getApplypml(0)){
             if(ix >= ix0 && ix < (ix0 + nx)){
                // Left
                Pml->Sxx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix-ix0,iz)];
                Vx[I2D(ix-ix0,iz)] -= dt*Rx[I2D(ix-ix0,iz)]*(Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix-ix0,iz)]);
             }
          }

          if(Pml->getApplypml(1)){
             i = ix + nxo - lpml;
             if(i >= ix0 && i < (ix0 + nx)){
                // Right
                Pml->Sxx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i-ix0,iz)];
                Vx[I2D(i-ix0,iz)] -= dt*Rx[I2D(i-ix0,iz)]*(Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i-ix0,iz)]);
             }
          }
       }
    }
    
    
    // Derivate Sxz backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I2D(ix,iz)] += dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
       for(ix=0; ix < nx; ix++){
          if(Pml->getApplypml(4)){
             if(iz >= iz0 && iz < (iz0 + nz)){
                // Top
                Pml->Sxzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz-iz0)];
                Vx[I2D(ix,iz-iz0)] -= dt*Rx[I2D(ix,iz-iz0)]*(Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz-iz0)]);
             }
          }
          if(Pml->getApplypml(5)){
             i = iz + nzo - lpml;
             if(i >= iz0 && i < (iz0 + nz)){
                //Bottom
                Pml->Sxzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i-iz0)];
                Vx[I2D(ix,i-iz0)] -= dt*Rx[I2D(ix,i-iz0)]*(Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i-iz0)]);
             }
          }
       }
    }

    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I2D(ix,iz)] += dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
       for(ix=0; ix < lpml; ix++){
          if(Pml->getApplypml(0)){
             if(ix >= ix0 && ix < (ix0 + nx)){
                // Left
                Pml->Sxzx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix-ix0,iz)];
                Vz[I2D(ix-ix0,iz)] -= dt*Rz[I2D(ix-ix0,iz)]*(Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix-ix0,iz)]);
             }
          }
          if(Pml->getApplypml(1)){
             i = ix + nxo - lpml;
             if(i >= ix0 && i < (ix0 + nx)){
                // Right
                Pml->Sxzx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i-ix0,iz)];
                Vz[I2D(i-ix0,iz)] -= dt*Rz[I2D(i-ix0,iz)]*(Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i-ix0,iz)]);
             }
          }
       }
    }
    
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I2D(ix,iz)] += dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
           if(Pml->getApplypml(4)){
              if(iz >= iz0 && iz < (iz0 + nz)){
                 // Top
                 Pml->Szz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz-iz0)];
                 Vz[I2D(ix,iz-iz0)] -= dt*Rz[I2D(ix,iz-iz0)]*(Pml->Szz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz-iz0)]);
              }
           }
           if(Pml->getApplypml(5)){
              i = iz + nzo - lpml;
              if(i >= iz0 && i < (iz0 + nz)){
                 //Bottom
                 Pml->Szz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i-iz0)];
                 Vz[I2D(ix,i-iz0)] -= dt*Rz[I2D(ix,i-iz0)]*(Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i-iz0)]);
              }
           }
        }
    }
    
} // End of forwardstepVelocity


template<typename T>
void WavesElastic2D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iz0, iz, nx, nz, lpml, nxo, nzo;
    T dt;
    lpml = model->getLpml();
    dt = this->getDt();
    T *L, *L2M, *M, *df;
    L = model->getL();
    L2M = model->getL2M();
    M = model->getM();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
       // Domain decomposition 
       nx = model->getNx();
       nz = model->getNz();
       ix0 = (model->getDomain())->getIx0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nzo = (model->getDomain())->getNz_orig();
    }else{
       nx = model->getNx() + 2*lpml;
       nz = model->getNz() + 2*lpml;
       ix0 = 0;
       iz0 = 0;
       nxo = nx;
       nzo = nz;
    }

    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I2D(ix,iz)] += dt*L2M[I2D(ix,iz)]*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += dt*L[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    

    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
       for(ix=0; ix < lpml; ix++){
          if(Pml->getApplypml(0)){
             if(ix >= ix0 && ix < (ix0 + nx)){
                // Left
                Pml->Vxx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix-ix0,iz)];
                Sxx[I2D(ix-ix0,iz)] -= dt*L2M[I2D(ix-ix0,iz)]*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix-ix0,iz)]);
                Szz[I2D(ix-ix0,iz)] -= dt*L[I2D(ix-ix0,iz)]*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix-ix0,iz)]);
             }
          }
          if(Pml->getApplypml(1)){
             i = ix + nxo - lpml;
             if(i >= ix0 && i < (ix0 + nx)){
                // Right
                Pml->Vxx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i-ix0,iz)];
                Sxx[I2D(i-ix0,iz)] -= dt*L2M[I2D(i-ix0,iz)]*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i-ix0,iz)]);
                Szz[I2D(i-ix0,iz)] -= dt*L[I2D(i-ix0,iz)]*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i-ix0,iz)]);
             }
          }
       }
    }
    
    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I2D(ix,iz)] += dt*L[I2D(ix,iz)]*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += dt*L2M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
       for(ix=0; ix < nx; ix++){
          if(Pml->getApplypml(4)){
             if(iz >= iz0 && iz < (iz0 + nz)){
                // Top
                Pml->Vzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz-iz0)];
                Sxx[I2D(ix,iz-iz0)] -= dt*L[I2D(ix,iz-iz0)]*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz-iz0)]);
                Szz[I2D(ix,iz-iz0)] -= dt*L2M[I2D(ix,iz-iz0)]*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz-iz0)]);
             }
          }
          if(Pml->getApplypml(5)){
             i = iz + nzo - lpml;
             if(i >= iz0 && i < (iz0 + nz)){
                //Bottom
                Pml->Vzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i-iz0)];
                Sxx[I2D(ix,i-iz0)] -= dt*L[I2D(ix,i-iz0)]*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i-iz0)]);
                Szz[I2D(ix,i-iz0)] -= dt*L2M[I2D(ix,i-iz0)]*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i-iz0)]);
             }
          }
       }
    }

    // Free surface conditions
    if(model->getFs()){
       iz = lpml;
       if(iz >= iz0 && iz < (iz0 + nz)){
          for(ix=0; ix < nx; ix++){
             Szz[I2D(ix,iz-iz0)] = 0.0;
          }
       }
    }
   
    // Derivate Vz forward with respect to x
    der->ddx_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I2D(ix,iz)] += dt*M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
       for(ix=0; ix < lpml; ix++){
          if(Pml->getApplypml(0)){
             if(ix >= ix0 && ix < (ix0 + nx)){
                // Left
                Pml->Vzx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix-ix0,iz)];
                Sxz[I2D(ix-ix0,iz)] -= dt*M[I2D(ix-ix0,iz)]*(Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix-ix0,iz)]);
             }
          }
          if(Pml->getApplypml(1)){
             i = ix + nxo - lpml;
             if(i >= ix0 && i < (ix0 + nx)){
                // Right
                Pml->Vzx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i-ix0,iz)];
                Sxz[I2D(i-ix0,iz)] -= dt*M[I2D(i-ix0,iz)]*(Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i-ix0,iz)]);
             }
          }
       }
    }
    
    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I2D(ix,iz)] += dt*M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
       for(ix=0; ix < nx; ix++){
          if(Pml->getApplypml(4)){
             if(iz >= iz0 && iz < (iz0 + nz)){
                // Top
                Pml->Vxz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz-iz0)];
                Sxz[I2D(ix,iz-iz0)] -= dt*M[I2D(ix,iz-iz0)]*(Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz-iz0)]);
             }
          }
          if(Pml->getApplypml(5)){
             i = iz + nzo - lpml;
             if(i >= iz0 && i < (iz0 + nz)){
                //Bottom
                Pml->Vxz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i-iz0)];
                Sxz[I2D(ix,i-iz0)] -= dt*M[I2D(ix,i-iz0)]*(Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i-iz0)]);
             }
          }
       }
    }
}

template<typename T>
void WavesElastic2D<T>::insertSource(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, lpml;
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
    }
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i, i1, i2;
    switch(sourcetype)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Sxx[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            Szz[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Vx[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Vz[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}


template<typename T>
void WavesElastic2D<T>::recordData(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, lpml;
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
    }
    rs_field field = data->getField();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    dataarray = data->getData();
    int i,i1,i2;
    Index Im(this->getNx(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {

                    // Factor to make reciprocity work
                    factor = R[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)] - R[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)];
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                        }
                    }
                }
            }
	    break;
        case VX:
            Fielddata1 = this->getVx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VZ:
            Fielddata1 = this->getVz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}


template<typename T>
WavesElastic2D<T>::~WavesElastic2D() {
    /* Free allocated variables */
    free(Sxx);
    free(Szz);
    free(Sxz);
    free(Vx);
    free(Vz);
}

// =============== 3D ELASTIC WAVES CLASS =============== //
/** The 3D Elastic WAVES model class
 *
 */

template<typename T>
WavesElastic3D<T>::WavesElastic3D(){
    // Nothing here
}

template<typename T>
WavesElastic3D<T>::WavesElastic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot): Waves<T>(3, _nx, _ny, _nz, _nt, _lpml, _dx, _dy, _dz, _dt, _ox, _oy, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesElastic3D<T>::WavesElastic3D(std::shared_ptr<rockseis::ModelElastic3D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    if((model->getDomain())->getStatus()){
       // Domain decomposition mode
       this->setDomain(model->getDomain());
       int ix0, nxo, iy0, nyo, iz0, nzo;
       bool low[3],high[3];
       for (int i=0; i<3; i++){
          low[i] = false;
          high[i] = false;
       }

       ix0 = (model->getDomain())->getIx0();
       iy0 = (model->getDomain())->getIy0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nyo = (model->getDomain())->getNy_orig();
       nzo = (model->getDomain())->getNz_orig();
       if(ix0 < _lpml) low[0]=true;
       if(ix0 + _nx >= nxo-_lpml) high[0]=true;
       if(iy0 < _lpml) low[1]=true;
       if(iy0 + _ny >= nyo-_lpml) high[1]=true;
       if(iz0 < _lpml) low[2]=true;
       if(iz0 + _nz >= nzo-_lpml) high[2]=true;

       /* Create associated PML class */
       Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt, &low[0], &high[0]);

       /* Allocate memory variables */
       int nx_pml, ny_pml, nz_pml;
       nx_pml = _nx;
       ny_pml = _ny;
       nz_pml = _nz;

       Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));

    }else{

       /* Create associated PML class */
       Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);

       /* Allocate memory variables */
       int nx_pml, ny_pml, nz_pml;
       nx_pml = _nx + 2*_lpml;
       ny_pml = _ny + 2*_lpml;
       nz_pml = _nz + 2*_lpml;
       Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
       Vz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    }
}

template<typename T>
WavesElastic3D<T>::~WavesElastic3D() {
    /* Free allocated variables */
    free(Sxx);
    free(Syy);
    free(Szz);
    free(Sxz);
    free(Syz);
    free(Sxy);
    free(Vx);
    free(Vy);
    free(Vz);
}

template<typename T>
void WavesElastic3D<T>::forwardstepVelocity(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iy0, iy, iz0, iz, nx, ny, nz, lpml, nxo, nyo, nzo;
    lpml = model->getLpml();
    T *Rx, *Ry, *Rz, *df;
    T dt;
    dt = this->getDt();
    Rx = model->getRx();
    Ry = model->getRx();
    Rz = model->getRz();
    df = der->getDf();

if((model->getDomain())->getStatus()){
        // Domain decomposition 
        nx = model->getNx();
        ny = model->getNy();
        nz = model->getNz();
        ix0 = (model->getDomain())->getIx0();
        iy0 = (model->getDomain())->getIy0();
        iz0 = (model->getDomain())->getIz0();
        nxo = (model->getDomain())->getNx_orig();
        nyo = (model->getDomain())->getNy_orig();
        nzo = (model->getDomain())->getNz_orig();
    }else{
        nx = model->getNx() + 2*lpml;
        ny = model->getNy() + 2*lpml;
        nz = model->getNz() + 2*lpml;
        ix0 = 0;
        iy0 = 0;
        iz0 = 0;
        nxo = nx;
        nyo = ny;
        nzo = nz;
    }
    
   //////////////////////////// VX //////////////////////////
   //
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                   if(ix >= ix0 && ix < (ix0 + nx)){
                      // Left
                      Pml->Sxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)];

                      Vx[I3D(ix-ix0,iy,iz)] -= dt*Rx[I3D(ix-ix0,iy,iz)]*(Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)]);
                   }
                }

                if(Pml->getApplypml(1)){
                   i = ix + nxo - lpml;
                   if(i >= ix0 && i < (ix0 + nx)){
                   // Right
                      Pml->Sxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)];
                      Vx[I3D(i-ix0,iy,iz)] -= dt*Rx[I3D(i-ix0,iy,iz)]*(Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)]);
                   }
                }
             }
          }
       }
    }

    // Derivate Sxy backward with respect to y
    der->ddy_bw(Sxy);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Sxyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy-iy0,iz)];

                      Vx[I3D(ix,iy-iy0,iz)] -= dt*Rx[I3D(ix,iy-iy0,iz)]*(Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Sxyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i-iy0,iz)];
                      Vx[I3D(ix,i-iy0,iz)] -= dt*Rx[I3D(ix,i-iy0,iz)]*(Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }
 

    // Derivate Sxy backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                   if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->Sxzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz-iz0)];

                      Vx[I3D(ix,iy,iz-iz0)] -= dt*Rx[I3D(ix,iy,iz-iz0)]*(Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                   }
                }
                if(Pml->getApplypml(5)){
                   i = iz + nzo - lpml;
                   if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i-iz0)];
                      Vx[I3D(ix,iy,i-iz0)] -= dt*Rx[I3D(ix,iy,i-iz0)]*(Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                   }
                }
             }
          }
       }
    }
 
   
   //////////////////////////// VY //////////////////////////
    

    // Derivate Sxy backward with respect to x
    der->ddx_bw(Sxy);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using non-staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                   if(ix >= ix0 && ix < (ix0 + nx)){
                      // Left
                      Pml->Sxyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix-ix0,iy,iz)];

                      Vy[I3D(ix-ix0,iy,iz)] -= dt*Ry[I3D(ix-ix0,iy,iz)]*(Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                   }
                }
                if(Pml->getApplypml(1)){
                   i = ix + nxo - lpml;
                   if(i >= ix0 && i < (ix0 + nx)){
                      // Right
                      Pml->Sxyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i-ix0,iy,iz)];
                      Vy[I3D(i-ix0,iy,iz)] -= dt*Ry[I3D(i-ix0,iy,iz)]*(Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                   }
                }
             }
          }
       }
    }


    // Derivate Syy forward with respect to y
    der->ddy_fw(Syy);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Syy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)];

                      Vy[I3D(ix,iy-iy0,iz)] -= dt*Ry[I3D(ix,iy-iy0,iz)]*(Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Syy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)];
                      Vy[I3D(ix,i-iy0,iz)] -= dt*Ry[I3D(ix,i-iy0,iz)]*(Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }
 

    // Derivate Syz backward with respect to z
    der->ddz_bw(Syz);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                   if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->Syzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz-iz0)];

                      Vy[I3D(ix,iy,iz-iz0)] -= dt*Ry[I3D(ix,iy,iz-iz0)]*(Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                   }
                }
                if(Pml->getApplypml(5)){
                   i = iz + nzo - lpml;
                   if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i-iz0)];
                      Vy[I3D(ix,iy,i-iz0)] -= dt*Ry[I3D(ix,iy,i-iz0)]*(Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                   }
                }
             }
          }
       }
    }
 
 
   //////////////////////////// VZ //////////////////////////
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iz=0; iz < nz; iz++){
             for(iy=0; iy < ny; iy++){
                for(ix=0; ix < lpml; ix++){
                   if(Pml->getApplypml(0)){
                      if(ix >= ix0 && ix < (ix0 + nx)){
                         // Left
                         Pml->Sxzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix-ix0,iy,iz)];

                         Vz[I3D(ix-ix0,iy,iz)] -= dt*Rz[I3D(ix-ix0,iy,iz)]*(Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                      }
                   }
                   if(Pml->getApplypml(1)){
                      i = ix + nxo - lpml;
                      if(i >= ix0 && i < (ix0 + nx)){
                         // Right
                         Pml->Sxzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i-ix0,iy,iz)];
                         Vz[I3D(i-ix0,iy,iz)] -= dt*Rz[I3D(i-ix0,iy,iz)]*(Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                      }
                   }
                }
             }
          }
       }
    }

    // Derivate Syz backward with respect to y
    der->ddy_bw(Syz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Syzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy-iy0,iz)];

                      Vz[I3D(ix,iy-iy0,iz)] -= dt*Rz[I3D(ix,iy-iy0,iz)]*(Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Syzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i-iy0,iz)];
                      Vz[I3D(ix,i-iy0,iz)] -= dt*Rz[I3D(ix,i-iy0,iz)]*(Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }

    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                   if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->Szz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)];
                      Vz[I3D(ix,iy,iz-iz0)] -= dt*Rz[I3D(ix,iy,iz-iz0)]*(Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)]);
                   }
                }
                if(Pml->getApplypml(5)){
                   i = iz + nzo - lpml;
                   if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->Szz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)];
                      Vz[I3D(ix,iy,i-iz0)] -= dt*Rz[I3D(ix,iy,i-iz0)]*(Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)]);
                   }
                }
             }
          }
       }
    }

} // End of forwardstepVelocity

template<typename T>
void WavesElastic3D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix0, ix, iy0, iy, iz0, iz, nx, ny, nz, lpml, nxo, nyo, nzo;
    T dt;
    lpml = model->getLpml();
    dt = this->getDt();
    T *L, *L2M, *M_xz, *M_yz, *M_xy, *df;
    L = model->getL();
    L2M = model->getL2M();
    M_xy = model->getM_xy();
    M_xz = model->getM_xz();
    M_yz = model->getM_yz();
    df = der->getDf();

    if((model->getDomain())->getStatus()){
       // Domain decomposition 
       nx = model->getNx();
       ny = model->getNy();
       nz = model->getNz();
       ix0 = (model->getDomain())->getIx0();
       iy0 = (model->getDomain())->getIy0();
       iz0 = (model->getDomain())->getIz0();
       nxo = (model->getDomain())->getNx_orig();
       nyo = (model->getDomain())->getNy_orig();
       nzo = (model->getDomain())->getNz_orig();
    }else{
       nx = model->getNx() + 2*lpml;
       ny = model->getNy() + 2*lpml;
       nz = model->getNz() + 2*lpml;
       ix0 = 0;
       iy0 = 0;
       iz0 = 0;
       nxo = nx;
       nyo = ny;
       nzo = nz;
    }

    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] += dt*L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] += dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] += dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }


    // Attenuate left and right using non-staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                   if(ix >= ix0 && ix < (ix0 + nx)){
                      // Left
                      Pml->Vxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix-ix0,iy,iz)];

                      Sxx[I3D(ix-ix0,iy,iz)] -= dt*L2M[I3D(ix-ix0,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                      Syy[I3D(ix-ix0,iy,iz)] -= dt*L[I3D(ix-ix0,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                      Szz[I3D(ix-ix0,iy,iz)] -= dt*L[I3D(ix-ix0,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix-ix0,iy,iz)]);
                   }
                }
                if(Pml->getApplypml(1)){
                   i = ix + nxo - lpml;
                   if(i >= ix0 && i < (ix0 + nx)){
                      // Right
                      Pml->Vxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i-ix0,iy,iz)];
                      Sxx[I3D(i-ix0,iy,iz)] -= dt*L2M[I3D(i-ix0,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                      Syy[I3D(i-ix0,iy,iz)] -= dt*L[I3D(i-ix0,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                      Szz[I3D(i-ix0,iy,iz)] -= dt*L[I3D(i-ix0,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i-ix0,iy,iz)]);
                   }
                }
             }
          }
       }
    }
    
    // Derivate Vy backward with respect to y
    der->ddy_bw(Vy);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] +=  dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] +=  dt*L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] +=  dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using non-staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Vyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy-iy0,iz)];

                      Sxx[I3D(ix,iy-iy0,iz)] -= dt*L[I3D(ix,iy-iy0,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                      Syy[I3D(ix,iy-iy0,iz)] -= dt*L2M[I3D(ix,iy-iy0,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                      Szz[I3D(ix,iy-iy0,iz)] -= dt*L[I3D(ix,iy-iy0,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Vyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i-iy0,iz)];
                      Sxx[I3D(ix,i-iy0,iz)] -= dt*L[I3D(ix,i-iy0,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                      Syy[I3D(ix,i-iy0,iz)] -= dt*L2M[I3D(ix,i-iy0,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                      Szz[I3D(ix,i-iy0,iz)] -= dt*L[I3D(ix,i-iy0,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }

    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] +=  dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] +=  dt*L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] +=  dt*L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
       iz = lpml;
       if(iz >= iz0 && iz < (iz0 + nz)){
          for(ix=0; ix < nx; ix++){
             for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz-iz0)] = 0.0;
             }
          }
       }
    }

    // Attenuate top and bottom using non-staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                   if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->Vzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz-iz0)];

                      Sxx[I3D(ix,iy,iz-iz0)] -= dt*L[I3D(ix,iy,iz-iz0)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                      Syy[I3D(ix,iy,iz-iz0)] -= dt*L[I3D(ix,iy,iz-iz0)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                      Szz[I3D(ix,iy,iz-iz0)] -= dt*L2M[I3D(ix,iy,iz-iz0)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz-iz0)]);
                   }
                }
                if(Pml->getApplypml(5)){
                   i = iz + nzo - lpml;
                   if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i-iz0)];
                      Sxx[I3D(ix,iy,i-iz0)] -= dt*L[I3D(ix,iy,i-iz0)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                      Syy[I3D(ix,iy,i-iz0)] -= dt*L[I3D(ix,iy,i-iz0)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                      Szz[I3D(ix,iy,i-iz0)] -= dt*L2M[I3D(ix,iy,i-iz0)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i-iz0)]);
                   }
                }
             }
          }
       }
    }
    
   //////////////////////////// Sxz //////////////////////////
  
    // Derivate Vz forward with respect to x
    der->ddx_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxz[I3D(ix,iy,iz)] +=  dt*M_xz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                   if(ix >= ix0 && ix < (ix0 + nx)){
                      // Left
                      Pml->Vzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)];

                      Sxz[I3D(ix-ix0,iy,iz)] -= dt*M_xz[I3D(ix-ix0,iy,iz)]*(Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)]);
                   }
                }
                if(Pml->getApplypml(1)){
                   i = ix + nxo - lpml;
                   if(i >= ix0 && i < (ix0 + nx)){
                      // Right
                      Pml->Vzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)];
                      Sxz[I3D(i-ix0,iy,iz)] -= dt*M_xz[I3D(i-ix0,iy,iz)]*(Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)]);
                   }
                }
             }
          }
       }
    }
 

    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxz[I3D(ix,iy,iz)] +=  dt*M_xz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iz=0; iz < lpml; iz++){
             for(iy=0; iy < ny; iy++){
                for(ix=0; ix < nx; ix++){
                   if(Pml->getApplypml(4)){
                      if(iz >= iz0 && iz < (iz0 + nz)){
                         // Top
                         Pml->Vxz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)];

                         Sxz[I3D(ix,iy,iz-iz0)] -= dt*M_xz[I3D(ix,iy,iz-iz0)]*(Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)]);
                      }
                   }
                   if(Pml->getApplypml(5)){
                      i = iz + nzo - lpml;
                      if(i >= iz0 && i < (iz0 + nz)){
                         //Bottom
                         Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)];
                         Sxz[I3D(ix,iy,i-iz0)] -= dt*M_xz[I3D(ix,iy,i-iz0)]*(Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)]);
                      }
                   }
                }
             }
          }
       }
    }


   //////////////////////////// Syz //////////////////////////
    // Derivate Vy forward with respect to z
    der->ddz_fw(Vy);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
       for(iy=0; iy < ny; iy++){
          for(ix=0; ix < nx; ix++){
             Syz[I3D(ix,iy,iz)] +=  dt*M_yz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
          }
       }
    }

    // Attenuate top and bottom using staggered variables
    if(Pml->getApplypml(4) || Pml->getApplypml(5)){
       for(iz=0; iz < lpml; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(4)){
                   if(iz >= iz0 && iz < (iz0 + nz)){
                      // Top
                      Pml->Vyz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)];

                      Syz[I3D(ix,iy,iz-iz0)] -= dt*M_yz[I3D(ix,iy,iz-iz0)]*(Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz-iz0)]);
                   }
                }
                if(Pml->getApplypml(5)){
                   i = iz + nzo - lpml;
                   if(i >= iz0 && i < (iz0 + nz)){
                      //Bottom
                      Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)];
                      Syz[I3D(ix,iy,i-iz0)] -= dt*M_yz[I3D(ix,iy,i-iz0)]*(Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i-iz0)]);
                   }
                }
             }
          }
       }
    }

    // Derivate Vz forward with respect to y
    der->ddy_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Syz[I3D(ix,iy,iz)] +=  dt*M_yz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Vzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)];

                      Syz[I3D(ix,iy-iy0,iz)] -= dt*M_yz[I3D(ix,iy-iy0,iz)]*(Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Vzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)];
                      Syz[I3D(ix,i-iy0,iz)] -= dt*M_yz[I3D(ix,i-iy0,iz)]*(Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }


   //////////////////////////// Sxy //////////////////////////
    // Derivate Vx forward with respect to y
    der->ddy_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxy[I3D(ix,iy,iz)] +=  dt*M_xy[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    if(Pml->getApplypml(2) || Pml->getApplypml(3)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < lpml; iy++){
             for(ix=0; ix < nx; ix++){
                if(Pml->getApplypml(2)){
                   if(iy >= iy0 && iy < (iy0 + ny)){
                      // Front
                      Pml->Vxy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)];

                      Sxy[I3D(ix,iy-iy0,iz)] -= dt*M_xy[I3D(ix,iy-iy0,iz)]*(Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy-iy0,iz)]);
                   }
                }
                if(Pml->getApplypml(3)){
                   i = iy + nyo - lpml;
                   if(i >= iy0 && i < (iy0 + ny)){
                      //Back
                      Pml->Vxy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)];
                      Sxy[I3D(ix,i-iy0,iz)] -= dt*M_xy[I3D(ix,i-iy0,iz)]*(Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i-iy0,iz)]);
                   }
                }
             }
          }
       }
    }

    // Derivate Vy forward with respect to x
    der->ddx_fw(Vy);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxy[I3D(ix,iy,iz)] +=  dt*M_xy[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    if(Pml->getApplypml(0) || Pml->getApplypml(1)){
       for(iz=0; iz < nz; iz++){
          for(iy=0; iy < ny; iy++){
             for(ix=0; ix < lpml; ix++){
                if(Pml->getApplypml(0)){
                   if(ix >= ix0 && ix < (ix0 + nx)){
                      // Left
                      Pml->Vyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)];

                      Sxy[I3D(ix-ix0,iy,iz)] -= dt*M_xy[I3D(ix-ix0,iy,iz)]*(Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix-ix0,iy,iz)]);
                   }
                }
                if(Pml->getApplypml(1)){
                   i = ix + nxo - lpml;
                   if(i >= ix0 && i < (ix0 + nx)){
                      // Right
                      Pml->Vyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)];
                      Sxy[I3D(i-ix0,iy,iz)] -= dt*M_xy[I3D(i-ix0,iy,iz)]*(Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i-ix0,iy,iz)]);
                   }
                }
             }
          }
       }
    }
}

template<typename T>
void WavesElastic3D<T>::insertSource(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, lpml;
    lpml = this->getLpml();
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
        ny = this->getNy();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
        ny = this->getNy() + 2*lpml;
    }

    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }
    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2,i3;
    T xtri[2], ytri[2], ztri[2];
    switch(sourcetype)
    {
       case PRESSURE:
          for (i=0; i < ntrace; i++) 
          { 
             if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
             {

                xtri[0] = 1 - shift[i].x;
                xtri[1] = shift[i].x;

                ytri[0] = 1 - shift[i].y;
                ytri[1] = shift[i].y;

                ztri[0] = 1 - shift[i].z;
                ztri[1] = shift[i].z;

                for(i1 = 0; i1 < 2; i1++){ 
                   for(i2 = 0; i2 < 2; i2++){ 
                      for(i3 = 0; i3 < 2; i3++){ 
                         Sxx[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] -= dt*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                         Syy[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] -= dt*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                         Szz[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] -= dt*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                      }
                   }
                }
             }
          }
          break;
       case VX:
          Mod = model->getRx();
          for (i=0; i < ntrace; i++) 
          { 
             if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
             {

                xtri[0] = 1 - shift[i].x;
                xtri[1] = shift[i].x;

                ytri[0] = 1 - shift[i].y;
                ytri[1] = shift[i].y;

                ztri[0] = 1 - shift[i].z;
                ztri[1] = shift[i].z;
                for(i1=0; i1<2; i1++){
                   for(i2=0; i2<2; i2++){
                      for(i3=0; i3<2; i3++){
                         Vx[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] += dt*Mod[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                      }
                   }
                }
             }
          }
          break;
       case VY:
          Mod = model->getRy();
          for (i=0; i < ntrace; i++) 
          { 
             if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
             {

                xtri[0] = 1 - shift[i].x;
                xtri[1] = shift[i].x;

                ytri[0] = 1 - shift[i].y;
                ytri[1] = shift[i].y;

                ztri[0] = 1 - shift[i].z;
                ztri[1] = shift[i].z;
                for(i1=0; i1<2; i1++){
                   for(i2=0; i2<2; i2++){
                      for(i3=0; i3<2; i3++){
                         Vy[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] += dt*Mod[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                      }
                   }
                }
             }
          }
          break;
       case VZ:
          Mod = model->getRz();
          for (i=0; i < ntrace; i++) 
          { 
             if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
             {

                xtri[0] = 1 - shift[i].x;
                xtri[1] = shift[i].x;

                ytri[0] = 1 - shift[i].y;
                ytri[1] = shift[i].y;

                ztri[0] = 1 - shift[i].z;
                ztri[1] = shift[i].z;
                for(i1=0; i1<2; i1++){
                   for(i2=0; i2<2; i2++){
                      for(i3=0; i3<2; i3++){
                         Vz[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)] += dt*Mod[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]*wav[Idat(it,i)]*xtri[i3]*ytri[i2]*ztri[i1];
                      }
                   }
                }
             }
          }
          break;
       default:
            break;
    }
}

template<typename T>
void WavesElastic3D<T>::recordData(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Fielddata3 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, ny, lpml;
    lpml = this->getLpml();
    if(this->getDomdec()){
        lpml = 0;
        nx = this->getNx();
        ny = this->getNy();
    }else{
        lpml = this->getLpml();
        nx = this->getNx() + 2*lpml;
        ny = this->getNy() + 2*lpml;
    }

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i,i1,i2,i3;
    T xtri[2], ytri[2], ztri[2];
    Index Im(this->getNx(), this->getNy(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSyy();
            Fielddata3 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {

                   // Factor to make reciprocity work
                   factor = (3.*R[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)] - 4.*R[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)])/3.;
                   xtri[0] = 1 - shift[i].x;
                   xtri[1] = shift[i].x;

                   ytri[0] = 1 - shift[i].y;
                   ytri[1] = shift[i].y;

                   ztri[0] = 1 - shift[i].z;
                   ztri[1] = shift[i].z;

                   for(i1 = 0; i1 < 2; i1++){ 
                      for(i2 = 0; i2 < 2; i2++){ 
                         for(i3 = 0; i3 < 2; i3++){ 
                            dataarray[Idat(it,i)] += (1./3.)*xtri[i3]*ytri[i2]*ztri[i1]*Fielddata1[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]/factor;
                            dataarray[Idat(it,i)] += (1./3.)*xtri[i3]*ytri[i2]*ztri[i1]*Fielddata2[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]/factor;
                            dataarray[Idat(it,i)] += (1./3.)*xtri[i3]*ytri[i2]*ztri[i1]*Fielddata3[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)]/factor;
                         }
                      }
                   }
                }
            }
            break;
        case VX:
            Fielddata1 = this->getVx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                   xtri[0] = 1 - shift[i].x;
                   xtri[1] = shift[i].x;

                   ytri[0] = 1 - shift[i].y;
                   ytri[1] = shift[i].y;

                   ztri[0] = 1 - shift[i].z;
                   ztri[1] = shift[i].z;

                   for(i1 = 0; i1 < 2; i1++){ 
                      for(i2 = 0; i2 < 2; i2++){ 
                         for(i3 = 0; i3 < 2; i3++){ 
                            dataarray[Idat(it,i)] += xtri[i3]*ytri[i2]*ztri[i1]*Fielddata1[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)];
                         }
                      }
                   }
                }
            }
            break;
        case VY:
            Fielddata1 = this->getVy();
            for (i=0; i < ntrace; i++) 
            { 
               if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
               {

                  xtri[0] = 1 - shift[i].x;
                  xtri[1] = shift[i].x;

                  ytri[0] = 1 - shift[i].y;
                  ytri[1] = shift[i].y;

                  ztri[0] = 1 - shift[i].z;
                  ztri[1] = shift[i].z;

                  for(i1 = 0; i1 < 2; i1++){ 
                     for(i2 = 0; i2 < 2; i2++){ 
                        for(i3 = 0; i3 < 2; i3++){ 
                           dataarray[Idat(it,i)] += xtri[i3]*ytri[i2]*ztri[i1]*Fielddata1[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)];
                        }
                     }
                  }
               }
            }
            break;

        case VZ:
            Fielddata1 = this->getVz();
            for (i=0; i < ntrace; i++) 
            { 
               if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
               {

                  xtri[0] = 1 - shift[i].x;
                  xtri[1] = shift[i].x;

                  ytri[0] = 1 - shift[i].y;
                  ytri[1] = shift[i].y;

                  ztri[0] = 1 - shift[i].z;
                  ztri[1] = shift[i].z;

                  for(i1 = 0; i1 < 2; i1++){ 
                     for(i2 = 0; i2 < 2; i2++){ 
                        for(i3 = 0; i3 < 2; i3++){ 
                           dataarray[Idat(it,i)] += xtri[i3]*ytri[i2]*ztri[i1]*Fielddata1[I3D(lpml + map[i].x + i3, lpml + map[i].y + i2, lpml + map[i].z + i1)];
                        }
                     }
                  }
               }
            }
            break;
        default:
            break;
    }
}


// =============== 2D ELASTIC DISPLACEMENT-STRESS WAVES CLASS =============== //
/** The 2D Elastic WAVES model class
 *
 */

template<typename T>
WavesElastic2D_DS<T>::WavesElastic2D_DS()	///< Constructor
{
    // Do nothing
}

template<typename T>
WavesElastic2D_DS<T>::WavesElastic2D_DS(const int _nx, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot): Waves<T>(2, _nx, 1, _nz, _nt, _lpml, _dx, 1.0, _dz, _dt, _ox, 0.0, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ux1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ux2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Uz1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Uz2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesElastic2D_DS<T>::WavesElastic2D_DS(std::shared_ptr<rockseis::ModelElastic2D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);
   
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ux1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Ux2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Uz1 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Uz2 = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}


template<typename T>
void WavesElastic2D_DS<T>::forwardstepDisplacement(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *Rx, *Rz, *df;
    Rx = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Ux
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Ux2[I2D(ix,iz)] = 2.0*Ux1[I2D(ix,iz)] - Ux2[I2D(ix,iz)] + dt*dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix,iz)];
            Ux2[I2D(ix,iz)] -= dt*dt*Rx[I2D(ix,iz)]*(Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i,iz)];
            Ux2[I2D(i,iz)] -= dt*dt*Rx[I2D(i,iz)]*(Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i,iz)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Ux
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Ux2[I2D(ix,iz)] += dt*dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Sxzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz)];
            
            Ux2[I2D(ix,iz)] -= dt*dt*Rx[I2D(ix,iz)]*(Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Sxzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i)];
            Ux2[I2D(ix,i)] -= dt*dt*Rx[I2D(ix,i)]*(Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Uz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Uz2[I2D(ix,iz)] = 2.0*Uz1[I2D(ix,iz)] - Uz2[I2D(ix,iz)] +dt*dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxzx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix,iz)];
            Uz2[I2D(ix,iz)] -= dt*dt*Rz[I2D(ix,iz)]*(Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxzx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i,iz)];
            Uz2[I2D(i,iz)] -= dt*dt*Rz[I2D(i,iz)]*(Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Uz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Uz2[I2D(ix,iz)] += dt*dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Szz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz)];
            
            Uz2[I2D(ix,iz)] -= dt*dt*Rz[I2D(ix,iz)]*(Pml->Szz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Szz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i)];
            Uz2[I2D(ix,i)] -= dt*dt*Rz[I2D(ix,i)]*(Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i)]);
        }
    }
    
} // End of forwardstepDisplacement


template<typename T>
void WavesElastic2D_DS<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *L, *L2M, *M, *df;
    L = model->getL();
    L2M = model->getL2M();
    M = model->getM();
    df = der->getDf();
    
    // Derivate Ux backward with respect to x
    der->ddx_bw(Ux1);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I2D(ix,iz)] = L2M[I2D(ix,iz)]*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] = L[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }


    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            Szz[I2D(ix,iz)] = 0.0;
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Vxx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix,iz)];
            
            Sxx[I2D(ix,iz)] -= L2M[I2D(ix,iz)]*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            Szz[I2D(ix,iz)] -= L[I2D(ix,iz)]*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Vxx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i,iz)];
            Sxx[I2D(i,iz)] -= L2M[I2D(i,iz)]*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
            Szz[I2D(i,iz)] -= L[I2D(i,iz)]*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Uz backward with respect to z
    der->ddz_bw(Uz1);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I2D(ix,iz)] += L[I2D(ix,iz)]*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += L2M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            Szz[I2D(ix,iz)] = 0.0;
        }
    }

    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Vzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz)];
            
            Sxx[I2D(ix,iz)] -= L[I2D(ix,iz)]*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            Szz[I2D(ix,iz)] -= L2M[I2D(ix,iz)]*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Vzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i)];
            Sxx[I2D(ix,i)] -= L[I2D(ix,i)]*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
            Szz[I2D(ix,i)] -= L2M[I2D(ix,i)]*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
        }
    }
    
    
    // Derivate Uz forward with respect to x
    der->ddx_fw(Uz1);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I2D(ix,iz)] = M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Vzx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] -= M[I2D(ix,iz)]*(Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Vzx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i,iz)];
            Sxz[I2D(i,iz)] -= M[I2D(i,iz)]*(Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Ux forward with respect to z
    der->ddz_fw(Ux1);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I2D(ix,iz)] += M[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Vxz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] -= M[I2D(ix,iz)]*(Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz)]);
            //Bottom
            i = iz + nz - lpml;
            Pml->Vxz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i)];
            Sxz[I2D(ix,i)] -= M[I2D(ix,i)]*(Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i)]);
        }
    }
}

template<typename T>
void WavesElastic2D_DS<T>::insertForcesource(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2;
    switch(sourcetype)
    {
        case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Ux2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Uz2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesElastic2D_DS<T>::insertPressuresource(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    int nt = this->getNt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2;
    switch(sourcetype)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Sxx[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            Szz[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesElastic2D_DS<T>::recordData(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    rs_field field = data->getField();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    dataarray = data->getData();
    int i,i1,i2;
    Index Im(this->getNx(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {

                    // Factor to make reciprocity work
                    factor = R[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)] - R[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)];
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                        }
                    }
                }
            }
	    break;
        case VX:
            Fielddata1 = this->getUx1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VZ:
            Fielddata1 = this->getUz1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        default:
            break;
    }

}

// Roll the displacement pointers
template<typename T>
void WavesElastic2D_DS<T>::roll()
{
    T *tmp;
    tmp = Ux2;
    Ux2 = Ux1;
    Ux1 = tmp;

    tmp = Uz2;
    Uz2 = Uz1;
    Uz1 = tmp;
}


template<typename T>
WavesElastic2D_DS<T>::~WavesElastic2D_DS() {
    /* Free allocated variables */
    free(Sxx);
    free(Szz);
    free(Sxz);
    free(Ux1);
    free(Ux2);
    free(Uz1);
    free(Uz2);
}

// =============== 3D ELASTIC DISPLACEMENT STRESS WAVES CLASS =============== //
/** The 3D Elastic WAVES model class
 *
 */

template<typename T>
WavesElastic3D_DS<T>::WavesElastic3D_DS(){
    // Nothing here
}

template<typename T>
WavesElastic3D_DS<T>::WavesElastic3D_DS(const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot): Waves<T>(3, _nx, _ny, _nz, _nt, _lpml, _dx, _dy, _dz, _dt, _ox, _oy, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ux1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ux2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uy1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uy2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uz1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uz2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesElastic3D_DS<T>::WavesElastic3D_DS(std::shared_ptr<rockseis::ModelElastic3D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ux1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ux2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uy1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uy2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uz1 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Uz2 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesElastic3D_DS<T>::~WavesElastic3D_DS() {
    /* Free allocated variables */
    free(Sxx);
    free(Syy);
    free(Szz);
    free(Sxz);
    free(Syz);
    free(Sxy);
    free(Ux1);
    free(Ux2);
    free(Uy1);
    free(Uy2);
    free(Uz1);
    free(Uz2);
}

template<typename T>
void WavesElastic3D_DS<T>::forwardstepDisplacement(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *Rx, *Ry, *Rz, *df;
    T dt;
    dt = this->getDt();
    Rx = model->getRx();
    Ry = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    
   //////////////////////////// VX //////////////////////////
   //
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ux2[I3D(ix,iy,iz)] = 2.0*Ux1[I3D(ix,iy,iz)] - Ux2[I3D(ix,iy,iz)] + dt*dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                
                Ux2[I3D(ix,iy,iz)] -= dt*dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Ux2[I3D(i,iy,iz)] -= dt*dt*Rx[I3D(i,iy,iz)]*(Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }

    // Derivate Sxy backward with respect to y
    der->ddy_bw(Sxy);
    // Compute Ux
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ux2[I3D(ix,iy,iz)] += dt*dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Sxyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                Ux2[I3D(ix,iy,iz)] -= dt*dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Sxyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Ux2[I3D(ix,i,iz)] -= dt*dt*Rx[I3D(ix,i,iz)]*(Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Sxz backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Ux
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ux2[I3D(ix,iy,iz)] += dt*dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Sxzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];
                
                Ux2[I3D(ix,iy,iz)] -= dt*dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Ux2[I3D(ix,iy,i)] -= dt*dt*Rx[I3D(ix,iy,i)]*(Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
 
   
   //////////////////////////// VY //////////////////////////
    

    // Derivate Sxy backward with respect to x
    der->ddx_bw(Sxy);
    // Compute Uy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uy2[I3D(ix,iy,iz)] = 2.0*Uy1[I3D(ix,iy,iz)] - Uy2[I3D(ix,iy,iz)] + dt*dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                
                Uy2[I3D(ix,iy,iz)] -= dt*dt*Ry[I3D(ix,iy,iz)]*(Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Uy2[I3D(i,iy,iz)] -= dt*dt*Ry[I3D(i,iy,iz)]*(Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }


    // Derivate Syy forward with respect to y
    der->ddy_fw(Syy);
    // Compute Uy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uy2[I3D(ix,iy,iz)] += dt*dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                Uy2[I3D(ix,iy,iz)] -= dt*dt*Ry[I3D(ix,iy,iz)]*(Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Syy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Uy2[I3D(ix,i,iz)] -= dt*dt*Ry[I3D(ix,i,iz)]*(Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Syz backward with respect to z
    der->ddz_bw(Syz);
    // Compute Uy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uy2[I3D(ix,iy,iz)] += dt*dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Syzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];
                
                Uy2[I3D(ix,iy,iz)] -= dt*dt*Ry[I3D(ix,iy,iz)]*(Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Uy2[I3D(ix,iy,i)] -= dt*dt*Ry[I3D(ix,iy,i)]*(Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
 
 
   //////////////////////////// VZ //////////////////////////
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Uz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uz2[I3D(ix,iy,iz)] = 2.0*Uz1[I3D(ix,iy,iz)] - Uz2[I3D(ix,iy,iz)] + dt*dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                Uz2[I3D(ix,iy,iz)] -= dt*dt*Rz[I3D(ix,iy,iz)]*(Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);

                // Right
                i = ix + nx - lpml;
                Pml->Sxzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Uz2[I3D(i,iy,iz)] -= dt*dt*Rz[I3D(i,iy,iz)]*(Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }

    // Derivate Syz backward with respect to y
    der->ddy_bw(Syz);
    // Compute Uz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uz2[I3D(ix,iy,iz)] += dt*dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                Uz2[I3D(ix,iy,iz)] -= dt*dt*Rz[I3D(ix,iy,iz)]*(Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Syzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Uz2[I3D(ix,i,iz)] -= dt*dt*Rz[I3D(ix,i,iz)]*(Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Uz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Uz2[I3D(ix,iy,iz)] += dt*dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Szz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                Uz2[I3D(ix,iy,iz)] -= dt*dt*Rz[I3D(ix,iy,iz)]*(Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Szz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Uz2[I3D(ix,iy,i)] -= dt*dt*Rz[I3D(ix,iy,i)]*(Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
    
    
} // End of forwardstepVelocity

template<typename T>
void WavesElastic3D_DS<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *L, *L2M, *M_xz, *M_yz, *M_xy, *df;
    L = model->getL();
    L2M = model->getL2M();
    M_xy = model->getM_xy();
    M_xz = model->getM_xz();
    M_yz = model->getM_yz();
    df = der->getDf();
    
    // Derivate Ux1 backward with respect to x
    der->ddx_bw(Ux1);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] = L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] = L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] = L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                
                Sxx[I3D(ix,iy,iz)] -= L2M[I3D(ix,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Vxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Sxx[I3D(i,iy,iz)] -= L2M[I3D(i,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Syy[I3D(i,iy,iz)] -= L[I3D(i,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Szz[I3D(i,iy,iz)] -= L[I3D(i,iy,iz)]*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
    
    // Derivate Uy1 backward with respect to y
    der->ddy_bw(Uy1);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] +=  L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] +=  L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] +=  L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                Sxx[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= L2M[I3D(ix,iy,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Vyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Sxx[I3D(ix,i,iz)] -= L[I3D(ix,i,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Syy[I3D(ix,i,iz)] -= L2M[I3D(ix,i,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Szz[I3D(ix,i,iz)] -= L[I3D(ix,i,iz)]*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
    
    // Derivate Uz1 backward with respect to z
    der->ddz_bw(Uz1);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I3D(ix,iy,iz)] +=  L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] +=  L[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] +=  L2M[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];
                
                Sxx[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= L[I3D(ix,iy,iz)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= L2M[I3D(ix,iy,iz)]*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Sxx[I3D(ix,iy,i)] -= L[I3D(ix,iy,i)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Syy[I3D(ix,iy,i)] -= L[I3D(ix,iy,i)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Szz[I3D(ix,iy,i)] -= L2M[I3D(ix,iy,i)]*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
    
   //////////////////////////// Sxz //////////////////////////
  
    // Derivate Uz1 forward with respect to x
    der->ddx_fw(Uz1);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxz[I3D(ix,iy,iz)] =  M_xz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                
                Sxz[I3D(ix,iy,iz)] -= M_xz[I3D(ix,iy,iz)]*(Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Vzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Sxz[I3D(i,iy,iz)] -= M_xz[I3D(i,iy,iz)]*(Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
 

    // Derivate Ux1 forward with respect to z
    der->ddz_fw(Ux1);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxz[I3D(ix,iy,iz)] +=  M_xz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vxz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                Sxz[I3D(ix,iy,iz)] -= M_xz[I3D(ix,iy,iz)]*(Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);

                //Bottom
                i = iz + nz - lpml;
                Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Sxz[I3D(ix,iy,i)] -= M_xz[I3D(ix,iy,i)]*(Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }


   //////////////////////////// Syz //////////////////////////
    // Derivate Uy forward with respect to z
    der->ddz_fw(Uy1);
    // Compute Syz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Syz[I3D(ix,iy,iz)] =  M_yz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vyz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                Syz[I3D(ix,iy,iz)] -= M_yz[I3D(ix,iy,iz)]*(Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);

                //Bottom
                i = iz + nz - lpml;
                Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Syz[I3D(ix,iy,i)] -= M_yz[I3D(ix,iy,i)]*(Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }

    // Derivate Uz1 forward with respect to y
    der->ddy_fw(Uz1);
    // Compute Syz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Syz[I3D(ix,iy,iz)] +=  M_yz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                Syz[I3D(ix,iy,iz)] -= M_yz[I3D(ix,iy,iz)]*(Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Vzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Syz[I3D(ix,i,iz)] -= M_yz[I3D(ix,i,iz)]*(Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }


   //////////////////////////// Sxy //////////////////////////
    // Derivate Ux1 forward with respect to y
    der->ddy_fw(Ux1);
    // Compute Sxy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxy[I3D(ix,iy,iz)] =  M_xy[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vxy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                Sxy[I3D(ix,iy,iz)] -= M_xy[I3D(ix,iy,iz)]*(Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);

                //Back
                i = iy + ny - lpml;
                Pml->Vxy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Sxy[I3D(ix,i,iz)] -= M_xy[I3D(ix,i,iz)]*(Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }

    // Derivate Uy1 forward with respect to x
    der->ddx_fw(Uy1);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxy[I3D(ix,iy,iz)] +=  M_xy[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                Sxy[I3D(ix,iy,iz)] -= M_xy[I3D(ix,iy,iz)]*(Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);

                // Right
                i = ix + nx - lpml;
                Pml->Vyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Sxy[I3D(i,iy,iz)] -= M_xy[I3D(i,iy,iz)]*(Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
}

template<typename T>
void WavesElastic3D_DS<T>::insertForcesource(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }
    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2,i3;
    switch(sourcetype)
    {
	    case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Ux2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        case VY:
            Mod = model->getRy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Uy2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Uz2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesElastic3D_DS<T>::insertPressuresource(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    int nt = this->getNt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }
    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2,i3;
    switch(sourcetype)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Sxx[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                                Syy[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                                Szz[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesElastic3D_DS<T>::recordData(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Fielddata3 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, ny, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i,i1,i2,i3;
    Index Im(this->getNx(), this->getNy(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSyy();
            Fielddata3 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {

                    // Factor to make reciprocity work
                    factor = (3.*R[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)] - 4.*R[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)])/3.;
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata3[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                            }
                        }
                    }
                }
            }
            break;
        case VX:
            Fielddata1 = this->getUx1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        case VY:
            Fielddata1 = this->getUy1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;

        case VZ:
            Fielddata1 = this->getUz1();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

// Roll the displacement pointers
template<typename T>
void WavesElastic3D_DS<T>::roll()
{
    T *tmp;
    tmp = Ux2;
    Ux2 = Ux1;
    Ux1 = tmp;

    tmp = Uy2;
    Uy2 = Uy1;
    Uy1 = tmp;

    tmp = Uz2;
    Uz2 = Uz1;
    Uz1 = tmp;
}

// =============== 2D VISCOELASTIC VELOCITY-STRESS WAVES CLASS =============== //
/** The 2D Viscoelastic WAVES model class
 *
 */

template<typename T>
WavesViscoelastic2D<T>::WavesViscoelastic2D()	///< Constructor
{
    // Do nothing
}


template<typename T>
WavesViscoelastic2D<T>::WavesViscoelastic2D(const int _nx, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot): Waves<T>(2, _nx, 1, _nz, _nt, _lpml, _dx, 1.0, _dz, _dt, _ox, 0.0, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mzz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesViscoelastic2D<T>::WavesViscoelastic2D(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);
   
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic2D<T>>(_nx, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mxx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mzz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Mxz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}


template<typename T>
void WavesViscoelastic2D<T>::forwardstepVelocity(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *Rx, *Rz, *df;
    Rx = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I2D(ix,iz)] += dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix,iz)];
            Vx[I2D(ix,iz)] -= dt*Rx[I2D(ix,iz)]*(Pml->Sxx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i,iz)];
            Vx[I2D(i,iz)] -= dt*Rx[I2D(i,iz)]*(Pml->Sxx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i,iz)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I2D(ix,iz)] += dt*Rx[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Sxzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz)];
            
            Vx[I2D(ix,iz)] -= dt*Rx[I2D(ix,iz)]*(Pml->Sxzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Sxzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i)];
            Vx[I2D(ix,i)] -= dt*Rx[I2D(ix,i)]*(Pml->Sxzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I2D(ix,iz)] += dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxzx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix,iz)];
            Vz[I2D(ix,iz)] -= dt*Rz[I2D(ix,iz)]*(Pml->Sxzx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxzx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i,iz)];
            Vz[I2D(i,iz)] -= dt*Rz[I2D(i,iz)]*(Pml->Sxzx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I2D(ix,iz)] += dt*Rz[I2D(ix,iz)]*df[I2D(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Szz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz)];
            
            Vz[I2D(ix,iz)] -= dt*Rz[I2D(ix,iz)]*(Pml->Szz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Szz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i)];
            Vz[I2D(ix,i)] -= dt*Rz[I2D(ix,i)]*(Pml->Szz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i)]);
        }
    }
    
} // End of forwardstepVelocity


template<typename T>
void WavesViscoelastic2D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *L2M, *M, *M_xz, *Tp, *Ts, *Ts_xz, *df;
    T dt = this->getDt();
    T to = dt*(2*PI*model->getF0());
    M = model->getM();
    L2M = model->getL2M();
    M_xz = model->getM_xz();
    Tp = model->getTp();
    Ts = model->getTs();
    Ts_xz = model->getTs_xz();
    df = der->getDf();
    T C11, C13, C44; 
    
    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            C11 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)];
            if(iz == lpml && model->getFs()){
                // Free surface conditions
                C13 = 0.0;
            }else{
                C13 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)] - 2*M[I2D(ix,iz)]*Ts[I2D(ix,iz)];
            }
            Sxx[I2D(ix,iz)] += dt*C11*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += dt*C13*df[I2D(ix,iz)];

            Mxx[I2D(ix,iz)] = ( (1-0.5*to)*Mxx[I2D(ix,iz)]- to*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1))*df[I2D(ix,iz)] ) ) / (1 + 0.5*to);
            Mzz[I2D(ix,iz)] = ( (1-0.5*to)*Mzz[I2D(ix,iz)]- to*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1)-2*M[I2D(ix,iz)]*(Ts[I2D(ix,iz)]-1))*df[I2D(ix,iz)] ) ) / (1 + 0.5*to);
        }
    }


    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            Szz[I2D(ix,iz)] = 0.0;
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            C11 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)];
            C13 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)] - 2*M[I2D(ix,iz)]*Ts[I2D(ix,iz)];
            Pml->Vxx_left[I2D_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->A_ltf[ix]*df[I2D(ix,iz)];
            
            Sxx[I2D(ix,iz)] -= dt*C11*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            Mxx[I2D(ix,iz)] += (to/(1.+0.5*to))*(L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1))*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);

            Szz[I2D(ix,iz)] -= dt*C13*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            Mzz[I2D(ix,iz)] += (to/(1.+0.5*to))*(L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1)-2*M[I2D(ix,iz)]*(Ts[I2D(ix,iz)]-1))*(Pml->Vxx_left[I2D_lr(ix,iz)] + Pml->C_ltf[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            C11 = L2M[I2D(i,iz)]*Tp[I2D(i,iz)];
            C13 = L2M[I2D(i,iz)]*Tp[I2D(i,iz)] - 2*M[I2D(i,iz)]*Ts[I2D(i,iz)];
            Pml->Vxx_right[I2D_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->A_rbb[ix]*df[I2D(i,iz)];

            Sxx[I2D(i,iz)] -= dt*C11*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
            Mxx[I2D(i,iz)] += (to/(1.+0.5*to))*(L2M[I2D(i,iz)]*(Tp[I2D(i,iz)]-1))*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);

            Szz[I2D(i,iz)] -= dt*C13*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
            Mzz[I2D(i,iz)] += (to/(1.+0.5*to))*(L2M[I2D(i,iz)]*(Tp[I2D(i,iz)]-1)-2*M[I2D(i,iz)]*(Ts[I2D(i,iz)]-1))*(Pml->Vxx_right[I2D_lr(ix,iz)] + Pml->C_rbb[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            C11 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)];
            if(iz == lpml && model->getFs()){
                // Free surface conditions
                C13 = 0.0;
            }else{
                C13 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)] - 2*M[I2D(ix,iz)]*Ts[I2D(ix,iz)];
            }
            Sxx[I2D(ix,iz)] += dt*C13*df[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += dt*C11*df[I2D(ix,iz)];

            Mxx[I2D(ix,iz)] -= (to*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1)-2*M[I2D(ix,iz)]*(Ts[I2D(ix,iz)]-1))*df[I2D(ix,iz)] ) ) / (1 + 0.5*to);
            Mzz[I2D(ix,iz)] -= (to*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1))*df[I2D(ix,iz)] ) ) / (1 + 0.5*to);

            // Convolution sum
            Sxx[I2D(ix,iz)] += dt*Mxx[I2D(ix,iz)];
            Szz[I2D(ix,iz)] += dt*Mzz[I2D(ix,iz)];
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            Szz[I2D(ix,iz)] = 0.0;
        }
    }

    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            C11 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)];
            C13 = L2M[I2D(ix,iz)]*Tp[I2D(ix,iz)] - 2*M[I2D(ix,iz)]*Ts[I2D(ix,iz)];
            // Top
            Pml->Vzz_top[I2D_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->A_ltf[iz]*df[I2D(ix,iz)];
            
            Sxx[I2D(ix,iz)] -= dt*C13*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            Mxx[I2D(ix,iz)] += (to/(1.+0.5*to))*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1)-2*M[I2D(ix,iz)]*(Ts[I2D(ix,iz)]-1)))*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            Szz[I2D(ix,iz)] -= dt*C11*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);
            Mzz[I2D(ix,iz)] += (to/(1.+0.5*to))*((L2M[I2D(ix,iz)]*(Tp[I2D(ix,iz)]-1)))*(Pml->Vzz_top[I2D_tb(ix,iz)] + Pml->C_ltf[iz]*df[I2D(ix,iz)]);

            //Bottom
            i = iz + nz - lpml;
            C11 = L2M[I2D(ix,i)]*Tp[I2D(ix,i)];
            C13 = L2M[I2D(ix,i)]*Tp[I2D(ix,i)] - 2*M[I2D(ix,i)]*Ts[I2D(ix,i)];
            Pml->Vzz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb[iz]*df[I2D(ix,i)];
            Sxx[I2D(ix,i)] -= dt*C13*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
            Mxx[I2D(ix,i)] += (to/(1.+0.5*to))*(L2M[I2D(ix,i)]*(Tp[I2D(ix,i)]-1)-2*M[I2D(ix,i)]*(Ts[I2D(ix,i)]-1))*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
            Szz[I2D(ix,i)] -= dt*C11*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
            Mzz[I2D(ix,i)] += (to/(1.+0.5*to))*(L2M[I2D(ix,i)]*(Tp[I2D(ix,i)]-1))*(Pml->Vzz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb[iz]*df[I2D(ix,i)]);
        }
    }
    
    
    // Derivate Vz forward with respect to x
    der->ddx_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            C44 = M_xz[I2D(ix,iz)]*Ts_xz[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] += dt*C44*df[I2D(ix,iz)];
            Mxz[I2D(ix,iz)] = ( (1-0.5*to)*Mxz[I2D(ix,iz)]- to*M_xz[I2D(ix,iz)]*(Ts_xz[I2D(ix,iz)]-1)*df[I2D(ix,iz)] ) / (1 + 0.5*to);
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            C44 = M_xz[I2D(ix,iz)]*Ts_xz[I2D(ix,iz)];
            Pml->Vzx_left[I2D_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] -= dt*C44*(Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix,iz)]);
            Mxz[I2D(ix,iz)] += (to/(1.+0.5*to))*M_xz[I2D(ix,iz)]*(Ts_xz[I2D(ix,iz)]-1)*(Pml->Vzx_left[I2D_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I2D(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            C44 = M_xz[I2D(i,iz)]*Ts_xz[I2D(i,iz)];
            Pml->Vzx_right[I2D_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I2D(i,iz)];
            Sxz[I2D(i,iz)] -= dt*C44*(Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i,iz)]);
            Mxz[I2D(i,iz)] += (to/(1.+0.5*to))*M_xz[I2D(i,iz)]*(Ts_xz[I2D(i,iz)]-1)*(Pml->Vzx_right[I2D_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I2D(i,iz)]);
        }
    }
    
    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            C44 = M_xz[I2D(ix,iz)]*Ts_xz[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] += dt*C44*df[I2D(ix,iz)];

            Mxz[I2D(ix,iz)] -= (to*M_xz[I2D(ix,iz)]*(Ts_xz[I2D(ix,iz)]-1)*df[I2D(ix,iz)] ) / (1 + 0.5*to);
            Sxz[I2D(ix,iz)] += dt*Mxz[I2D(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            C44 = M_xz[I2D(ix,iz)]*Ts_xz[I2D(ix,iz)];
            Pml->Vxz_top[I2D_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I2D(ix,iz)];
            Sxz[I2D(ix,iz)] -= dt*C44*(Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz)]);
            Mxz[I2D(ix,iz)] += (to/(1.+0.5*to))*M_xz[I2D(ix,iz)]*(Ts_xz[I2D(ix,iz)]-1)*(Pml->Vxz_top[I2D_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I2D(ix,iz)]);

            //Bottom
            i = iz + nz - lpml;
            C44 = M_xz[I2D(ix,i)]*Ts_xz[I2D(ix,i)];
            Pml->Vxz_bottom[I2D_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I2D(ix,i)];
            Sxz[I2D(ix,i)] -= dt*C44*(Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i)]);
            Mxz[I2D(ix,i)] += (to/(1.+0.5*to))*M_xz[I2D(ix,i)]*(Ts_xz[I2D(ix,i)]-1)*(Pml->Vxz_bottom[I2D_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I2D(ix,i)]);
        }
    }
}
    
template<typename T>
void WavesViscoelastic2D<T>::insertSource(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }

    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i, i1, i2;
    switch(sourcetype)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Sxx[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            Szz[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Vx[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            Vz[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)] += dt*Mod[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesViscoelastic2D<T>::recordData(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it){
    Point2D<int> *map;
    Point2D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    rs_field field = data->getField();

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    dataarray = data->getData();
    int i,i1,i2;
    Index Im(this->getNx(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {

                    // Factor to make reciprocity work
                    factor = R[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)]*Vp[Im(map[i].x, map[i].y)] - R[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)]*Vs[Im(map[i].x, map[i].y)];
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                            dataarray[Idat(it,i)] += (1./2.)*Fielddata2[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                        }
                    }
                }
            }
	    break;
        case VX:
            Fielddata1 = this->getVx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        case VZ:
            Fielddata1 = this->getVz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >= 0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            dataarray[Idat(it,i)] += Fielddata1[I2D(lpml + map[i].x - (LANC_SIZE-1) + i2, lpml + map[i].y - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i2), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i1) ,LANC_SIZE);
                        }
                    }
                }
            }
            break;
        default:
            break;
    }

}

template<typename T>
WavesViscoelastic2D<T>::~WavesViscoelastic2D() {
    /* Free allocated variables */
    free(Sxx);
    free(Szz);
    free(Sxz);
    free(Mxx);
    free(Mzz);
    free(Mxz);
    free(Vx);
    free(Vz);
}

// =============== 3D VISCOELASTIC WAVES CLASS =============== //
/** The 3D Viscoelastic WAVES model class
 *
 */

template<typename T>
WavesViscoelastic3D<T>::WavesViscoelastic3D(){
    // Nothing here
}

template<typename T>
WavesViscoelastic3D<T>::WavesViscoelastic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot): Waves<T>(3, _nx, _ny, _nz, _nt, _lpml, _dx, _dy, _dz, _dt, _ox, _oy, _oz, _ot) {
    
    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Myy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mzz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Myz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesViscoelastic3D<T>::WavesViscoelastic3D(std::shared_ptr<rockseis::ModelViscoelastic3D<T>> model, int _nt, T _dt, T _ot): Waves<T>() {
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _dim, _lpml;

    /* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    _dim = model->getDim();
    _lpml = model->getLpml();
    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
    this->setDim(_dim);

    /* Create associated PML class */
    Pml = std::make_shared<PmlElastic3D<T>>(_nx, _ny, _nz, _lpml, _dt);
    
    /* Allocate memory variables */
    int nx_pml, ny_pml, nz_pml;
    this->setDim(2);
    nx_pml = _nx + 2*_lpml;
    ny_pml = _ny + 2*_lpml;
    nz_pml = _nz + 2*_lpml;
    Sxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Szz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Syz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Sxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Myy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mzz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Myz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Mxy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Vz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
WavesViscoelastic3D<T>::~WavesViscoelastic3D() {
    /* Free allocated variables */
    free(Sxx);
    free(Syy);
    free(Szz);
    free(Sxz);
    free(Syz);
    free(Sxy);
    free(Mxx);
    free(Myy);
    free(Mzz);
    free(Mxz);
    free(Myz);
    free(Mxy);
    free(Vx);
    free(Vy);
    free(Vz);
}

template<typename T>
void WavesViscoelastic3D<T>::forwardstepVelocity(std::shared_ptr<rockseis::ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *Rx, *Ry, *Rz, *df;
    T dt;
    dt = this->getDt();
    Rx = model->getRx();
    Ry = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    
   //////////////////////////// VX //////////////////////////
   //
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                
                Vx[I3D(ix,iy,iz)] -= dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Vx[I3D(i,iy,iz)] -= dt*Rx[I3D(i,iy,iz)]*(Pml->Sxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }

    // Derivate Sxy backward with respect to y
    der->ddy_bw(Sxy);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Sxyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                
                Vx[I3D(ix,iy,iz)] -= dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                i = iy + ny - lpml;
                //Back
                Pml->Sxyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Vx[I3D(ix,i,iz)] -= dt*Rx[I3D(ix,i,iz)]*(Pml->Sxyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Sxy backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I3D(ix,iy,iz)] += dt*Rx[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Sxzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];
                
                Vx[I3D(ix,iy,iz)] -= dt*Rx[I3D(ix,iy,iz)]*(Pml->Sxzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Vx[I3D(ix,iy,i)] -= dt*Rx[I3D(ix,iy,i)]*(Pml->Sxzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
 
   
   //////////////////////////// VY //////////////////////////
    

    // Derivate Sxy backward with respect to x
    der->ddx_bw(Sxy);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                
                Vy[I3D(ix,iy,iz)] -= dt*Ry[I3D(ix,iy,iz)]*(Pml->Sxyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Vy[I3D(i,iy,iz)] -= dt*Ry[I3D(i,iy,iz)]*(Pml->Sxyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }


    // Derivate Syy forward with respect to y
    der->ddy_fw(Syy);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                
                Vy[I3D(ix,iy,iz)] -= dt*Ry[I3D(ix,iy,iz)]*(Pml->Syy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Syy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Vy[I3D(ix,i,iz)] -= dt*Ry[I3D(ix,i,iz)]*(Pml->Syy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Syz backward with respect to z
    der->ddz_bw(Syz);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I3D(ix,iy,iz)] += dt*Ry[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Syzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];
                
                Vy[I3D(ix,iy,iz)] -= dt*Ry[I3D(ix,iy,iz)]*(Pml->Syzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Vy[I3D(ix,iy,i)] -= dt*Ry[I3D(ix,iy,i)]*(Pml->Syzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
 
 
   //////////////////////////// VZ //////////////////////////
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                
                Vz[I3D(ix,iy,iz)] -= dt*Rz[I3D(ix,iy,iz)]*(Pml->Sxzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Vz[I3D(i,iy,iz)] -= dt*Rz[I3D(i,iy,iz)]*(Pml->Sxzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }

    // Derivate Syz backward with respect to y
    der->ddy_bw(Syz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                
                Vz[I3D(ix,iy,iz)] -= dt*Rz[I3D(ix,iy,iz)]*(Pml->Syzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Syzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Vz[I3D(ix,i,iz)] -= dt*Rz[I3D(ix,i,iz)]*(Pml->Syzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
 
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I3D(ix,iy,iz)] += dt*Rz[I3D(ix,iy,iz)]*df[I3D(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Szz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                Vz[I3D(ix,iy,iz)] -= dt*Rz[I3D(ix,iy,iz)]*(Pml->Szz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Szz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Vz[I3D(ix,iy,i)] -= dt*Rz[I3D(ix,iy,i)]*(Pml->Szz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
    
    
} // End of forwardstepVelocity

template<typename T>
void WavesViscoelastic3D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T to = dt*(2*PI*model->getF0());
    T *L2M, *M, *M_xz, *M_yz, *M_xy, *Tp, *Ts, *Ts_xz, *Ts_yz, *Ts_xy, *df;

    L2M = model->getL2M();
    M = model->getM();
    M_xz = model->getM_xz();
    M_xy = model->getM_xy();
    M_yz = model->getM_yz();

    Tp = model->getTp();
    Ts = model->getTs();
    Ts_xz = model->getTs_xz();
    Ts_xy = model->getTs_xy();
    Ts_yz = model->getTs_yz();
    df = der->getDf();
    T C11, C13, C44, T11, T13, T44;
    
    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                if(iz == lpml && model->getFs()){
                    // Free surface conditions
                    C13 = 0.0;
                }else{
                    C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                }
                Sxx[I3D(ix,iy,iz)] += dt*C11*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];

                Mxx[I3D(ix,iy,iz)] = ( (1-0.5*to)*Mxx[I3D(ix,iy,iz)]- to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Myy[I3D(ix,iy,iz)] = ( (1-0.5*to)*Myy[I3D(ix,iy,iz)]- to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Mzz[I3D(ix,iy,iz)] = ( (1-0.5*to)*Mzz[I3D(ix,iy,iz)]- to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }

    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1));
                Pml->Vxx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I3D(ix,iy,iz)];
                
                Sxx[I3D(ix,iy,iz)] -= dt*C11*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Mxx[I3D(ix,iy,iz)] += T11*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Myy[I3D(ix,iy,iz)] += T13*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);
                Mzz[I3D(ix,iy,iz)] += T13*(Pml->Vxx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I3D(ix,iy,iz)]);

                // Right
                i = ix + nx - lpml;
                C11 = L2M[I3D(i,iy,iz)]*Tp[I3D(i,iy,iz)];
                C13 = L2M[I3D(i,iy,iz)]*Tp[I3D(i,iy,iz)] - 2*M[I3D(i,iy,iz)]*Ts[I3D(i,iy,iz)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(i,iy,iz)]*(Tp[I3D(i,iy,iz)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(i,iy,iz)]*(Tp[I3D(i,iy,iz)]-1)-2*M[I3D(i,iy,iz)]*(Ts[I3D(i,iy,iz)]-1));
                Pml->Vxx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I3D(i,iy,iz)];
                Sxx[I3D(i,iy,iz)] -= dt*C11*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Syy[I3D(i,iy,iz)] -= dt*C13*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Szz[I3D(i,iy,iz)] -= dt*C13*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Mxx[I3D(i,iy,iz)] += T11*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Myy[I3D(i,iy,iz)] += T13*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
                Mzz[I3D(i,iy,iz)] += T13*(Pml->Vxx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
    

    // Derivate Vy backward with respect to y
    der->ddy_bw(Vy);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                if(iz == lpml && model->getFs()){
                    // Free surface conditions
                    C13 = 0.0;
                }else{
                    C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                }
                Sxx[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] += dt*C11*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];
                Mxx[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Myy[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Mzz[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
            }
        }
    }
    
    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1));
                Pml->Vyy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I3D(ix,iy,iz)];
                
                Sxx[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= dt*C11*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Mxx[I3D(ix,iy,iz)] += T13*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Myy[I3D(ix,iy,iz)] += T11*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                Mzz[I3D(ix,iy,iz)] += T13*(Pml->Vyy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I3D(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                C11 = L2M[I3D(ix,i,iz)]*Tp[I3D(ix,i,iz)];
                C13 = L2M[I3D(ix,i,iz)]*Tp[I3D(ix,i,iz)] - 2*M[I3D(ix,i,iz)]*Ts[I3D(ix,i,iz)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(ix,i,iz)]*(Tp[I3D(ix,i,iz)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(ix,i,iz)]*(Tp[I3D(ix,i,iz)]-1)-2*M[I3D(ix,i,iz)]*(Ts[I3D(ix,i,iz)]-1));
                Pml->Vyy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I3D(ix,i,iz)];
                Sxx[I3D(ix,i,iz)] -= dt*C13*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Syy[I3D(ix,i,iz)] -= dt*C11*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Szz[I3D(ix,i,iz)] -= dt*C13*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Mxx[I3D(ix,i,iz)] += T13*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Myy[I3D(ix,i,iz)] += T11*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
                Mzz[I3D(ix,i,iz)] += T13*(Pml->Vyy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }
    
    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                if(iz == lpml && model->getFs()){
                    // Free surface conditions
                    C13 = 0.0;
                }else{
                    C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                }
                Sxx[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] += dt*C13*df[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] += dt*C11*df[I3D(ix,iy,iz)];

                Mxx[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Myy[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);
                Mzz[I3D(ix,iy,iz)] -= (to*((L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1))*df[I3D(ix,iy,iz)] ) ) / (1 + 0.5*to);


                // Convolution sum
                Sxx[I3D(ix,iy,iz)] += dt*Mxx[I3D(ix,iy,iz)];
                Syy[I3D(ix,iy,iz)] += dt*Myy[I3D(ix,iy,iz)];
                Szz[I3D(ix,iy,iz)] += dt*Mzz[I3D(ix,iy,iz)];
            }
        }
    }

    // Free surface conditions
    if(model->getFs()){
        iz = lpml;
        for(ix=0; ix < nx; ix++){
            for(iy=0; iy < ny; iy++){
                Szz[I3D(ix,iy,iz)] = 0.0;
            }
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                C11 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)];
                C13 = L2M[I3D(ix,iy,iz)]*Tp[I3D(ix,iy,iz)] - 2*M[I3D(ix,iy,iz)]*Ts[I3D(ix,iy,iz)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,iz)]*(Tp[I3D(ix,iy,iz)]-1)-2*M[I3D(ix,iy,iz)]*(Ts[I3D(ix,iy,iz)]-1));
                Pml->Vzz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I3D(ix,iy,iz)];

                Sxx[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Syy[I3D(ix,iy,iz)] -= dt*C13*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Szz[I3D(ix,iy,iz)] -= dt*C11*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Mxx[I3D(ix,iy,iz)] += T13*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Myy[I3D(ix,iy,iz)] += T13*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                Mzz[I3D(ix,iy,iz)] += T11*(Pml->Vzz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I3D(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                C11 = L2M[I3D(ix,iy,i)]*Tp[I3D(ix,iy,i)];
                C13 = L2M[I3D(ix,iy,i)]*Tp[I3D(ix,iy,i)] - 2*M[I3D(ix,iy,i)]*Ts[I3D(ix,iy,i)];
                T11 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,i)]*(Tp[I3D(ix,iy,i)]-1));
                T13 = (to/(1.+0.5*to))*(L2M[I3D(ix,iy,i)]*(Tp[I3D(ix,iy,i)]-1)-2*M[I3D(ix,iy,i)]*(Ts[I3D(ix,iy,i)]-1));
                Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I3D(ix,iy,i)];
                Sxx[I3D(ix,iy,i)] -= dt*C13*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Syy[I3D(ix,iy,i)] -= dt*C13*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Szz[I3D(ix,iy,i)] -= dt*C11*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Mxx[I3D(ix,iy,i)] += T13*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Myy[I3D(ix,iy,i)] += T13*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
                Mzz[I3D(ix,iy,i)] += T11*(Pml->Vzz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }
    
   //////////////////////////// Sxz //////////////////////////
  
    // Derivate Vz forward with respect to x
    der->ddx_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_xz[I3D(ix,iy,iz)]*Ts_xz[I3D(ix,iy,iz)];
                Sxz[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Mxz[I3D(ix,iy,iz)] = ( (1-0.5*to)*Mxz[I3D(ix,iy,iz)]- to*M_xz[I3D(ix,iy,iz)]*(Ts_xz[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                C44 = M_xz[I3D(ix,iy,iz)]*Ts_xz[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xz[I3D(ix,iy,iz)]*(Ts_xz[I3D(ix,iy,iz)]-1);
                Pml->Vzx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                
                Sxz[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                Mxz[I3D(ix,iy,iz)] += T44*(Pml->Vzx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                C44 = M_xz[I3D(i,iy,iz)]*Ts_xz[I3D(i,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xz[I3D(i,iy,iz)]*(Ts_xz[I3D(i,iy,iz)]-1);
                Pml->Vzx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Sxz[I3D(i,iy,iz)] -= dt*C44*(Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
                Mxz[I3D(i,iy,iz)] += T44*(Pml->Vzx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
 

    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_xz[I3D(ix,iy,iz)]*Ts_xz[I3D(ix,iy,iz)];
                Sxz[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Mxz[I3D(ix,iy,iz)] -= (to*M_xz[I3D(ix,iy,iz)]*(Ts_xz[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);

                // Convolution sum
                Sxz[I3D(ix,iy,iz)] += dt*Mxz[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                C44 = M_xz[I3D(ix,iy,iz)]*Ts_xz[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xz[I3D(ix,iy,iz)]*(Ts_xz[I3D(ix,iy,iz)]-1);
                Pml->Vxz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                
                Sxz[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                Mxz[I3D(ix,iy,iz)] += T44*(Pml->Vxz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                C44 = M_xz[I3D(ix,iy,i)]*Ts_xz[I3D(ix,iy,i)];
                T44 = (to/(1.+0.5*to))*M_xz[I3D(ix,iy,i)]*(Ts_xz[I3D(ix,iy,i)]-1);
                Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Sxz[I3D(ix,iy,i)] -= dt*C44*(Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
                Mxz[I3D(ix,iy,i)] += T44*(Pml->Vxz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }


   //////////////////////////// Syz //////////////////////////
    // Derivate Vy forward with respect to z
    der->ddz_fw(Vy);
    // Compute Syz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_yz[I3D(ix,iy,iz)]*Ts_yz[I3D(ix,iy,iz)];
                Syz[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Myz[I3D(ix,iy,iz)] = ( (1-0.5*to)*Myz[I3D(ix,iy,iz)]- to*M_yz[I3D(ix,iy,iz)]*(Ts_yz[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                C44 = M_yz[I3D(ix,iy,iz)]*Ts_yz[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_yz[I3D(ix,iy,iz)]*(Ts_yz[I3D(ix,iy,iz)]-1);
                Pml->Vyz_top[I3D_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I3D(ix,iy,iz)];
                
                Syz[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                Myz[I3D(ix,iy,iz)] += T44*(Pml->Vyz_top[I3D_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I3D(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                C44 = M_yz[I3D(ix,iy,i)]*Ts_yz[I3D(ix,iy,i)];
                T44 = (to/(1.+0.5*to))*M_yz[I3D(ix,iy,i)]*(Ts_yz[I3D(ix,iy,i)]-1);
                Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I3D(ix,iy,i)];
                Syz[I3D(ix,iy,i)] -= dt*C44*(Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
                Myz[I3D(ix,iy,i)] += T44*(Pml->Vyz_bottom[I3D_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I3D(ix,iy,i)]);
            }
        }
    }

    // Derivate Vz forward with respect to y
    der->ddy_fw(Vz);
    // Compute Syz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_yz[I3D(ix,iy,iz)]*Ts_yz[I3D(ix,iy,iz)];
                Syz[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Myz[I3D(ix,iy,iz)] -= (to*M_yz[I3D(ix,iy,iz)]*(Ts_yz[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);

                // Convolution sum
                Syz[I3D(ix,iy,iz)] += dt*Myz[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                C44 = M_yz[I3D(ix,iy,iz)]*Ts_yz[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_yz[I3D(ix,iy,iz)]*(Ts_yz[I3D(ix,iy,iz)]-1);
                Pml->Vzy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                
                Syz[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);
                Myz[I3D(ix,iy,iz)] += T44*(Pml->Vzy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                C44 = M_yz[I3D(ix,i,iz)]*Ts_yz[I3D(ix,i,iz)];
                T44 = (to/(1.+0.5*to))*M_yz[I3D(ix,i,iz)]*(Ts_yz[I3D(ix,i,iz)]-1);
                Pml->Vzy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Syz[I3D(ix,i,iz)] -= dt*C44*(Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
                Myz[I3D(ix,i,iz)] += T44*(Pml->Vzy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }


   //////////////////////////// Sxy //////////////////////////
    // Derivate Vx forward with respect to y
    der->ddy_fw(Vx);
    // Compute Sxy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_xy[I3D(ix,iy,iz)]*Ts_xy[I3D(ix,iy,iz)];
                Sxy[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Mxy[I3D(ix,iy,iz)] = ( (1-0.5*to)*Mxy[I3D(ix,iy,iz)]- to*M_xy[I3D(ix,iy,iz)]*(Ts_xy[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                C44 = M_xy[I3D(ix,iy,iz)]*Ts_xy[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xy[I3D(ix,iy,iz)]*(Ts_xy[I3D(ix,iy,iz)]-1);
                Pml->Vxy_front[I3D_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I3D(ix,iy,iz)];
                
                Sxy[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);
                Mxy[I3D(ix,iy,iz)] += T44*(Pml->Vxy_front[I3D_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I3D(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                C44 = M_xy[I3D(ix,i,iz)]*Ts_xy[I3D(ix,i,iz)];
                T44 = (to/(1.+0.5*to))*M_xy[I3D(ix,i,iz)]*(Ts_xy[I3D(ix,i,iz)]-1);
                Pml->Vxy_back[I3D_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I3D(ix,i,iz)];
                Sxy[I3D(ix,i,iz)] -= dt*C44*(Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
                Mxy[I3D(ix,i,iz)] += T44*(Pml->Vxy_back[I3D_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I3D(ix,i,iz)]);
            }
        }
    }

    // Derivate Vy forward with respect to x
    der->ddx_fw(Vy);
    // Compute Sxy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                C44 = M_xy[I3D(ix,iy,iz)]*Ts_xy[I3D(ix,iy,iz)];
                Sxy[I3D(ix,iy,iz)] += dt*C44*df[I3D(ix,iy,iz)];
                Mxy[I3D(ix,iy,iz)] -= (to*M_xy[I3D(ix,iy,iz)]*(Ts_xy[I3D(ix,iy,iz)]-1)*df[I3D(ix,iy,iz)] ) / (1 + 0.5*to);

                // Convolution sum
                Sxy[I3D(ix,iy,iz)] += dt*Mxy[I3D(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                C44 = M_xy[I3D(ix,iy,iz)]*Ts_xy[I3D(ix,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xy[I3D(ix,iy,iz)]*(Ts_xy[I3D(ix,iy,iz)]-1);
                Pml->Vyx_left[I3D_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I3D(ix,iy,iz)];
                
                Sxy[I3D(ix,iy,iz)] -= dt*C44*(Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                Mxy[I3D(ix,iy,iz)] += T44*(Pml->Vyx_left[I3D_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I3D(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                C44 = M_xy[I3D(i,iy,iz)]*Ts_xy[I3D(i,iy,iz)];
                T44 = (to/(1.+0.5*to))*M_xy[I3D(i,iy,iz)]*(Ts_xy[I3D(i,iy,iz)]-1);
                Pml->Vyx_right[I3D_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I3D(i,iy,iz)];
                Sxy[I3D(i,iy,iz)] -= dt*C44*(Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
                Mxy[I3D(i,iy,iz)] += T44*(Pml->Vyx_right[I3D_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I3D(i,iy,iz)]);
            }
        }
    }
}

template<typename T>
void WavesViscoelastic3D<T>::insertSource(std::shared_ptr<rockseis::ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        map = (source->getGeom())->getSmap();
        shift = (source->getGeom())->getSshift();
    }else{
        map = (source->getGeom())->getGmap();
        shift = (source->getGeom())->getGshift();
    }
    rs_field sourcetype = source->getField();
    wav = source->getData();
    int i,i1,i2,i3;
    switch(sourcetype)
    {
        case PRESSURE:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Sxx[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                                Syy[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                                Szz[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] -= dt*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        case VX:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Vx[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        case VY:
            Mod = model->getRy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Vy[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        case VZ:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                Vz[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)] += dt*Mod[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*wav[Idat(it,i)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE); 
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}

template<typename T>
void WavesViscoelastic3D<T>::recordData(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it){
    Point3D<int> *map;
    Point3D<T> *shift;
    T *dataarray; 
    T *Fielddata1 = NULL;
    T *Fielddata2 = NULL;
    T *Fielddata3 = NULL;
    T *Vp, *Vs, *R;
    T factor;
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, ny, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;

    // Get correct map (data or receiver mapping)
    if(maptype == SMAP) {
        map = (data->getGeom())->getSmap();
        shift = (data->getGeom())->getSshift();
    }else{
        map = (data->getGeom())->getGmap();
        shift = (data->getGeom())->getGshift();
    }

    rs_field field = data->getField();
    dataarray = data->getData();
    int i,i1,i2,i3;
    Index Im(this->getNx(), this->getNy(), this->getNz());
    switch(field)
    {
        case PRESSURE:
            Fielddata1 = this->getSxx();
            Fielddata2 = this->getSyy();
            Fielddata3 = this->getSzz();
            Vp = model->getVp();
            Vs = model->getVs();
            R = model->getR();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {

                    // Factor to make reciprocity work
                    factor = (3.*R[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)]*Vp[Im(map[i].x, map[i].y, map[i].z)] - 4.*R[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)]*Vs[Im(map[i].x, map[i].y, map[i].z)])/3.;
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata2[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                                dataarray[Idat(it,i)] += (1./3.)*Fielddata3[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE)/factor;
                            }
                        }
                    }
                }
            }
            break;
        case VX:
            Fielddata1 = this->getVx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        case VY:
            Fielddata1 = this->getVy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;

        case VZ:
            Fielddata1 = this->getVz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    for(i1=0; i1<2*LANC_SIZE; i1++){
                        for(i2=0; i2<2*LANC_SIZE; i2++){
                            for(i3=0; i3<2*LANC_SIZE; i3++){
                                dataarray[Idat(it,i)] += Fielddata1[I3D(lpml + map[i].x - (LANC_SIZE-1) + i3, lpml + map[i].y - (LANC_SIZE-1) + i2, lpml + map[i].z - (LANC_SIZE-1) + i1)]*LANC(shift[i].x + (LANC_SIZE-1-i3), LANC_SIZE)*LANC(shift[i].y + (LANC_SIZE-1-i2) ,LANC_SIZE)*LANC(shift[i].z + (LANC_SIZE-1-i1) ,LANC_SIZE);
                            }
                        }
                    }
                }
            }
            break;
        default:
            break;
    }
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Waves<float>;
template class Waves<double>;

template class WavesAcoustic2D<float>;
template class WavesAcoustic2D<double>;
template class WavesAcoustic3D<float>;
template class WavesAcoustic3D<double>;
template class WavesElastic2D<float>;
template class WavesElastic2D<double>;
template class WavesElastic3D<float>;
template class WavesElastic3D<double>;

template class WavesElastic2D_DS<float>;
template class WavesElastic2D_DS<double>;
template class WavesElastic3D_DS<float>;
template class WavesElastic3D_DS<double>;

template class WavesViscoelastic2D<float>;
template class WavesViscoelastic2D<double>;
template class WavesViscoelastic3D<float>;
template class WavesViscoelastic3D<double>;
}


