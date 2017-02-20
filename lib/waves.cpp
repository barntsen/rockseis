// Include statements
#include "waves.h"

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
}

// =============== 2D ACOUSTIC MODEL CLASS =============== //
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
WavesAcoustic2D<T>::~WavesAcoustic2D() {
    /* Free allocated variables */
    free(P1);
    free(P2);
    free(Ax);
    free(Az);
}

template<typename T>
void WavesAcoustic2D<T>::forwardstepAcceleration(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *Rx, *Rz, *df;
    Rx = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    rockseis::Index I(nx,nz);
    rockseis::Index I_lr(lpml,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,lpml); // index for top and bottom pml zone
    
    // Derivate P forward with respect to x
    der->ddx_fw(P1);
    // Compute Ax
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Ax[I(ix,iz)] = Rx[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->P_left[I_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->P_left[I_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iz)];
            
            Ax[I(ix,iz)] -= Rx[I(ix,iz)]*(Pml->P_left[I_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->P_right[I_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->P_right[I_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iz)];
            Ax[I(i,iz)] -= Rx[I(i,iz)]*(Pml->P_right[I_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iz)]);
        }
    }
    
    
    // Derivate P forward with respect to z
    der->ddz_fw(P1);
    // Compute Az
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Az[I(ix,iz)] = Rz[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->P_top[I_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->P_top[I_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iz)];
            
            Az[I(ix,iz)] -= Rz[I(ix,iz)]*(Pml->P_top[I_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->P_bottom[I_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->P_bottom[I_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,i)];
            Az[I(ix,i)] -= Rz[I(ix,i)]*(Pml->P_bottom[I_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,i)]);
        }
    }
    
    
} // End of forwardstepAcceleration

template<typename T>
void WavesAcoustic2D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *L, *df;
    L = model->getL();
    df = der->getDf();
    rockseis::Index I(nx,nz);
    rockseis::Index I_lr(lpml,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,lpml); // index for top and bottom pml zone
    
    // Derivate Ax backward with respect to x
    der->ddx_bw(Ax);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            P2[I(ix,iz)] = 2.0 * P1[I(ix,iz)] - P2[I(ix,iz)] + dt*dt*L[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Axx_left[I_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Axx_left[I_lr(ix,iz)] + Pml->A_ltf[ix]*df[I(ix,iz)];
            
            P2[I(ix,iz)] -= dt*dt*L[I(ix,iz)]*(Pml->Axx_left[I_lr(ix,iz)] + Pml->C_ltf[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Axx_right[I_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Axx_right[I_lr(ix,iz)] + Pml->A_rbb[ix]*df[I(i,iz)];
            P2[I(i,iz)] -= dt*dt*L[I(i,iz)]*(Pml->Axx_right[I_lr(ix,iz)] + Pml->C_rbb[ix]*df[I(i,iz)]);
        }
    }
    
    
    // Derivate Az backward with respect to z
    der->ddz_bw(Az);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            P2[I(ix,iz)] +=  dt*dt*L[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Azz_top[I_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Azz_top[I_tb(ix,iz)] + Pml->A_ltf[iz]*df[I(ix,iz)];
            
            P2[I(ix,iz)] -= dt*dt*L[I(ix,iz)]*(Pml->Azz_top[I_tb(ix,iz)] + Pml->C_ltf[iz]*df[I(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Azz_bottom[I_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Azz_bottom[I_tb(ix,iz)] + Pml->A_rbb[iz]*df[I(ix,i)];
            P2[I(ix,i)] -= dt*dt*L[I(ix,i)]*(Pml->Azz_bottom[I_tb(ix,iz)] + Pml->C_rbb[iz]*df[I(ix,i)]);
        }
    }
    
}

template<typename T>
void WavesAcoustic2D<T>::insertSource(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, int sourcetype, int maptype, int it){
    Point2D<int> *map;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    nz = this->getNz() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == 0) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    wav = source->getData();
    int i;
    //Indexes 
    Index I(nx, nz); //Model and Field indexes
    Index Idat(nt, ntrace); // Data indexes
    switch(sourcetype)
    {
        case 0:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    Ax[I(lpml + map[i].x, lpml + map[i].y)] += Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 1:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    Az[I(lpml + map[i].x, lpml + map[i].y)] += Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 2:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            {
                if(map[i].x >= 0 && map[i].y >=0)
                { 
                    P2[I(lpml + map[i].x, lpml + map[i].y)] += dt*dt*Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        default:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    P2[I(lpml + map[i].x, lpml + map[i].y)] += dt*dt*Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
    }
}

template<typename T>
void WavesAcoustic2D<T>::recordData(std::shared_ptr<rockseis::Data2D<T>> data, int field, int maptype, int it){
    Point2D<int> *map;
    T *dataarray; 
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    nz = this->getNz() + 2*lpml;

    // Get correct map (data or receiver mapping)
    if(maptype == 0) {
        map = (data->getGeom())->getSmap();
    }else{
        map = (data->getGeom())->getGmap();
    }

    dataarray = data->getData();
    int i;
    Index I(nx, nz);
    Index Idat(nt, ntrace);
    switch(field)
    {
        case 0:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = P2[I(lpml + map[i].x, lpml + map[i].y)];
                }
            }
        case 1:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = Ax[I(lpml + map[i].x, lpml + map[i].y)];
                }
            }
            break;
        case 3:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = Az[I(lpml + map[i].x, lpml + map[i].y)];
                }
            }
            break;
        default:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0)
                {
                    dataarray[Idat(it,i)] = P2[I(lpml + map[i].x, lpml + map[i].y)];
                }
            }
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
void WavesAcoustic3D<T>::forwardstepAcceleration(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    T *Rx, *Ry, *Rz, *df;
    Rx = model->getRx();
    Ry = model->getRx();
    Rz = model->getRz();
    df = der->getDf();
    rockseis::Index I(nx,ny,nz);
    rockseis::Index I_lr(lpml,ny,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,ny,lpml); // index for top and bottom pml zone
    rockseis::Index I_fb(nx,lpml,nz); // index for front and back pml zone
    
    // Derivate P forward with respect to x
    der->ddx_fw(P1);
    // Compute Ax
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ax[I(ix,iy,iz)] = Rx[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->P_left[I_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->P_left[I_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iy,iz)];
                
                Ax[I(ix,iy,iz)] -= Rx[I(ix,iy,iz)]*(Pml->P_left[I_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->P_right[I_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->P_right[I_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iy,iz)];
                Ax[I(i,iy,iz)] -= Rx[I(i,iy,iz)]*(Pml->P_right[I_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iy,iz)]);
            }
        }
    }
    
    // Derivate P forward with respect to y
    der->ddy_fw(P1);
    // Compute Ay
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Ay[I(ix,iy,iz)] = Ry[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->P_front[I_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->P_front[I_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I(ix,iy,iz)];
                
                Ay[I(ix,iy,iz)] -= Ry[I(ix,iy,iz)]*(Pml->P_front[I_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I(ix,iy,iz)]);
                i = iy + ny - lpml;
                //Back
                Pml->P_back[I_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->P_back[I_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I(ix,i,iz)];
                Ay[I(ix,i,iz)] -= Ry[I(ix,i,iz)]*(Pml->P_back[I_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I(ix,i,iz)]);
            }
        }
    }
    
    // Derivate P forward with respect to z
    der->ddz_fw(P1);
    // Compute Az
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Az[I(ix,iy,iz)] = Rz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->P_top[I_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->P_top[I_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iy,iz)];
                
                Az[I(ix,iy,iz)] -= Rz[I(ix,iy,iz)]*(Pml->P_top[I_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->P_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->P_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,iy,i)];
                Az[I(ix,iy,i)] -= Rz[I(ix,iy,i)]*(Pml->P_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,iy,i)]);
            }
        }
    }
    
    
} // End of forwardstepAcceleration

template<typename T>
void WavesAcoustic3D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *L, *df;
    L = model->getL();
    df = der->getDf();
    rockseis::Index I(nx,ny,nz);
    rockseis::Index I_lr(lpml,ny,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,ny,lpml); // index for top and bottom pml zone
    rockseis::Index I_fb(nx,lpml,nz); // index for front and back pml zone
    
    // Derivate Ax backward with respect to x
    der->ddx_bw(Ax);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I(ix,iy,iz)] = 2.0 * P1[I(ix,iy,iz)] - P2[I(ix,iy,iz)] + dt*dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Axx_left[I_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Axx_left[I_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I(ix,iy,iz)];
                
                P2[I(ix,iy,iz)] -= dt*dt*L[I(ix,iy,iz)]*(Pml->Axx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Axx_right[I_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Axx_right[I_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I(i,iy,iz)];
                P2[I(i,iy,iz)] -= dt*dt*L[I(i,iy,iz)]*(Pml->Axx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
            }
        }
    }
    
    // Derivate Ay backward with respect to y
    der->ddy_bw(Ay);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I(ix,iy,iz)] +=  dt*dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Ayy_front[I_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Ayy_front[I_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I(ix,iy,iz)];
                
                P2[I(ix,iy,iz)] -= dt*dt*L[I(ix,iy,iz)]*(Pml->Ayy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                i = iy + ny - lpml;
                //Back
                Pml->Ayy_back[I_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Ayy_back[I_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I(ix,i,iz)];
                P2[I(ix,i,iz)] -= dt*dt*L[I(ix,i,iz)]*(Pml->Ayy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
            }
        }
    }
    
    
    // Derivate Az backward with respect to z
    der->ddz_bw(Az);
    // Compute P
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                P2[I(ix,iy,iz)] +=  dt*dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Azz_top[I_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Azz_top[I_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I(ix,iy,iz)];
                
                P2[I(ix,iy,iz)] -= dt*dt*L[I(ix,iy,iz)]*(Pml->Azz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Azz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Azz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I(ix,iy,i)];
                P2[I(ix,iy,i)] -= dt*dt*L[I(ix,iy,i)]*(Pml->Azz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
            }
        }
    }
    
}

template<typename T>
void WavesAcoustic3D<T>::insertSource(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, int sourcetype, int maptype, int it){
    Point3D<int> *map;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    nz = this->getNz() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == 0 ) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }
    wav = source->getData();
    int i;
    Index I(nx, ny, nz); //Model and Field indexes
    Index Idat(nt, ntrace); // Data indexes
    switch(sourcetype)
    {
        case 0:
            Mod = model->getL();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    P2[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*dt*Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
	    case 1:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Ax[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 2:
            Mod = model->getRy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Ay[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 3:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Az[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
    }
}

template<typename T>
void WavesAcoustic3D<T>::recordData(std::shared_ptr<rockseis::Data3D<T>> data, int field, int maptype, int it){
    Point3D<int> *map;
    T *dataarray; 
    int ntrace = data->getNtrace();
    int nt = data->getNt();
    int nx, ny, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    nz = this->getNz() + 2*lpml;

    // Get correct map (data or receiver mapping)
    if(maptype == 0) {
        map = (data->getGeom())->getSmap();
    }else{
        map = (data->getGeom())->getGmap();
    }

    dataarray = data->getData();
    int i;
    Index I(nx, ny, nz);
    Index Idat(nt, ntrace);
    switch(field)
    {
        case 0:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    dataarray[Idat(it,i)] = P2[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)];
                }
            }
        case 1:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    dataarray[Idat(it,i)] = Ax[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)];
                }
            }
            break;
        case 2:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    dataarray[Idat(it,i)] = Ay[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)];
                }
            }
            break;

        case 3:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    dataarray[Idat(it,i)] = Az[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)];
                }
            }
            break;
        default:
            for (i=0; i < ntrace; i++) 
            { 

                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    dataarray[Idat(it,i)] = P2[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)];
                }
            }
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

// =============== 2D ELASTIC WAVES CLASS =============== //
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
void WavesElastic2D<T>::forwardstepVelocity(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
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
    rockseis::Index I(nx,nz);
    rockseis::Index I_lr(lpml,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,lpml); // index for top and bottom pml zone
    
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I(ix,iz)] += dt*Rx[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxx_left[I_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iz)];
            Vx[I(ix,iz)] -= dt*Rx[I(ix,iz)]*(Pml->Sxx_left[I_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxx_right[I_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iz)];
            Vx[I(i,iz)] -= dt*Rx[I(i,iz)]*(Pml->Sxx_right[I_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iz)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vx[I(ix,iz)] += dt*Rx[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Sxzz_top[I_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I_tb(ix,iz)] + Pml->A_ltf[iz]*df[I(ix,iz)];
            
            Vx[I(ix,iz)] -= dt*Rx[I(ix,iz)]*(Pml->Sxzz_top[I_tb(ix,iz)] + Pml->C_ltf[iz]*df[I(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Sxzz_bottom[I_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I_tb(ix,iz)] + Pml->A_rbb[iz]*df[I(ix,i)];
            Vx[I(ix,i)] -= dt*Rx[I(ix,i)]*(Pml->Sxzz_bottom[I_tb(ix,iz)] + Pml->C_rbb[iz]*df[I(ix,i)]);
        }
    }
    
    
    // Derivate Sxz backward with respect to x
    der->ddx_bw(Sxz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I(ix,iz)] += dt*Rz[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Sxzx_left[I_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I_lr(ix,iz)] + Pml->A_ltf[ix]*df[I(ix,iz)];
            Vz[I(ix,iz)] -= dt*Rz[I(ix,iz)]*(Pml->Sxzx_left[I_lr(ix,iz)] + Pml->C_ltf[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Sxzx_right[I_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I_lr(ix,iz)] + Pml->A_rbb[ix]*df[I(i,iz)];
            Vz[I(i,iz)] -= dt*Rz[I(i,iz)]*(Pml->Sxzx_right[I_lr(ix,iz)] + Pml->C_rbb[ix]*df[I(i,iz)]);
        }
    }
    
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Vz[I(ix,iz)] += dt*Rz[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Szz_top[I_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iz)];
            
            Vz[I(ix,iz)] -= dt*Rz[I(ix,iz)]*(Pml->Szz_top[I_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Szz_bottom[I_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,i)];
            Vz[I(ix,i)] -= dt*Rz[I(ix,i)]*(Pml->Szz_bottom[I_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,i)]);
        }
    }
    
} // End of forwardstepVelocity


template<typename T>
void WavesElastic2D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iz, nx, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *L, *L2M, *M, *df;
    L = model->getL();
    L2M = model->getL2M();
    M = model->getM();
    df = der->getDf();
    rockseis::Index I(nx,nz);
    rockseis::Index I_lr(lpml,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,lpml); // index for top and bottom pml zone
    
    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I(ix,iz)] += dt*L2M[I(ix,iz)]*df[I(ix,iz)];
            Szz[I(ix,iz)] += dt*L[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Vxx_left[I_lr(ix,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I_lr(ix,iz)] + Pml->A_ltf[ix]*df[I(ix,iz)];
            
            Sxx[I(ix,iz)] -= dt*L2M[I(ix,iz)]*(Pml->Vxx_left[I_lr(ix,iz)] + Pml->C_ltf[ix]*df[I(ix,iz)]);
            Szz[I(ix,iz)] -= dt*L[I(ix,iz)]*(Pml->Vxx_left[I_lr(ix,iz)] + Pml->C_ltf[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Vxx_right[I_lr(ix,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I_lr(ix,iz)] + Pml->A_rbb[ix]*df[I(i,iz)];
            Sxx[I(i,iz)] -= dt*L2M[I(i,iz)]*(Pml->Vxx_right[I_lr(ix,iz)] + Pml->C_rbb[ix]*df[I(i,iz)]);
            Szz[I(i,iz)] -= dt*L[I(i,iz)]*(Pml->Vxx_right[I_lr(ix,iz)] + Pml->C_rbb[ix]*df[I(i,iz)]);
        }
    }
    
    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx and Szz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxx[I(ix,iz)] += dt*L[I(ix,iz)]*df[I(ix,iz)];
            Szz[I(ix,iz)] += dt*L2M[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Vzz_top[I_tb(ix,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I_tb(ix,iz)] + Pml->A_ltf[iz]*df[I(ix,iz)];
            
            Sxx[I(ix,iz)] -= dt*L[I(ix,iz)]*(Pml->Vzz_top[I_tb(ix,iz)] + Pml->C_ltf[iz]*df[I(ix,iz)]);
            Szz[I(ix,iz)] -= dt*L2M[I(ix,iz)]*(Pml->Vzz_top[I_tb(ix,iz)] + Pml->C_ltf[iz]*df[I(ix,iz)]);
            i = iz + nz - lpml;
            //Bottom
            Pml->Vzz_bottom[I_tb(ix,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I_tb(ix,iz)] + Pml->A_rbb[iz]*df[I(ix,i)];
            Sxx[I(ix,i)] -= dt*L[I(ix,i)]*(Pml->Vzz_bottom[I_tb(ix,iz)] + Pml->C_rbb[iz]*df[I(ix,i)]);
            Szz[I(ix,i)] -= dt*L2M[I(ix,i)]*(Pml->Vzz_bottom[I_tb(ix,iz)] + Pml->C_rbb[iz]*df[I(ix,i)]);
        }
    }
    
    
    // Derivate Vz forward with respect to x
    der->ddx_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I(ix,iz)] += dt*M[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < lpml; ix++){
            // Left
            Pml->Vzx_left[I_lr(ix,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I_lr(ix,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iz)];
            Sxz[I(ix,iz)] -= dt*M[I(ix,iz)]*(Pml->Vzx_left[I_lr(ix,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iz)]);
            // Right
            i = ix + nx - lpml;
            Pml->Vzx_right[I_lr(ix,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I_lr(ix,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iz)];
            Sxz[I(i,iz)] -= dt*M[I(i,iz)]*(Pml->Vzx_right[I_lr(ix,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iz)]);
        }
    }
    
    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(ix=0; ix < nx; ix++){
            Sxz[I(ix,iz)] += dt*M[I(ix,iz)]*df[I(ix,iz)];
        }
    }
    
    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(ix=0; ix < nx; ix++){
            // Top
            Pml->Vxz_top[I_tb(ix,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I_tb(ix,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iz)];
            Sxz[I(ix,iz)] -= dt*M[I(ix,iz)]*(Pml->Vxz_top[I_tb(ix,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iz)]);
            //Bottom
            i = iz + nz - lpml;
            Pml->Vxz_bottom[I_tb(ix,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I_tb(ix,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,i)];
            Sxz[I(ix,i)] -= dt*M[I(ix,i)]*(Pml->Vxz_bottom[I_tb(ix,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,i)]);
        }
    }
}


template<typename T>
void WavesElastic2D<T>::insertSource(std::shared_ptr<rockseis::ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, int sourcetype, int maptype, int it){
    Point2D<int> *map;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    nz = this->getNz() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype == 0) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }

    wav = source->getData();
    int i;
    Index I(nx, nz); //Model and Field indexes
    Index Idat(nt, ntrace); // Data indexes
    switch(sourcetype)
    {
        case 0:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    Sxx[I(lpml + map[i].x, lpml + map[i].y)] += dt*wav[Idat(it,i)]; 
                    Szz[I(lpml + map[i].x, lpml + map[i].y)] += dt*wav[Idat(it,i)]; 
                }
            }
            break;
	    case 1:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    Vx[I(lpml + map[i].x, lpml + map[i].y)] += dt*Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 3:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    Vz[I(lpml + map[i].x, lpml + map[i].y)] += dt*Mod[I(lpml + map[i].x, lpml + map[i].y)]*wav[Idat(it,i)]; 
                }
            }
            break;
        default:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0){
                    Sxx[I(lpml + map[i].x, lpml + map[i].y)] += dt*wav[Idat(it,i)]; 
                    Szz[I(lpml + map[i].x, lpml + map[i].y)] += dt*wav[Idat(it,i)]; 
                }
            }
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
    rockseis::Index I(nx,ny,nz);
    rockseis::Index I_lr(lpml,ny,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,ny,lpml); // index for top and bottom pml zone
    rockseis::Index I_fb(nx,lpml,nz); // index for front and back pml zone
    
   //////////////////////////// VX //////////////////////////
   //
    // Derivate Sxx forward with respect to x
    der->ddx_fw(Sxx);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I(ix,iy,iz)] += dt*Rx[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxx_left[I_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Sxx_left[I_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iy,iz)];
                
                Vx[I(ix,iy,iz)] -= dt*Rx[I(ix,iy,iz)]*(Pml->Sxx_left[I_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxx_right[I_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Sxx_right[I_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iy,iz)];
                Vx[I(i,iy,iz)] -= dt*Rx[I(i,iy,iz)]*(Pml->Sxx_right[I_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iy,iz)]);
            }
        }
    }

    // Derivate Sxy backward with respect to y
    der->ddy_bw(Sxy);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I(ix,iy,iz)] += dt*Rx[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Sxyy_front[I_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Sxyy_front[I_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I(ix,iy,iz)];
                
                Vx[I(ix,iy,iz)] -= dt*Rx[I(ix,iy,iz)]*(Pml->Sxyy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                i = iy + ny - lpml;
                //Back
                Pml->Sxyy_back[I_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Sxyy_back[I_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I(ix,i,iz)];
                Vx[I(ix,i,iz)] -= dt*Rx[I(ix,i,iz)]*(Pml->Sxyy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Sxy backward with respect to z
    der->ddz_bw(Sxz);
    // Compute Vx
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vx[I(ix,iy,iz)] += dt*Rx[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Sxzz_top[I_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Sxzz_top[I_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I(ix,iy,iz)];
                
                Vx[I(ix,iy,iz)] -= dt*Rx[I(ix,iy,iz)]*(Pml->Sxzz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Sxzz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Sxzz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I(ix,iy,i)];
                Vx[I(ix,iy,i)] -= dt*Rx[I(ix,iy,i)]*(Pml->Sxzz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
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
                Vy[I(ix,iy,iz)] += dt*Ry[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxyx_left[I_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxyx_left[I_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I(ix,iy,iz)];
                
                Vy[I(ix,iy,iz)] -= dt*Ry[I(ix,iy,iz)]*(Pml->Sxyx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxyx_right[I_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxyx_right[I_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I(i,iy,iz)];
                Vy[I(i,iy,iz)] -= dt*Ry[I(i,iy,iz)]*(Pml->Sxyx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
            }
        }
    }


    // Derivate Syy forward with respect to y
    der->ddy_fw(Syy);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I(ix,iy,iz)] += dt*Ry[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syy_front[I_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Syy_front[I_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I(ix,iy,iz)];
                
                Vy[I(ix,iy,iz)] -= dt*Ry[I(ix,iy,iz)]*(Pml->Syy_front[I_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Syy_back[I_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Syy_back[I_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I(ix,i,iz)];
                Vy[I(ix,i,iz)] -= dt*Ry[I(ix,i,iz)]*(Pml->Syy_back[I_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I(ix,i,iz)]);
            }
        }
    }
 

    // Derivate Syz backward with respect to z
    der->ddz_bw(Syz);
    // Compute Vy
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vy[I(ix,iy,iz)] += dt*Ry[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate bottom and top using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Syzz_top[I_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Syzz_top[I_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I(ix,iy,iz)];
                
                Vy[I(ix,iy,iz)] -= dt*Ry[I(ix,iy,iz)]*(Pml->Syzz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                //Bottom
                i = iz + nz - lpml;
                Pml->Syzz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Syzz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I(ix,iy,i)];
                Vy[I(ix,iy,i)] -= dt*Ry[I(ix,iy,i)]*(Pml->Syzz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
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
                Vz[I(ix,iy,iz)] += dt*Rz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Sxzx_left[I_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Sxzx_left[I_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I(ix,iy,iz)];
                
                Vz[I(ix,iy,iz)] -= dt*Rz[I(ix,iy,iz)]*(Pml->Sxzx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Sxzx_right[I_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Sxzx_right[I_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I(i,iy,iz)];
                Vz[I(i,iy,iz)] -= dt*Rz[I(i,iy,iz)]*(Pml->Sxzx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
            }
        }
    }

    // Derivate Syz backward with respect to y
    der->ddy_bw(Syz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I(ix,iy,iz)] += dt*Rz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Syzy_front[I_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Syzy_front[I_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I(ix,iy,iz)];
                
                Vz[I(ix,iy,iz)] -= dt*Rz[I(ix,iy,iz)]*(Pml->Syzy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Syzy_back[I_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Syzy_back[I_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I(ix,i,iz)];
                Vz[I(ix,i,iz)] -= dt*Rz[I(ix,i,iz)]*(Pml->Syzy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
            }
        }
    }
 
    // Derivate Szz forward with respect to z
    der->ddz_fw(Szz);
    // Compute Vz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Vz[I(ix,iy,iz)] += dt*Rz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate bottom and top using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Szz_top[I_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Szz_top[I_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iy,iz)];
                Vz[I(ix,iy,iz)] -= dt*Rz[I(ix,iy,iz)]*(Pml->Szz_top[I_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Szz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Szz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,iy,i)];
                Vz[I(ix,iy,i)] -= dt*Rz[I(ix,iy,i)]*(Pml->Szz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,iy,i)]);
            }
        }
    }
    
    
} // End of forwardstepVelocity

template<typename T>
void WavesElastic3D<T>::forwardstepStress(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Der<T>> der){
    int i, ix, iy, iz, nx, ny, nz, lpml;
    T dt;
    lpml = model->getLpml();
    nx = model->getNx() + 2*lpml;
    ny = model->getNy() + 2*lpml;
    nz = model->getNz() + 2*lpml;
    dt = this->getDt();
    T *L, *L2M, *M_xz, *M_yz, *M_xy, *df;
    L = model->getL();
    L2M = model->getL2M();
    M_xy = model->getM_xy();
    M_xz = model->getM_xz();
    M_yz = model->getM_yz();
    df = der->getDf();
    rockseis::Index I(nx,ny,nz);
    rockseis::Index I_lr(lpml,ny,nz); // index for left and right pml zone
    rockseis::Index I_tb(nx,ny,lpml); // index for top and bottom pml zone
    rockseis::Index I_fb(nx,lpml,nz); // index for front and back pml zone
    
    // Derivate Vx backward with respect to x
    der->ddx_bw(Vx);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I(ix,iy,iz)] += dt*L2M[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Syy[I(ix,iy,iz)] += dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Szz[I(ix,iy,iz)] += dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate left and right using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vxx_left[I_lr(ix,iy,iz)] = Pml->B_ltf[ix]*Pml->Vxx_left[I_lr(ix,iy,iz)] + Pml->A_ltf[ix]*df[I(ix,iy,iz)];
                
                Sxx[I(ix,iy,iz)] -= dt*L2M[I(ix,iy,iz)]*(Pml->Vxx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                Syy[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vxx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                Szz[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vxx_left[I_lr(ix,iy,iz)] + Pml->C_ltf[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Vxx_right[I_lr(ix,iy,iz)] = Pml->B_rbb[ix]*Pml->Vxx_right[I_lr(ix,iy,iz)] + Pml->A_rbb[ix]*df[I(i,iy,iz)];
                Sxx[I(i,iy,iz)] -= dt*L2M[I(i,iy,iz)]*(Pml->Vxx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
                Syy[I(i,iy,iz)] -= dt*L[I(i,iy,iz)]*(Pml->Vxx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
                Szz[I(i,iy,iz)] -= dt*L[I(i,iy,iz)]*(Pml->Vxx_right[I_lr(ix,iy,iz)] + Pml->C_rbb[ix]*df[I(i,iy,iz)]);
            }
        }
    }
    
    // Derivate Vy backward with respect to y
    der->ddy_bw(Vy);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I(ix,iy,iz)] +=  dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Syy[I(ix,iy,iz)] +=  dt*L2M[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Szz[I(ix,iy,iz)] +=  dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate front and back using non-staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vyy_front[I_fb(ix,iy,iz)] = Pml->B_ltf[iy]*Pml->Vyy_front[I_fb(ix,iy,iz)] + Pml->A_ltf[iy]*df[I(ix,iy,iz)];
                
                Sxx[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vyy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                Syy[I(ix,iy,iz)] -= dt*L2M[I(ix,iy,iz)]*(Pml->Vyy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                Szz[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vyy_front[I_fb(ix,iy,iz)] + Pml->C_ltf[iy]*df[I(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Vyy_back[I_fb(ix,iy,iz)] = Pml->B_rbb[iy]*Pml->Vyy_back[I_fb(ix,iy,iz)] + Pml->A_rbb[iy]*df[I(ix,i,iz)];
                Sxx[I(ix,i,iz)] -= dt*L[I(ix,i,iz)]*(Pml->Vyy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
                Syy[I(ix,i,iz)] -= dt*L2M[I(ix,i,iz)]*(Pml->Vyy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
                Szz[I(ix,i,iz)] -= dt*L[I(ix,i,iz)]*(Pml->Vyy_back[I_fb(ix,iy,iz)] + Pml->C_rbb[iy]*df[I(ix,i,iz)]);
            }
        }
    }
    
    // Derivate Vz backward with respect to z
    der->ddz_bw(Vz);
    // Compute Sxx,Syy,Szz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxx[I(ix,iy,iz)] +=  dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Syy[I(ix,iy,iz)] +=  dt*L[I(ix,iy,iz)]*df[I(ix,iy,iz)];
                Szz[I(ix,iy,iz)] +=  dt*L2M[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }
    
    // Attenuate top and bottom using non-staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vzz_top[I_tb(ix,iy,iz)] = Pml->B_ltf[iz]*Pml->Vzz_top[I_tb(ix,iy,iz)] + Pml->A_ltf[iz]*df[I(ix,iy,iz)];
                
                Sxx[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vzz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                Syy[I(ix,iy,iz)] -= dt*L[I(ix,iy,iz)]*(Pml->Vzz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                Szz[I(ix,iy,iz)] -= dt*L2M[I(ix,iy,iz)]*(Pml->Vzz_top[I_tb(ix,iy,iz)] + Pml->C_ltf[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Vzz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb[iz]*Pml->Vzz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb[iz]*df[I(ix,iy,i)];
                Sxx[I(ix,iy,i)] -= dt*L[I(ix,iy,i)]*(Pml->Vzz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
                Syy[I(ix,iy,i)] -= dt*L[I(ix,iy,i)]*(Pml->Vzz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
                Szz[I(ix,iy,i)] -= dt*L2M[I(ix,iy,i)]*(Pml->Vzz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb[iz]*df[I(ix,iy,i)]);
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
                Sxz[I(ix,iy,iz)] +=  dt*M_xz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vzx_left[I_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vzx_left[I_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iy,iz)];
                
                Sxz[I(ix,iy,iz)] -= dt*M_xz[I(ix,iy,iz)]*(Pml->Vzx_left[I_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Vzx_right[I_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vzx_right[I_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iy,iz)];
                Sxz[I(i,iy,iz)] -= dt*M_xz[I(i,iy,iz)]*(Pml->Vzx_right[I_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iy,iz)]);
            }
        }
    }
 

    // Derivate Vx forward with respect to z
    der->ddz_fw(Vx);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxz[I(ix,iy,iz)] +=  dt*M_xz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vxz_top[I_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vxz_top[I_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iy,iz)];
                
                Sxz[I(ix,iy,iz)] -= dt*M_xz[I(ix,iy,iz)]*(Pml->Vxz_top[I_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Vxz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vxz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,iy,i)];
                Sxz[I(ix,iy,i)] -= dt*M_xz[I(ix,iy,i)]*(Pml->Vxz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,iy,i)]);
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
                Syz[I(ix,iy,iz)] +=  dt*M_yz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate top and bottom using staggered variables
    for(iz=0; iz < lpml; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                // Top
                Pml->Vyz_top[I_tb(ix,iy,iz)] = Pml->B_ltf_stag[iz]*Pml->Vyz_top[I_tb(ix,iy,iz)] + Pml->A_ltf_stag[iz]*df[I(ix,iy,iz)];
                
                Syz[I(ix,iy,iz)] -= dt*M_yz[I(ix,iy,iz)]*(Pml->Vyz_top[I_tb(ix,iy,iz)] + Pml->C_ltf_stag[iz]*df[I(ix,iy,iz)]);
                i = iz + nz - lpml;
                //Bottom
                Pml->Vyz_bottom[I_tb(ix,iy,iz)] = Pml->B_rbb_stag[iz]*Pml->Vyz_bottom[I_tb(ix,iy,iz)] + Pml->A_rbb_stag[iz]*df[I(ix,iy,i)];
                Syz[I(ix,iy,i)] -= dt*M_yz[I(ix,iy,i)]*(Pml->Vyz_bottom[I_tb(ix,iy,iz)] + Pml->C_rbb_stag[iz]*df[I(ix,iy,i)]);
            }
        }
    }

    // Derivate Vz forward with respect to y
    der->ddy_fw(Vz);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Syz[I(ix,iy,iz)] +=  dt*M_yz[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vzy_front[I_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vzy_front[I_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I(ix,iy,iz)];
                
                Syz[I(ix,iy,iz)] -= dt*M_yz[I(ix,iy,iz)]*(Pml->Vzy_front[I_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Vzy_back[I_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vzy_back[I_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I(ix,i,iz)];
                Syz[I(ix,i,iz)] -= dt*M_yz[I(ix,i,iz)]*(Pml->Vzy_back[I_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I(ix,i,iz)]);
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
                Sxy[I(ix,iy,iz)] +=  dt*M_xy[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate front and back using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < lpml; iy++){
            for(ix=0; ix < nx; ix++){
                // Front
                Pml->Vxy_front[I_fb(ix,iy,iz)] = Pml->B_ltf_stag[iy]*Pml->Vxy_front[I_fb(ix,iy,iz)] + Pml->A_ltf_stag[iy]*df[I(ix,iy,iz)];
                
                Sxy[I(ix,iy,iz)] -= dt*M_xy[I(ix,iy,iz)]*(Pml->Vxy_front[I_fb(ix,iy,iz)] + Pml->C_ltf_stag[iy]*df[I(ix,iy,iz)]);
                //Back
                i = iy + ny - lpml;
                Pml->Vxy_back[I_fb(ix,iy,iz)] = Pml->B_rbb_stag[iy]*Pml->Vxy_back[I_fb(ix,iy,iz)] + Pml->A_rbb_stag[iy]*df[I(ix,i,iz)];
                Sxy[I(ix,i,iz)] -= dt*M_xy[I(ix,i,iz)]*(Pml->Vxy_back[I_fb(ix,iy,iz)] + Pml->C_rbb_stag[iy]*df[I(ix,i,iz)]);
            }
        }
    }

    // Derivate Vy forward with respect to x
    der->ddx_fw(Vy);
    // Compute Sxz
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < nx; ix++){
                Sxy[I(ix,iy,iz)] +=  dt*M_xy[I(ix,iy,iz)]*df[I(ix,iy,iz)];
            }
        }
    }

    // Attenuate left and right using staggered variables
    for(iz=0; iz < nz; iz++){
        for(iy=0; iy < ny; iy++){
            for(ix=0; ix < lpml; ix++){
                // Left
                Pml->Vyx_left[I_lr(ix,iy,iz)] = Pml->B_ltf_stag[ix]*Pml->Vyx_left[I_lr(ix,iy,iz)] + Pml->A_ltf_stag[ix]*df[I(ix,iy,iz)];
                
                Sxy[I(ix,iy,iz)] -= dt*M_xy[I(ix,iy,iz)]*(Pml->Vyx_left[I_lr(ix,iy,iz)] + Pml->C_ltf_stag[ix]*df[I(ix,iy,iz)]);
                // Right
                i = ix + nx - lpml;
                Pml->Vyx_right[I_lr(ix,iy,iz)] = Pml->B_rbb_stag[ix]*Pml->Vyx_right[I_lr(ix,iy,iz)] + Pml->A_rbb_stag[ix]*df[I(i,iy,iz)];
                Sxy[I(i,iy,iz)] -= dt*M_xy[I(i,iy,iz)]*(Pml->Vyx_right[I_lr(ix,iy,iz)] + Pml->C_rbb_stag[ix]*df[I(i,iy,iz)]);
            }
        }
    }
 

}

template<typename T>
void WavesElastic3D<T>::insertSource(std::shared_ptr<rockseis::ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, int sourcetype, int maptype, int it){
    Point3D<int> *map;
    T *wav; 
    int ntrace = source->getNtrace();
    int nx, ny, nz, lpml;
    lpml = this->getLpml();
    nx = this->getNx() + 2*lpml;
    ny = this->getNy() + 2*lpml;
    nz = this->getNz() + 2*lpml;
    T *Mod;
    int nt = this->getNt();
    T dt = this->getDt();

    // Get correct map (source or receiver mapping)
    if(maptype ==0) {
        map = (source->getGeom())->getSmap();
    }else{
        map = (source->getGeom())->getGmap();
    }
    wav = source->getData();
    int i;
    Index I(nx, ny, nz); //Model and Field indexes
    Index Idat(nt, ntrace); // Data indexes
    switch(sourcetype)
    {
        case 0:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Sxx[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                    Syy[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                    Szz[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                }
            }
            break;
	    case 1:
            Mod = model->getRx();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Vx[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 2:
            Mod = model->getRy();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Vy[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        case 3:
            Mod = model->getRz();
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Vz[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*Mod[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)]*wav[Idat(it,i)]; 
                }
            }
            break;
        default:
            for (i=0; i < ntrace; i++) 
            { 
                if(map[i].x >= 0 && map[i].y >=0 && map[i].z >=0)
                {
                    Sxx[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                    Syy[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                    Szz[I(lpml + map[i].x, lpml + map[i].y, lpml + map[i].z)] += dt*wav[Idat(it,i)]; 
                }
            }
            break;

    }
}

}


