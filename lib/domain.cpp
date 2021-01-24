#include "domain.h"

namespace rockseis {
// constructor
template<typename T>
Domain<T>::Domain()
{
    lpad = 0;
    low = -1;
    high = -1;
    padl = 0;
    padh = 0;
    nd = 1;
    d = 0;
    status = false;
    allocated = false;
}

template<typename T>
void Domain<T>::setupDomain(const int nx, const int ny, const int nz, const int _d, const int _nd, const int order)
{
    if(_d > _nd-1) rs_error("Domain<T>::setupDomain:Domain number larger than number of domains");
    if(_d < 0) rs_error("Domain<T>::setupDomain:Domain number less than zero.");
    this->setNd(_nd);
    this->setD(_d);
    this->setLpad(order+1);
    if(_d > 0){
        this->setLow(_d-1);
        this->setPadl(this->getLpad());
    }else{
        this->setLow(-1);
        this->setPadl(0);
    }
    if(_d < _nd-1){
        this->setHigh(_d+1);
        this->setPadh(this->getLpad());
    }else{
        this->setHigh(-1);
        this->setPadh(0);
    }
    this->setNx_orig(nx);
    this->setNy_orig(ny);
    this->setNz_orig(nz);

    int tmp[] = {nx,ny,nz}; 
    std::sort(tmp, tmp+3);
    if (tmp[2] == nz) this->setDim(2);
    if (tmp[2] == ny) this->setDim(1);
    if (tmp[2] == nx) this->setDim(0);

    int nxdom=0, nydom=0, nzdom=0;
    int nxpad=0, nypad=0, nzpad=0;
    int ix0=0, iy0=0, iz0=0;
    switch(this->getDim()){
        case 0:
            nxdom = std::ceil((T) nx/((T) _nd));

            //Find global coordinates of the domain
            ix0 = _d*nxdom - this->getPadl();
            iy0 = 0;
            iz0 = 0;

            // Treat special case of last domain which can be smaller than the other domains
            if (_d == _nd-1){ 
                nxdom = nx - nxdom*(_nd-1);
            }
            if(nxdom <= 0){
                rs_error("Domain<T>::setupDomain:Number of domains cannot exceed the number of grid points");
            }
            nydom = ny;
            nzdom = nz;
            nxpad = nxdom+this->getPadl()+this->getPadh();
            nypad = ny;
            nzpad = nz;
            wrk = (T *) calloc(nz*ny*(this->getLpad()), sizeof(T));
            if(wrk == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk array");
            break;
        case 1:
            nydom = std::ceil((T) ny/((T) _nd));

            //Find global coordinates of the domain
            ix0 = 0;
            iy0 = _d*nydom - this->getPadl();
            iz0 = 0;
            // Treat special case of last domain which can be smaller than the other domains
            if (_d == _nd-1){ 
                nydom = ny - nydom*(_nd-1);
            }
            if(nydom <= 0){
                rs_error("Domain<T>::setupDomain:Number of domains cannot exceed the number of grid points");
            }
            nxdom = nx;
            nzdom = nz;

            nxpad = nxdom;
            nypad = nydom+this->getPadl()+this->getPadh();
            nzpad = nzdom;
            wrk = (T *) calloc(nx*nz*(this->getLpad()), sizeof(T));
            if(wrk == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk array");
            break;
        case 2:
            nzdom = std::ceil((T) nz/((T) _nd));

            //Find global coordinates of the domain
            ix0 = 0;
            iy0 = 0;
            iz0 = _d*nzdom - this->getPadl();
            // Treat special case of last domain which can be smaller than the other domains
            if (_d == _nd-1){ 
                nzdom = nz - nzdom*(_nd-1);
            }
            if(nzdom <= 0){
                rs_error("Domain<T>::setupDomain:Number of domains cannot exceed the number of grid points");
            }
            nxdom = nx;
            nydom = ny;

            nxpad = nxdom;
            nypad = nydom;
            nzpad = nzdom+this->getPadl()+this->getPadh();
            wrk = (T *) calloc(nx*ny*(this->getLpad()), sizeof(T));
            if(wrk == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk array");
            break;
        default:
            rs_error("Domain<T>::setupDomain:Invalid dimension");
            break;
    }

    this->setNx_pad(nxpad);
    this->setNy_pad(nypad);
    this->setNz_pad(nzpad);

    this->setIx0(ix0);
    this->setIy0(iy0);
    this->setIz0(iz0);
    this->status = true;
    this->allocated = true;
}


template<typename T>
void Domain<T>::copyFromboundary(const bool side, const T *array){
    int ix,iy,iz,io;
    if(side){
        if(this->getHigh() < 0) rs_error("Domain<T>::copyFromboundary: Invalid side.");
    }else{
        if(this->getLow() < 0) rs_error("Domain<T>::copyFromboundary: Invalid side.");
    }
    int pad=0, n0=0, n1=0, n2=0, nwrk=0;

    n0 = this->getNx_pad();
    n1 = this->getNy_pad();
    n2 = this->getNz_pad();
    pad = this->getLpad();

    switch(this->getDim()){
        case 0:
            nwrk = n1;
            if(side){
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk[IDWRK(io,iy,iz)] = array[IDARRAY(n0-2*pad + io, iy, iz)];
                        }
                    }
                }
            }else{
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk[IDWRK(io,iy,iz)] = array[IDARRAY(io + pad, iy, iz)];
                        }
                    }
                }
            }
            break;
        case 1:
            nwrk = n0;
            if(side){
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk[IDWRK(io,ix,iz)] = array[IDARRAY(ix, n1-2*pad + io, iz)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk[IDWRK(io,ix,iz)] = array[IDARRAY(ix, io + pad, iz)];
                        }
                    }
                }
            }
            break;
        case 2:
            nwrk = n0;
            if(side){
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            wrk[IDWRK(io,ix,iy)] = array[IDARRAY(ix, iy, n2-2*pad + io)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            wrk[IDWRK(io,ix,iy)] = array[IDARRAY(ix, iy, io + pad)];
                        }
                    }
                }
            }
            break;
        default:
            rs_error("Domain<T>::copyBoundary:Invalid dimension");
            break;
            break;
    }
}

template<typename T>
void Domain<T>::copyToboundary(const bool side, T *array){
    int ix,iy,iz,io;
    if(side){
        if(this->getHigh() < 0) rs_error("Domain<T>::copyToboundary: Invalid side.");
    }else{
        if(this->getLow() < 0) rs_error("Domain<T>::copyToboundary: Invalid side.");
    }
    int pad=0, n0=0, n1=0, n2=0, nwrk=0;

    n0 = this->getNx_pad();
    n1 = this->getNy_pad();
    n2 = this->getNz_pad();
    pad = this->getLpad();

    switch(this->getDim()){
        case 0:
            nwrk = n1;
            if(side){
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                           array[IDARRAY(n0-pad + io, iy, iz)] = wrk[IDWRK(io,iy,iz)];
                        }
                    }
                }
            }else{
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                            array[IDARRAY(io, iy, iz)] = wrk[IDWRK(io,iy,iz)];
                        }
                    }
                }
            }
            break;
        case 1:
            nwrk = n0;
            if(side){
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                            array[IDARRAY(ix, n1-pad + io, iz)] = wrk[IDWRK(io,ix,iz)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                          array[IDARRAY(ix, io, iz)] = wrk[IDWRK(io,ix,iz)];
                        }
                    }
                }
            }
            break;
        case 2:
            nwrk = n0;
            if(side){
                pad = this->getPadh();
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            array[IDARRAY(ix, iy, n2-pad + io)] = wrk[IDWRK(io,ix,iy)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            array[IDARRAY(ix, iy, io)] = wrk[IDWRK(io,ix,iy)];
                        }
                    }
                }
            }
            break;
        default:
            rs_error("Domain<T>::copyBoundary:Invalid dimension");
            break;
            break;
    }

}



template<typename T>
Domain<T>::~Domain() {
    if(allocated){
        free(wrk);
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
