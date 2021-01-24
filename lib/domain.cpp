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
            wrk = (T *) calloc(nz*ny*(order+1), sizeof(T));
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
            wrk = (T *) calloc(nx*nz*(order+1), sizeof(T));
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
            wrk = (T *) calloc(nx*ny*(order+1), sizeof(T));
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
Domain<T>::~Domain() {
    if(allocated){
 //       free(wrk);
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
