#include "domain.h"

namespace rockseis {
// constructor
template<typename T>
Domain<T>::Domain()
{
    geometry = std::make_shared<Geometry<T>>(); 
    lpad = 0;
    low = -1;
    high = -1;
    padl = 0;
    padh = 0;
    nd = 1;
    d = 0;
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

    int tmp[] = {nx,ny,nz}; 
    std::sort(tmp, tmp+3);
    if (tmp[2] == nz) this->setDim(2);
    if (tmp[2] == ny) this->setDim(1);
    if (tmp[2] == nx) this->setDim(0);

    int nxdom, nydom, nzdom;
    int nxpad, nypad, nzpad;
    int ix0, iy0, iz0;
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

            break;
        default:
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

            break;
    }

    this->setNx(nxdom);
    this->setNx_pad(nxpad);

    this->setNy(nydom);
    this->setNy_pad(nypad);

    this->setNz(nzdom);
    this->setNz_pad(nzpad);

    this->setIx0(ix0);
    this->setIy0(iy0);
    this->setIz0(iz0);
}


template<typename T>
Domain<T>::~Domain() {
   // Do nothing
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
