#include "domain.h"

namespace rockseis {
// constructor
template<typename T>
Domain<T>::Domain()
{
    lpad = 0;
    low.x = -1;
    low.y = -1;
    low.z = -1;
    high.x = -1;
    high.y = -1;
    high.z = -1;
    padl.x = 0;
    padl.y = 0;
    padl.z = 0;
    padh.x = 0;
    padh.y = 0;
    padh.z = 0;
    nd = 1;
    nd0 = 1;
    nd1 = 1;
    nd2 = 1;
    d = 0;
    status = false;
    for (int i=0; i< 12; i++) edge[i] = -1;
    for (int i=0; i< 7; i++) corner[i] = -1;
    for (int i=0; i< 7; i++) allocated[i] = false;
    mpiset = false;
}
template<typename T>
int Domain<T>::getId(const int d, const int dim)
{
   int x,y,z,tmp;
   z = (int) (d /(this->getNd0()*this->getNd1()));
   tmp = d - z*(this->getNd0()*this->getNd1());
   y = (int) (tmp/this->getNd0());
   x = tmp % this->getNd0();

   switch (dim){
      case 0:
         return x; 
         break;
      case 1:
         return y;
         break;
      case 2:
         return z;
         break;
      default:
         rs_error("Domain<T>::setLow: Invalid dimension.");
         return 0;
         break;
   }
}

template<typename T>
void Domain<T>::setLow(const int val, const int dim)
{
   switch (dim){
      case 0:
         low.x = val;
         break;
      case 1:
         low.y = val;
         break;
      case 2:
         low.z = val;
         break;
      default:
         rs_error("Domain<T>::setLow: Invalid dimension.");
         break;
   }
}

template<typename T>
int Domain<T>::getLow(const int dim)
{
   switch (dim){
      case 0:
         return low.x;
         break;
      case 1:
         return low.y;
         break;
      case 2:
         return low.z;
         break;
      default:
         rs_error("Domain<T>::getLow: Invalid dimension.");
         return 0;
         break;
   }
}

template<typename T>
void Domain<T>::setHigh(const int val, const int dim)
{
   switch (dim){
      case 0:
         high.x = val;
         break;
      case 1:
         high.y = val;
         break;
      case 2:
         high.z = val;
         break;
      default:
         rs_error("Domain<T>::setHigh: Invalid dimension.");
         break;
   }
}

template<typename T>
int Domain<T>::getHigh(const int dim)
{
   switch (dim){
      case 0:
         return high.x;
         break;
      case 1:
         return high.y;
         break;
      case 2:
         return high.z;
         break;
      default:
         rs_error("Domain<T>::getHigh: Invalid dimension.");
         return 0;
         break;
   }
}

template<typename T>
void Domain<T>::setPadl(const int val, const int dim)
{
   switch (dim){
      case 0:
         padl.x = val;
         break;
      case 1:
         padl.y = val;
         break;
      case 2:
         padl.z = val;
         break;
      default:
         rs_error("Domain<T>::setPadl: Invalid dimension.");
         break;
   }
}

template<typename T>
int Domain<T>::getPadl(const int dim)
{
   switch (dim){
      case 0:
         return padl.x;
         break;
      case 1:
         return padl.y;
         break;
      case 2:
         return padl.z;
         break;
      default:
         rs_error("Domain<T>::getPadl: Invalid dimension.");
         return 0;
         break;
   }
}

template<typename T>
void Domain<T>::setPadh(const int val, const int dim)
{
   switch (dim){
      case 0:
         padh.x = val;
         break;
      case 1:
         padh.y = val;
         break;
      case 2:
         padh.z = val;
         break;
      default:
         rs_error("Domain<T>::setPadh: Invalid dimension.");
         break;
   }
}

template<typename T>
int Domain<T>::getPadh(const int dim)
{
   switch (dim){
      case 0:
         return padh.x;
         break;
      case 1:
         return padh.y;
         break;
      case 2:
         return padh.z;
         break;
      default:
         rs_error("Domain<T>::getPadh: Invalid dimension.");
         return 0;
         break;
   }
}

template<typename T>
void Domain<T>::setupDomain1D(const int nx, const int ny, const int nz, const int _d, const int _nd, const int order)
{
    if(_d > _nd-1) rs_error("Domain<T>::setupDomain:Domain number larger than number of domains");
    if(_d < 0) rs_error("Domain<T>::setupDomain:Domain number less than zero.");
    this->setNd(_nd);
    this->setD(_d);
    this->setLpad(order);
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
    if(tmp[2] == nz) {
       this->setDim(2);
       nd0 = 1;
       nd1 = 1;
       nd2 = _nd;
    }
    if(tmp[2] == ny) {
       this->setDim(1);
       nd0 = 1;
       nd1 = _nd;
       nd2 = 1;
    }
    if(tmp[2] == nx){
       this->setDim(0);
       nd0 = _nd;
       nd1 = 1;
       nd2 = 1;
    }


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
            if(nxdom < this->getLpad()){
                rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
            }
            nydom = ny;
            nzdom = nz;
            nxpad = nxdom+this->getPadl()+this->getPadh();
            nypad = ny;
            nzpad = nz;
            wrksize[0] = nz*ny*(this->getLpad());
            wrk0 = (T *) calloc(wrksize[0], sizeof(T));
            if(wrk0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk0 array");
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
            if(nydom < this->getLpad()){
                rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
            }
            nxdom = nx;
            nzdom = nz;

            nxpad = nxdom;
            nypad = nydom+this->getPadl()+this->getPadh();
            nzpad = nzdom;
            wrksize[0] = nx*nz*(this->getLpad());
            wrk0 = (T *) calloc(wrksize[0], sizeof(T));
            if(wrk0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk0 array");
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

            if(nzdom < this->getLpad()){
                rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
            }
            nxdom = nx;
            nydom = ny;

            nxpad = nxdom;
            nypad = nydom;
            nzpad = nzdom+this->getPadl()+this->getPadh();
            wrksize[0] = nx*ny*(this->getLpad());
            wrk0 = (T *) calloc(wrksize[0], sizeof(T));
            if(wrk0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating wrk0 array");
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
    this->allocated[0] = true;
}

template<typename T>
void Domain<T>::setupDomain3D(const int nx, const int ny, const int nz, const int _d, const int _nd0, const int _nd1, const int _nd2, const int order)
{
   int nd = _nd0*_nd1*_nd2;
    if(_d > nd-1) rs_error("Domain<T>::setupDomain:Domain number larger than number of domains");
    if(_d < 0) rs_error("Domain<T>::setupDomain:Domain number less than zero.");
    this->setNd(nd);
    this->setNd0(_nd0);
    this->setNd1(_nd1);
    this->setNd2(_nd2);
    this->setD(_d);
    this->setLpad(order);
    this->setNx_orig(nx);
    this->setNy_orig(ny);
    this->setNz_orig(nz);

    int nxdom=0, nydom=0, nzdom=0;
    int nxpad=0, nypad=0, nzpad=0;
    int ix0=0, iy0=0, iz0=0;

    nxdom = std::ceil((T) nx/((T) _nd0));

    // Find side neighbours
    if(getId(_d, 0) == 0){
       this->setPadl(0,0);
       this->setLow(-1,0);
    }else{
       this->setPadl(this->getLpad(),0);
       this->setLow(IDDOM(getId(_d, 0)-1, getId(_d, 1), getId(_d, 2)), 0);
    }

    if(getId(_d, 0) == _nd0-1){
       this->setPadh(0,0);
       this->setHigh(-1,0);
    }else{
       this->setPadh(this->getLpad(),0);
       this->setHigh(IDDOM(getId(_d, 0)+1, getId(_d, 1), getId(_d, 2)), 0);
    }

   if(getId(_d, 1) == 0){
       this->setPadl(0,1);
       this->setLow(-1,1);
    }else{
       this->setPadl(this->getLpad(),1);
       this->setLow(IDDOM(getId(_d, 0), getId(_d, 1)-1, getId(_d, 2)), 1);
    }

    if(getId(_d, 1) == _nd1-1){
       this->setPadh(0,1);
       this->setHigh(-1,1);
    }else{
       this->setPadh(this->getLpad(),1);
       this->setHigh(IDDOM(getId(_d, 0), getId(_d, 1)+1, getId(_d, 2)), 1);
    }

   if(getId(_d, 2) == 0){
       this->setPadl(0,2);
       this->setLow(-1,2);
    }else{
       this->setPadl(this->getLpad(),2);
       this->setLow(IDDOM(getId(_d, 0), getId(_d, 1), getId(_d, 2)-1), 2);
    }

    if(getId(_d, 2) == _nd2-1){
       this->setPadh(0,2);
       this->setHigh(-1,2);
    }else{
       this->setPadh(this->getLpad(),2);
       this->setHigh(IDDOM(getId(_d, 0), getId(_d, 1), getId(_d, 2)+1), 2);
    }

    // Find edges
    int id0,id1,id2;
    id0 = getId(_d, 0);
    id1 = getId(_d, 1);
    id2 = getId(_d, 2);
    int d0,d1,d2;
    for(int i=0; i<2; i++){
       d0 = 2*i -1;
       for(int j=0; j<2; j++){
          d1 = 2*j -1;
             if((((id0 + d0) >= 0) && ((id0 + d0) < _nd0)) && ((((id1 + d1) >= 0) && ((id1 + d1) < _nd1)))) {
                this->setEdge(IDDOM(id0+d0, id1+d1, id2), IDEDGE(i,j,0));
             }else{
                this->setEdge(-1, IDEDGE(i,j,0));
             }

             if((((id0 + d0) >= 0) && ((id0 + d0) < _nd0)) && (((id2 + d1) >= 0) && ((id2 + d1) < _nd2))){
                this->setEdge(IDDOM(id0+d0, id1, id2+d1), IDEDGE(i,j,1));
             }else{
                this->setEdge(-1, IDEDGE(i,j,1));
             }

             if((((id1 + d0) >= 0) && ((id1 + d0) < _nd1)) && (((id2 + d1) >= 0) && ((id2 + d1) < _nd2))){
                this->setEdge(IDDOM(id0, id1+d0, id2+d1), IDEDGE(i,j,2));
             }else{
                this->setEdge(-1, IDEDGE(i,j,2));
             }
       }
    }

    // Find corners
    for(int i=0; i<2; i++){
       d0 = 2*i -1;
       for(int j=0; j<2; j++){
          d1 = 2*j -1;
          for(int k=0; k<2; k++){
             d2 = 2*k -1;
             if((((id0 + d0) >= 0) && ((id0 + d0) < _nd0)) && (((id1 + d1) >= 0) && ((id1 + d1) < _nd1)) && (((id2 + d2) >= 0) && ((id2 + d2) < _nd2))){
                this->setCorner(IDDOM(id0+d0, id1+d1, id2+d2), IDCRN(i,j,k));
             }else{
                this->setCorner(-1, IDCRN(i,j,k));
             }
          }
       }
    }

    //Find global coordinates of the domain
    ix0 = getId(_d, 0)*nxdom - this->getPadl(0);

    // Treat special case of last domain which can be smaller than the other domains
    if (getId(_d, 0) == _nd0-1){ 
       nxdom = nx - nxdom*(_nd0-1);
    }
    if(nxdom < this->getLpad()){
       rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
    }
    nxpad = nxdom+this->getPadl(0)+this->getPadh(0);

    nydom = std::ceil((T) ny/((T) _nd1));

    //Find global coordinates of the domain
    iy0 = getId(_d, 1)*nydom - this->getPadl(1);
    // Treat special case of last domain which can be smaller than the other domains
    if (getId(_d, 1) == _nd1-1){ 
       nydom = ny - nydom*(_nd1-1);
    }
    if(nydom < this->getLpad()){
       rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
    }
    nypad = nydom+this->getPadl(1)+this->getPadh(1);

    nzdom = std::ceil((T) nz/((T) _nd2));

    iz0 = getId(_d, 2)*nzdom - this->getPadl(2);
    // Treat special case of last domain which can be smaller than the other domains
    if (getId(_d, 2) == _nd2-1){ 
       nzdom = nz - nzdom*(_nd2-1);
    }

    if(nzdom < this->getLpad()){
       rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
    }
    nzpad = nzdom+this->getPadl(2)+this->getPadh(2);

    this->setNx_pad(nxpad);
    this->setNy_pad(nypad);
    this->setNz_pad(nzpad);

    this->setIx0(ix0);
    this->setIy0(iy0);
    this->setIz0(iz0);
    this->status = true;

    // Allocate wrk arrays
    int nxwrk,nywrk,nzwrk;
    nxwrk = nx-this->getPadl(0)-this->getPadh(0);
    nywrk = ny-this->getPadl(1)-this->getPadh(1);
    nzwrk = nz-this->getPadl(2)-this->getPadh(2);

    if(this->getPadl(0) > 0 || this->getPadh(0) > 0){
       wrksize[0] = nywrk*nzwrk*this->getLpad();
       wrk0 = (T *) calloc(wrksize[0], sizeof(T));
       if(wrk0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrk0 array");
       this->allocated[0] = true;
    }

    if(this->getPadl(1) > 0 || this->getPadh(1) > 0){
       wrksize[1] = nxwrk*nzwrk*this->getLpad();
       wrk1 = (T *) calloc(wrksize[1], sizeof(T));
       if(wrk1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrk1 array");
       this->allocated[1] = true;
    }

    if(this->getPadl(2) > 0 || this->getPadh(2) > 0){
       wrksize[2] = nxwrk*nywrk*this->getLpad();
       wrk2 = (T *) calloc(wrksize[2], sizeof(T));
       if(wrk2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrk2 array");
       this->allocated[2] = true;
    }

    int edgsum=0;
    for(int i=0; i<4; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
    wrkedgsize[0] = nzwrk*this->getLpad()*this->getLpad();
       wrkedg0 = (T *) calloc(wrkedgsize[0], sizeof(T));
       if(wrkedg0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrkedg0 array");
       this->allocated[3] = true;
    }

    edgsum=0;
    for(int i=4; i<8; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
    wrkedgsize[1] = nywrk*this->getLpad()*this->getLpad();
       wrkedg1 = (T *) calloc(wrkedgsize[1], sizeof(T));
       if(wrkedg1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrkedg1 array");
       this->allocated[4] = true;
    }

    edgsum=0;
    for(int i=8; i<12; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
    wrkedgsize[2] = nxwrk*this->getLpad()*this->getLpad();
       wrkedg2 = (T *) calloc(wrkedgsize[2], sizeof(T));
       if(wrkedg2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrkedg2 array");
       this->allocated[5] = true;
    }

    int crnsum=0;
    for(int i=0; i<8; i++) crnsum += this->getCorner(i);
    if(crnsum > -8){
       wrkcrnsize = this->getLpad()*this->getLpad();
       if(nywrk > 1) wrkcrnsize *= this->getLpad();
       wrkcrn = (T *) calloc(wrkcrnsize, sizeof(T));
       if(wrkcrn == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating wrkcrn array");
       this->allocated[6] = true;
    }

}

template<typename T>
void Domain<T>::copyFromboundary2(const T *array){
   size_t ix,iy,iz,io;
   size_t padl0=0, padl1, padl2, padh0=0, padh1=0, padh2=0, n0=0, n1=0, n2=0, nwrk=0;

   padl0 = this->getPadl(0);
   padl1 = this->getPadl(1);
   padl2 = this->getPadl(2);
   padh0 = this->getPadh(0);
   padh1 = this->getPadh(1);
   padh2 = this->getPadh(2);

   n0 = this->getNx_pad() - padl0 - padh0;
   n1 = this->getNy_pad() - padl1 - padh1;
   n2 = this->getNz_pad() - padl2 - padh2;

   if(this->getLow(0)>0){
      nwrk = n1-padl1-padh1;
      for (io=0; io < pad; io++) {
         for (iy=0; iy < n1; iy++) {
            for (iz=0; iz < n2; iz++) {
               wrk0[IDWRK(io,iy,iz)] = array[IDARRAY(io + padl0, iy+padl1, iz+padl2)];
            }
         }
      }
   }
   if(this->getHigh(0)>0){
      nwrk = n1-padl1-padh1;
      for (io=0; io < pad; io++) {
         for (iy=0; iy < n1; iy++) {
            for (iz=0; iz < n2; iz++) {
               wrk0[IDWRK(io,iy,iz)] = array[IDARRAY(n0-padh0 + io, iy+padl1, iz+padl2)];
            }
         }
      }
   }

   if(this->getLow(1)>0){
      nwrk = n0-padl0-padh0;
      for (ix=0; ix < n0; ix++) {
         for (io=0; io < pad; io++) {
            for (iz=0; iz < n2; iz++) {
                            wrk1[IDWRK(io,ix,iz)] = array[IDARRAY(ix+padl0, io + padl1, iz+padl2)];
            }
         }
      }
   }
   if(this->getHigh(1)>0){
      nwrk = n0-padl0-padh0;
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
               wrk1[IDWRK(io,ix,iz)] = array[IDARRAY(ix + padl0, n1-padh1 + io, iz+padl2)];
            }
         }
      }
   }


}


template<typename T>
void Domain<T>::copyFromboundary(const bool side, const T *array){
    size_t ix,iy,iz,io;
    if(side){
        if(this->getHigh() < 0) rs_error("Domain<T>::copyFromboundary: Invalid side.");
    }else{
        if(this->getLow() < 0) rs_error("Domain<T>::copyFromboundary: Invalid side.");
    }
    size_t pad=0, n0=0, n1=0, n2=0, nwrk=0;

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
                            wrk0[IDWRK(io,iy,iz)] = array[IDARRAY(n0-2*pad + io, iy, iz)];
                        }
                    }
                }
            }else{
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk0[IDWRK(io,iy,iz)] = array[IDARRAY(io + pad, iy, iz)];
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
                            wrk0[IDWRK(io,ix,iz)] = array[IDARRAY(ix, n1-2*pad + io, iz)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                            wrk0[IDWRK(io,ix,iz)] = array[IDARRAY(ix, io + pad, iz)];
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
                            wrk0[IDWRK(io,ix,iy)] = array[IDARRAY(ix, iy, n2-2*pad + io)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            wrk0[IDWRK(io,ix,iy)] = array[IDARRAY(ix, iy, io + pad)];
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
    size_t ix,iy,iz,io;
    if(side){
        if(this->getHigh() < 0) rs_error("Domain<T>::copyToboundary: Invalid side.");
    }else{
        if(this->getLow() < 0) rs_error("Domain<T>::copyToboundary: Invalid side.");
    }
    size_t pad=0, n0=0, n1=0, n2=0, nwrk=0;

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
                           array[IDARRAY(n0-pad + io, iy, iz)] = wrk0[IDWRK(io,iy,iz)];
                        }
                    }
                }
            }else{
                for (io=0; io < pad; io++) {
                    for (iy=0; iy < n1; iy++) {
                        for (iz=0; iz < n2; iz++) {
                            array[IDARRAY(io, iy, iz)] = wrk0[IDWRK(io,iy,iz)];
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
                            array[IDARRAY(ix, n1-pad + io, iz)] = wrk0[IDWRK(io,ix,iz)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (io=0; io < pad; io++) {
                        for (iz=0; iz < n2; iz++) {
                          array[IDARRAY(ix, io, iz)] = wrk0[IDWRK(io,ix,iz)];
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
                            array[IDARRAY(ix, iy, n2-pad + io)] = wrk0[IDWRK(io,ix,iy)];
                        }
                    }
                }
            }else{
                for (ix=0; ix < n0; ix++) {
                    for (iy=0; iy < n1; iy++) {
                        for (io=0; io < pad; io++) {
                            array[IDARRAY(ix, iy, io)] = wrk0[IDWRK(io,ix,iy)];
                        }
                    }
                }
            }
            break;
        default:
            rs_error("Domain<T>::copyToboundary:Invalid dimension");
            break;
            break;
    }

}


template<typename T>
void Domain<T>::shareEdges(T *array){
    if(!mpiset) rs_error("Domain<T>::shareEdges: MPI not set in domain.");
    if(mpi->getDomainrank() % 2){
        if(this->getLow() >= 0){
            mpi->receiveEdges(wrk0, wrksize[0], this->getLow());
            this->copyToboundary(0, array);
        }
        if(this->getHigh() >= 0){
            mpi->receiveEdges(wrk0, wrksize[0], this->getHigh());
            this->copyToboundary(1, array);
        }
        if(this->getLow() >= 0){
            this->copyFromboundary(0, array);
            mpi->sendEdges(wrk0, wrksize[0], this->getLow());
        }
        if(this->getHigh() >= 0){
            this->copyFromboundary(1, array);
            mpi->sendEdges(wrk0, wrksize[0], this->getHigh());
        }
    }else{
        if(this->getHigh() >= 0){
            this->copyFromboundary(1, array);
            mpi->sendEdges(wrk0, wrksize[0], this->getHigh());
        }
        if(this->getLow() >= 0){
            this->copyFromboundary(0, array);
            mpi->sendEdges(wrk0, wrksize[0], this->getLow());
        }
        if(this->getHigh() >= 0){
            mpi->receiveEdges(wrk0, wrksize[0], this->getHigh());
            this->copyToboundary(1, array);
        }
        if(this->getLow() >= 0){
            mpi->receiveEdges(wrk0, wrksize[0], this->getLow());
            this->copyToboundary(0, array);
        }
    }
}


template<typename T>
Domain<T>::~Domain() {
    if(allocated[0]){
        free(wrk0);
    }

    if(allocated[1]){
        free(wrk1);
    }

    if(allocated[2]){
        free(wrk2);
    }

    if(allocated[3]){
        free(wrkedg0);
    }

    if(allocated[4]){
        free(wrkedg1);
    }

    if(allocated[5]){
        free(wrkedg2);
    }

    if(allocated[6]){
        free(wrkcrn);
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
