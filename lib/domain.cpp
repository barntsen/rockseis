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
            sndwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(sndwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating sndwrkside0 array");
            rcvwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(rcvwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating rcvwrkside0 array");
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
            sndwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(sndwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating sndwrkside0 array");
            rcvwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(rcvwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating rcvwrkside0 array");
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
            sndwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(sndwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating sndwrkside0 array");
            rcvwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
            if(rcvwrkside0 == NULL) rs_error("Domain<T>::setupDomain:Error allocating rcvwrkside0 array");
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
    if(nxdom < this->getLpad() && _nd0 > 1){
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
    if(nydom < this->getLpad() && _nd1 > 1){
       rs_error("Domain<T>::setupDomain:Number of grid cells in a domain is less than the order of the finite difference stencil. Decrease the number of domains or the stencil.");
    }
    nypad = nydom+this->getPadl(1)+this->getPadh(1);

    nzdom = std::ceil((T) nz/((T) _nd2));

    iz0 = getId(_d, 2)*nzdom - this->getPadl(2);
    // Treat special case of last domain which can be smaller than the other domains
    if (getId(_d, 2) == _nd2-1){ 
       nzdom = nz - nzdom*(_nd2-1);
    }

    if(nzdom < this->getLpad() && _nd2 > 1){
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
       sndwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
       if(sndwrkside0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkside0 array");
       rcvwrkside0 = (T *) calloc(wrksize[0], sizeof(T));
       if(rcvwrkside0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkside0 array");
       this->allocated[0] = true;
    }

    if(this->getPadl(1) > 0 || this->getPadh(1) > 0){
       wrksize[1] = nxwrk*nzwrk*this->getLpad();
       sndwrkside1 = (T *) calloc(wrksize[1], sizeof(T));
       if(sndwrkside1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkside1 array");
       rcvwrkside1 = (T *) calloc(wrksize[1], sizeof(T));
       if(rcvwrkside1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkside1 array");
       this->allocated[1] = true;
    }

    if(this->getPadl(2) > 0 || this->getPadh(2) > 0){
       wrksize[2] = nxwrk*nywrk*this->getLpad();
       sndwrkside2 = (T *) calloc(wrksize[2], sizeof(T));
       if(sndwrkside2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkside2 array");
       rcvwrkside2 = (T *) calloc(wrksize[2], sizeof(T));
       if(rcvwrkside2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkside2 array");
       this->allocated[2] = true;
    }

    int edgsum=0;
    for(int i=0; i<4; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
       wrkedgsize[0] = nzwrk*this->getLpad()*this->getLpad();
       sndwrkedg0 = (T *) calloc(wrkedgsize[0], sizeof(T));
       if(sndwrkedg0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkedg0 array");
       rcvwrkedg0 = (T *) calloc(wrkedgsize[0], sizeof(T));
       if(rcvwrkedg0 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkedg0 array");
       this->allocated[3] = true;
    }

    edgsum=0;
    for(int i=4; i<8; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
    wrkedgsize[1] = nywrk*this->getLpad()*this->getLpad();
       sndwrkedg1 = (T *) calloc(wrkedgsize[1], sizeof(T));
       if(sndwrkedg1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkedg1 array");
       rcvwrkedg1 = (T *) calloc(wrkedgsize[1], sizeof(T));
       if(rcvwrkedg1 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkedg1 array");
       this->allocated[4] = true;
    }

    edgsum=0;
    for(int i=8; i<12; i++) edgsum += this->getEdge(i);
    if(edgsum > -4){
    wrkedgsize[2] = nxwrk*this->getLpad()*this->getLpad();
       sndwrkedg2 = (T *) calloc(wrkedgsize[2], sizeof(T));
       if(sndwrkedg2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkedg2 array");
       rcvwrkedg2 = (T *) calloc(wrkedgsize[2], sizeof(T));
       if(rcvwrkedg2 == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkedg2 array");
       this->allocated[5] = true;
    }

    int crnsum=0;
    for(int i=0; i<8; i++) crnsum += this->getCorner(i);
    if(crnsum > -8){
       wrkcrnsize = this->getLpad()*this->getLpad();
       if(nywrk > 1) wrkcrnsize *= this->getLpad();
       sndwrkcrn = (T *) calloc(wrkcrnsize, sizeof(T));
       if(sndwrkcrn == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating sndwrkcrn array");
       rcvwrkcrn = (T *) calloc(wrkcrnsize, sizeof(T));
       if(rcvwrkcrn == NULL) rs_error("Domain<T>::setupDomain3D:Error allocating rcvwrkcrn array");
       this->allocated[6] = true;
    }

}

template<typename T>
void Domain<T>::copyFromside(const int dim, const bool side, const T *array){
   size_t ix,iy,iz,io;
   size_t padl0=0, padl1, padl2, padh0=0, padh1=0, padh2=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   padl0 = this->getPadl(0);
   padl1 = this->getPadl(1);
   padl2 = this->getPadl(2);
   padh0 = this->getPadh(0);
   padh1 = this->getPadh(1);
   padh2 = this->getPadh(2);

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();

   int pad = this->getLpad();

   switch(dim){
      case 0:
         nwrk0 = pad;
         nwrk1 = n1-padl1-padh1;
         if(side == 0){
            if(this->getLow(0)>=0){
               for (iz=0; iz < n2- padl2-padh2; iz++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (io=0; io < pad; io++) {
                        sndwrkside0[IDWRK(io,iy,iz)] = array[IDARRAY(io + pad, iy+padl1, iz+padl2)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(0)>=0){
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (io=0; io < pad; io++) {
                        sndwrkside0[IDWRK(io,iy,iz)] = array[IDARRAY(n0-2*pad + io, iy+padl1, iz+padl2)];
                     }
                  }
               }
            }
         }
         break;
      case 1:
         nwrk0 = n0-padl0-padh0;
         nwrk1 = pad;
         if(side == 0){
            if(this->getLow(1)>=0){
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (io=0; io < pad; io++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        sndwrkside1[IDWRK(ix,io,iz)] = array[IDARRAY(ix+padl0, io+pad, iz+padl2)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(1)>=0){
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (io=0; io < pad; io++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        sndwrkside1[IDWRK(ix,io,iz)] = array[IDARRAY(ix+padl0, n1-2*pad+io, iz+padl2)];
                     }
                  }
               }
            }
         }
         break;
      case 2:
         nwrk0 = n0-padl0-padh0;
         nwrk1 = n1-padl1-padh1;
         if(side == 0){
            if(this->getLow(2)>=0){
               for (io=0; io < pad; io++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        sndwrkside2[IDWRK(ix,iy,io)] = array[IDARRAY(ix+padl0, iy+padl1, io+pad)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(2)>=0){
               for (io=0; io < pad; io++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        sndwrkside2[IDWRK(ix,iy,io)] = array[IDARRAY(ix+padl0, iy+padl1, n2-2*pad+io)];
                     }
                  }
               }
            }
         }
         break;
      default:
         rs_error("Domain::copyFromside:Invalid dimension");
         break;
   }
}

template<typename T>
void Domain<T>::copyFromedge(const int id, const int d0, const int d1, const T *array){
   size_t ix,iy,iz,io0,io1;
   size_t padl0=0, padl1, padl2, padh0=0, padh1=0, padh2=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   padl0 = this->getPadl(0);
   padl1 = this->getPadl(1);
   padl2 = this->getPadl(2);
   padh0 = this->getPadh(0);
   padh1 = this->getPadh(1);
   padh2 = this->getPadh(2);

   int pad = this->getLpad();

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();
   int start0, start1;


   if(id<4){
      if(d0<0){
         start0=pad;
      }else{
         start0=n0 - 2*pad;
      }

      if(d1<0){
         start1=pad;
      }else{
         start1=n1 - 2*pad;
      }
      nwrk0 = pad;
      nwrk1 = pad;
      for (iz=0; iz < n2-padl2-padh2; iz++) {
         for (io1=0; io1 < pad; io1++) {
            for (io0=0; io0 < pad; io0++) {
               sndwrkedg0[IDWRK(io0,io1,iz)]=array[IDARRAY(start0+io0, start1+io1,iz+padl2)];
            }
         }
      }
   }

   if(id>=4 && id < 8){
      if(d0<0){
         start0=pad;
      }else{
         start0=n0 - 2*pad;
      }

      if(d1<0){
         start1=pad;
      }else{
         start1=n2 - 2*pad;
      }
      nwrk0 = pad;
      nwrk1 = n1-padl1-padh1;
      for (io1=0; io1 < pad; io1++) {
         for (iy=0; iy < n1-padl1-padh1; iy++) {
            for (io0=0; io0 < pad; io0++) {
               sndwrkedg1[IDWRK(io0,iy,io1)]=array[IDARRAY(start0+io0, iy+padl1, start1+io1)];
            }
         }
      }
   }

   if(id>=8 && id < 12){
      if(d0<0){
         start0=pad;
      }else{
         start0=n1 - 2*pad;
      }

      if(d1<0){
         start1=pad;
      }else{
         start1=n2 - 2*pad;
      }
      nwrk0 = n0-padl0-padh0;
      nwrk1 = pad;
      for (io1=0; io1 < pad; io1++) {
         for (io0=0; io0 < pad; io0++) {
            for (ix=0; ix < n0-padl0-padh0; ix++) {
               sndwrkedg2[IDWRK(ix,io0,io1)]=array[IDARRAY(ix+padl0, start0+io0,  start1+io1)];
            }
         }
      }
   }
}

template<typename T>
void Domain<T>::copyFromcorner(const int d0, const int d1, const int d2, const T *array){
   size_t io0,io1,io2;
   size_t n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   int pad = this->getLpad();

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();
   int start0, start1, start2;


   if(d0<0){
      start0=pad;
   }else{
      start0=n0 - 2*pad;
   }

   if(d1<0){
      start1=pad;
   }else{
      start1=n1 - 2*pad;
   }

   if(d2<0){
      start2=pad;
   }else{
      start2=n2 - 2*pad;
   }
   nwrk0 = pad;
   nwrk1 = pad;
   for (io2=0; io2 < pad; io2++) {
      for (io1=0; io1 < pad; io1++) {
         for (io0=0; io0 < pad; io0++) {
            sndwrkcrn[IDWRK(io0,io1,io2)]=array[IDARRAY(start0+io0, start1+io1, start2+io2)];
         }
      }
   }
}

template<typename T>
void Domain<T>::copyToside(const int dim, const bool side, T *array){
   size_t ix,iy,iz,io;
   size_t padl0=0, padl1, padl2, padh0=0, padh1=0, padh2=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   padl0 = this->getPadl(0);
   padl1 = this->getPadl(1);
   padl2 = this->getPadl(2);
   padh0 = this->getPadh(0);
   padh1 = this->getPadh(1);
   padh2 = this->getPadh(2);

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();

   int pad = this->getLpad();

   switch(dim){
      case 0:
         if(side == 0){
            if(this->getLow(0)>=0){
               nwrk0 = pad;
               nwrk1 = n1-padl1-padh1;
               for (iz=0; iz < n2- padl2-padh2; iz++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (io=0; io < pad; io++) {
                        array[IDARRAY(io, iy+padl1, iz+padl2)] = rcvwrkside0[IDWRK(io,iy,iz)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(0)>=0){
               nwrk0 = pad;
               nwrk1 = n1-padl1-padh1;
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (io=0; io < pad; io++) {
                        array[IDARRAY(n0-pad+io, iy+padl1, iz+padl2)] = rcvwrkside0[IDWRK(io,iy,iz)];
                     }
                  }
               }
            }
         }
         break;
      case 1:
         if(side == 0){
            if(this->getLow(1)>=0){
               nwrk0 = n0-padl0-padh0;
               nwrk1 = pad;
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (io=0; io < pad; io++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        array[IDARRAY(ix+padl0, io, iz+padl2)] = rcvwrkside1[IDWRK(ix,io,iz)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(1)>=0){
               nwrk0 = n0-padl0-padh0;
               nwrk1 = pad;
               for (iz=0; iz < n2-padl2-padh2; iz++) {
                  for (io=0; io < pad; io++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        array[IDARRAY(ix+padl0, n1-pad+io, iz+padl2)] =  rcvwrkside1[IDWRK(ix,io,iz)];
                     }
                  }
               }
            }
         }
         break;
      case 2:
         if(side == 0){
            if(this->getLow(2)>=0){
               nwrk0 = n0-padl0-padh0;
               nwrk1 = n1-padl1-padh1;
               for (io=0; io < pad; io++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        array[IDARRAY(ix+padl0, iy+padl1, io)] = rcvwrkside2[IDWRK(ix,iy,io)];
                     }
                  }
               }
            }
         }else{
            if(this->getHigh(2)>=0){
               nwrk0 = n0-padl0-padh0;
               nwrk1 = n1-padl1-padh1;
               for (io=0; io < pad; io++) {
                  for (iy=0; iy < n1-padl1-padh1; iy++) {
                     for (ix=0; ix < n0-padl0-padh0; ix++) {
                        array[IDARRAY(ix+padl0, iy+padl1, n2-pad+io)] = rcvwrkside2[IDWRK(ix,iy,io)];
                     }
                  }
               }
            }
         }
         break;
      default:
         rs_error("Domain::copyToside:Invalid dimension");
         break;
   }
}

template<typename T>
void Domain<T>::copyToedge(const int id, const int d0, const int d1, T *array){
   size_t ix,iy,iz,io0,io1;
   size_t padl0=0, padl1, padl2, padh0=0, padh1=0, padh2=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   padl0 = this->getPadl(0);
   padl1 = this->getPadl(1);
   padl2 = this->getPadl(2);
   padh0 = this->getPadh(0);
   padh1 = this->getPadh(1);
   padh2 = this->getPadh(2);

   int pad = this->getLpad();

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();
   int start0, start1;


   if(id<4){
      if(d0<0){
         start0=0;
      }else{
         start0=n0 - pad;
      }

      if(d1<0){
         start1=0;
      }else{
         start1=n1 - pad;
      }
      nwrk0 = pad;
      nwrk1 = pad;
      for (iz=0; iz < n2-padl2-padh2; iz++) {
         for (io1=0; io1 < pad; io1++) {
            for (io0=0; io0 < pad; io0++) {
               array[IDARRAY(start0+io0, start1+io1,iz+padl2)] = rcvwrkedg0[IDWRK(io0,io1,iz)];
            }
         }
      }
   }

   if(id>=4 && id < 8){
      if(d0<0){
         start0=0;
      }else{
         start0=n0 - pad;
      }

      if(d1<0){
         start1=0;
      }else{
         start1=n2 - pad;
      }
      nwrk0 = pad;
      nwrk1 = n1-padl1-padh1;
      for (io1=0; io1 < pad; io1++) {
         for (iy=0; iy < n1-padl1-padh1; iy++) {
            for (io0=0; io0 < pad; io0++) {
               array[IDARRAY(start0+io0, iy+padl1, start1+io1)] =  rcvwrkedg1[IDWRK(io0,iy,io1)];
            }
         }
      }
   }

   if(id>=8 && id < 12){
      if(d0<0){
         start0=0;
      }else{
         start0=n1 - pad;
      }

      if(d1<0){
         start1=0;
      }else{
         start1=n2 - pad;
      }
      nwrk0 = n0-padl0-padh0;
      nwrk1 = pad;
      for (io1=0; io1 < pad; io1++) {
         for (io0=0; io0 < pad; io0++) {
            for (ix=0; ix < n0-padl0-padh0; ix++) {
               array[IDARRAY(ix+padl0, start0+io0,  start1+io1)] = rcvwrkedg2[IDWRK(ix,io0,io1)];
            }
         }
      }
   }
}

template<typename T>
void Domain<T>::copyTocorner(const int d0, const int d1, const int d2, T *array){
   size_t io0,io1,io2;
   size_t n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   int pad = this->getLpad();

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();
   int start0, start1, start2;

   if(d0<0){
      start0=0;
   }else{
      start0=n0 - pad;
   }

   if(d1<0){
      start1=0;
   }else{
      start1=n1 - pad;
   }

   if(d2<0){
      start2=0;
   }else{
      start2=n2 - pad;
   }
   nwrk0 = pad;
   nwrk1 = pad;
   for (io2=0; io2 < pad; io2++) {
      for (io1=0; io1 < pad; io1++) {
         for (io0=0; io0 < pad; io0++) {
            array[IDARRAY(start0+io0, start1+io1, start2+io2)] = rcvwrkcrn[IDWRK(io0,io1,io2)];
         }
      }
   }
}

template<typename T>
void Domain<T>::copyFromside(const bool side, const T *array){
   size_t ix,iy,iz,io;
   if(side){
      if(this->getHigh() < 0) rs_error("Domain<T>::copyFromside: Invalid side.");
   }else{
      if(this->getLow() < 0) rs_error("Domain<T>::copyFromside: Invalid side.");
   }
   size_t pad=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();

   n2 = this->getNz_pad();
   pad = this->getLpad();

   switch(this->getDim()){
      case 0:
         nwrk0 = pad;
         nwrk1 = n1;
         if(side){
            for (iz=0; iz < n2; iz++) {
               for (iy=0; iy < n1; iy++) {
                  for (io=0; io < pad; io++) {
                     sndwrkside0[IDWRK(io,iy,iz)] = array[IDARRAY(n0-2*pad + io, iy, iz)];
                  }
               }
            }
         }else{
            for (iz=0; iz < n2; iz++) {
               for (iy=0; iy < n1; iy++) {
                  for (io=0; io < pad; io++) {
                     sndwrkside0[IDWRK(io,iy,iz)] = array[IDARRAY(io + pad, iy, iz)];
                  }
               }
            }
         }
         break;
      case 1:
         nwrk0 = n0;
         nwrk1 = pad;
         if(side){
            for (iz=0; iz < n2; iz++) {
               for (io=0; io < pad; io++) {
                  for (ix=0; ix < n0; ix++) {
                     sndwrkside0[IDWRK(ix,io,iz)] = array[IDARRAY(ix, n1-2*pad + io, iz)];
                  }
               }
            }
         }else{
            for (iz=0; iz < n2; iz++) {
               for (io=0; io < pad; io++) {
                  for (ix=0; ix < n0; ix++) {
                     sndwrkside0[IDWRK(ix,io,iz)] = array[IDARRAY(ix, io + pad, iz)];
                  }
               }
            }
         }
         break;
      case 2:
         nwrk0 = n0;
         nwrk1 = n1;
         if(side){
            for (io=0; io < pad; io++) {
               for (iy=0; iy < n1; iy++) {
                  for (ix=0; ix < n0; ix++) {
                     sndwrkside0[IDWRK(ix,iy,io)] = array[IDARRAY(ix, iy, n2-2*pad + io)];
                  }
               }
            }
         }else{
            for (io=0; io < pad; io++) {
               for (iy=0; iy < n1; iy++) {
                  for (ix=0; ix < n0; ix++) {
                     sndwrkside0[IDWRK(ix,iy,io)] = array[IDARRAY(ix, iy, io + pad)];
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
void Domain<T>::copyToside(const bool side, T *array){
   size_t ix,iy,iz,io;
   if(side){
      if(this->getHigh() < 0) rs_error("Domain<T>::copyToside: Invalid side.");
   }else{
      if(this->getLow() < 0) rs_error("Domain<T>::copyToside: Invalid side.");
   }
   size_t pad=0, n0=0, n1=0, n2=0, nwrk0=0, nwrk1=0;

   n0 = this->getNx_pad();
   n1 = this->getNy_pad();
   n2 = this->getNz_pad();
   pad = this->getLpad();

   switch(this->getDim()){
      case 0:
         nwrk0 = pad;
         nwrk1 = n1;
         if(side){
            for (iz=0; iz < n2; iz++) {
               for (iy=0; iy < n1; iy++) {
                  for (io=0; io < pad; io++) {
                     array[IDARRAY(n0-pad + io, iy, iz)] = rcvwrkside0[IDWRK(io,iy,iz)];
                  }
               }
            }
         }else{
            for (iz=0; iz < n2; iz++) {
               for (iy=0; iy < n1; iy++) {
                  for (io=0; io < pad; io++) {
                     array[IDARRAY(io, iy, iz)] = rcvwrkside0[IDWRK(io,iy,iz)];
                  }
               }
            }
         }
         break;
      case 1:
         nwrk0 = n0;
         nwrk1 = pad;
         if(side){
            for (iz=0; iz < n2; iz++) {
               for (io=0; io < pad; io++) {
                  for (ix=0; ix < n0; ix++) {
                     array[IDARRAY(ix, n1-pad + io, iz)] = rcvwrkside0[IDWRK(ix,io,iz)];
                  }
               }
            }
         }else{
            for (iz=0; iz < n2; iz++) {
               for (io=0; io < pad; io++) {
                  for (ix=0; ix < n0; ix++) {
                     array[IDARRAY(ix, io, iz)] = rcvwrkside0[IDWRK(ix,io,iz)];
                  }
               }
            }
         }
         break;
      case 2:
         nwrk0 = n0;
         nwrk1 = n1;
         if(side){
            for (io=0; io < pad; io++) {
               for (iy=0; iy < n1; iy++) {
                  for (ix=0; ix < n0; ix++) {
                     array[IDARRAY(ix, iy, n2-pad + io)] = rcvwrkside0[IDWRK(ix,iy,io)];
                  }
               }
            }
         }else{
            for (io=0; io < pad; io++) {
               for (iy=0; iy < n1; iy++) {
                  for (ix=0; ix < n0; ix++) {
                     array[IDARRAY(ix, iy, io)] = rcvwrkside0[IDWRK(ix,iy,io)];
                  }
               }
            }
         }
         break;
      default:
         rs_error("Domain<T>::copyToside:Invalid dimension");
         break;
         break;
   }

}

template<typename T>
void Domain<T>::shareEdges1D(T *array){
    if(!mpiset) rs_error("Domain<T>::shareEdges1D: MPI not set in domain.");

    if(mpi->getDomainrank() % 2){
        if(this->getLow() >= 0){
            this->copyFromside(0, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getLow());
            this->copyToside(0, array);
        }
        if(this->getHigh() >= 0){
            this->copyFromside(1, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getHigh());
            this->copyToside(1, array);
        }
    }else{
        if(this->getHigh() >= 0){
            this->copyFromside(1, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getHigh());
            this->copyToside(1, array);
        }
        if(this->getLow() >= 0){
            this->copyFromside(0, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getLow());
            this->copyToside(0, array);
        }
    }
}

template<typename T>
void Domain<T>::shareEdges3D(T *array){
    if(!mpiset) rs_error("Domain<T>::shareEdges3D: MPI not set in domain.");

    //Sides
    if(this->getId(mpi->getDomainrank(),0) % 2){
        if(this->getLow(0) >= 0){
            this->copyFromside(0, 0, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getLow(0));
            this->copyToside(0, 0, array);
        }
        if(this->getHigh(0) >= 0){
            this->copyFromside(0, 1, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getHigh(0));
            this->copyToside(0, 1, array);
        }
    }else{
        if(this->getHigh(0) >= 0){
            this->copyFromside(0, 1, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getHigh(0));
            this->copyToside(0, 1, array);
        }
        if(this->getLow(0) >= 0){
            this->copyFromside(0, 0, array);
            mpi->sendrecvEdges(sndwrkside0, rcvwrkside0, wrksize[0], this->getLow(0));
            this->copyToside(0, 0, array);
        }
    }

    if(this->getId(mpi->getDomainrank(),1) % 2){
       if(this->getLow(1) >= 0){
            this->copyFromside(1, 0, array);
            mpi->sendrecvEdges(sndwrkside1, rcvwrkside1, wrksize[1], this->getLow(1));
            this->copyToside(1, 0, array);
        }
        if(this->getHigh(1) >= 0){
            this->copyFromside(1, 1, array);
            mpi->sendrecvEdges(sndwrkside1,rcvwrkside1, wrksize[1], this->getHigh(1));
            this->copyToside(1, 1, array);
        }
    }else{
        if(this->getHigh(1) >= 0){
            this->copyFromside(1, 1, array);
            mpi->sendrecvEdges(sndwrkside1, rcvwrkside1, wrksize[1], this->getHigh(1));
            this->copyToside(1, 1, array);
        }
        if(this->getLow(1) >= 0){
            this->copyFromside(1, 0, array);
            mpi->sendrecvEdges(sndwrkside1, rcvwrkside1, wrksize[1], this->getLow(1));
            this->copyToside(1, 0, array);
        }
    }

    if(this->getId(mpi->getDomainrank(),2) % 2){
       if(this->getLow(2) >= 0){
          this->copyFromside(2, 0, array);
          mpi->sendrecvEdges(sndwrkside2, rcvwrkside2, wrksize[2], this->getLow(2));
          this->copyToside(2, 0, array);
       }
       if(this->getHigh(2) >= 0){
          this->copyFromside(2, 1, array);
          mpi->sendrecvEdges(sndwrkside2, rcvwrkside2, wrksize[2], this->getHigh(2));
          this->copyToside(2, 1, array);
       }
    }else{
       if(this->getHigh(2) >= 0){
          this->copyFromside(2, 1, array);
          mpi->sendrecvEdges(sndwrkside2, rcvwrkside2, wrksize[2], this->getHigh(2));
          this->copyToside(2, 1, array);
       }
       if(this->getLow(2) >= 0){
          this->copyFromside(2, 0, array);
          mpi->sendrecvEdges(sndwrkside2, rcvwrkside2, wrksize[2], this->getLow(2));
          this->copyToside(2, 0, array);
       }
    }

    //Edges
    int d0,d1,d2;
    int k;
    for(int i = 0; i < 2; i++){
       d0 = i*2 - 1;
       for(int j = 0; j < 2; j++){
          d1 = j*2 - 1;
          k = IDEDGE(i,j,0);
          if(this->getEdge(k) >= 0){
             if(k % 2){
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg0, rcvwrkedg0, wrkedgsize[0], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }else{
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg0, rcvwrkedg0, wrkedgsize[0], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }
          }
       }
    }

    for(int i = 0; i < 2; i++){
       d0 = i*2 - 1;
       for(int j = 0; j < 2; j++){
          d1 = j*2 - 1;
          k = IDEDGE(i,j,1);
          if(this->getEdge(k) >= 0){
             if(k % 2){
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg1, rcvwrkedg1, wrkedgsize[1], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }else{
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg1, rcvwrkedg1, wrkedgsize[1], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }
          }
       }
    }
    for(int i = 0; i < 2; i++){
       d0 = i*2 - 1;
       for(int j = 0; j < 2; j++){
          d1 = j*2 - 1;
          k = IDEDGE(i,j,2);
          if(this->getEdge(k) >= 0){
             if(k % 2){
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg2, rcvwrkedg2, wrkedgsize[2], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }else{
                this->copyFromedge(k, d0, d1, array);
                mpi->sendrecvEdges(sndwrkedg2, rcvwrkedg2, wrkedgsize[2], this->getEdge(k));
                this->copyToedge(k, d0, d1, array);
             }
          }
       }
    }

    //Corners
    int l;
    for(int i = 0; i < 2; i++){
       d0 = i*2 - 1;
       for(int j = 0; j < 2; j++){
          d1 = j*2 - 1;
          for(int k = 0; k < 2; k++){
             d2 = k*2 - 1;
             l = IDCRN(i,j,k);
             if(this->getCorner(l) >= 0){
                if(l % 2){
                   this->copyFromcorner(d0, d1, d2, array);
                   mpi->sendrecvEdges(sndwrkcrn, rcvwrkcrn, wrkcrnsize, this->getCorner(l));
                   this->copyTocorner(d0, d1, d2, array);
                }else{
                   this->copyFromcorner(d0, d1, d2, array);
                   mpi->sendrecvEdges(sndwrkcrn, rcvwrkcrn, wrkcrnsize, this->getCorner(l));
                   this->copyTocorner(d0, d1, d2, array);
                }
             }
          }
       }
    }
}

template<typename T>
Domain<T>::~Domain() {
    if(allocated[0]){
        free(sndwrkside0);
        free(rcvwrkside0);
    }

    if(allocated[1]){
        free(sndwrkside1);
        free(rcvwrkside1);
    }

    if(allocated[2]){
        free(sndwrkside2);
        free(rcvwrkside2);
    }

    if(allocated[3]){
        free(sndwrkedg0);
        free(rcvwrkedg0);
    }

    if(allocated[4]){
        free(sndwrkedg1);
        free(rcvwrkedg1);
    }

    if(allocated[5]){
        free(sndwrkedg2);
        free(rcvwrkedg2);
    }

    if(allocated[6]){
        free(sndwrkcrn);
        free(rcvwrkcrn);
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
