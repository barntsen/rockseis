#ifndef DOMAIN_H
#define DOMAIN_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "utils.h"
#include "parallel.h"

#define DOMAIN_OK 1
#define DOMAIN_ERR 0

#define IDWRK(i,j,k) ((k)*pad*nwrk + (j)*pad + (i)) 
#define IDARRAY(i,j,k) ((k)*n0*n1 + (j)*n0 + (i)) 
#define IDDOM(i,j,k) ((k)*nd0*nd1 + (j)*nd0 + (i)) 
#define IDEDGE(i,j,k) ((k)*4 + (j)*2 + (i))
#define IDCRN(i,j,k) ((k)*4 + (j)*2 + (i))

namespace rockseis {

// =============== ABSTRACT DOMAIN CLASS =============== //
/** The abstract snap class
 *
 */
template<typename T>
class Domain {
public:
    Domain(); ///<Constructor
    virtual ~Domain();  ///< Destructor

    // Set functions
    void setLpad(const int _lpad) { lpad = _lpad; }     ///< Set lpad
    void setDim(const int _dim) { dim = _dim; } ///< Set dim
    void setLow(const int val, const int dim);   ///< Set rank of low side 
    void setLow(const int val) { setLow(val,0); }   ///< Set rank of low side 
    void setHigh(const int val, const int dim);   ///< Set rank of low side 
    void setHigh(const int val) { setHigh(val,0); }   ///< Set rank of low side 
    void setPadl(const int val, const int dim);   ///< Set rank of low side 
    void setPadl(const int val) { setPadl(val,0); }   ///< Set rank of low side 
    void setPadh(const int val, const int dim);   ///< Set rank of low side 
    void setPadh(const int val) { setPadh(val,0); }   ///< Set rank of low side 
    void setNx_orig(const int val) { nx_orig = val; }     ///< Set nx_orig
    void setNy_orig(const int val) { ny_orig = val; }     ///< Set ny_orig
    void setNz_orig(const int val) { nz_orig = val; }     ///< Set nz_orig
    void setNx_pad(const int val) { nx_pad = val; }     ///< Set nx_pad
    void setNy_pad(const int val) { ny_pad = val; }     ///< Set ny_pad
    void setNz_pad(const int val) { nz_pad = val; }     ///< Set nz_pad
    void setNd(const int _nd) { nd=_nd; }///< Set total number of domains
    void setNd0(const int _nd0) { nd0=_nd0; }///< Set number of domains in dimension 0
    void setNd1(const int _nd1) { nd1=_nd1; }///< Set number of domains in dimension 1
    void setNd2(const int _nd2) { nd2=_nd2; }///< Set number of domains in dimension 2
    void setD(const int _d) { d=_d; }///< Set domain number
    void setCorner(const int val, const int i) { if(i>=0 && i<8) corner[i] = val; else rs_error("Domain::setCorner: invalid corner index."); }
    void setEdge(const int val, const int i) { if(i>=0 && i<12) edge[i] = val; else rs_error("Domain::setEdge: invalid edge index."); }

    void setMpi(MPIdomaindecomp *_mpi) { mpi = _mpi; mpiset = true;}
    void copyFromboundary(const int dim, const bool side, const T *array); ///< Copy from array to wrk
    void copyFromedge(const int id, const int d0, const int d1, const T *array); ///< Get an edge from the array to send to neighour
    void copyFromcorner(const int d0, const int d1, const int d2, const T *array); ///< Get an edge from the array to send to neighour
    void copyToboundary(const int dim, const bool side, T *array); ///< Copy from array to wrk
    void copyToedge(const int id, const int d0, const int d1, T *array); ///< Get an edge from the array to send to neighour
    void copyTocorner(const int d0, const int d1, const int d2, T *array); ///< Get an edge from the array to send to neighour
    void copyFromboundary(const bool side, const T *array); ///< Copy from array to wrk
    void copyToboundary(const bool side, T *array); ///< Copy from wrk to array
    void shareEdges1D(T *array); ///< Share edges of array with neighbors
    void shareEdges3D(T *array); ///< Share edges of array with neighbors

    void setIx0(const int _ix0) { ix0 = _ix0; }     ///< Set ix0 (origin in global model)
    void setIy0(const int _iy0) { iy0 = _iy0; }     ///< Set iy0 (origin in global model)
    void setIz0(const int _iz0) { iz0 = _iz0; }     ///< Set iz0 (origin in global model)

    // Get functions
    int getNx_orig() { return nx_orig; }     ///< Get original model Nx 
    int getNy_orig() { return ny_orig ; }     ///< Get original model
    int getNz_orig() { return nz_orig ; }     ///< Get original model
    int getNx_pad() { return nx_pad; }     ///< Get Nx padded
    int getNy_pad() { return ny_pad ; }     ///< Get Ny padded
    int getNz_pad() { return nz_pad ; }     ///< Get Nz padded
    int getLpad() { return lpad; }              ///< Get lpad
    int getDim() { return dim; }                ///< Get dim
    int getLow(const int dim);                ///< Get low rank 
    int getLow() { return getLow(0); }                ///< Get low rank 
    int getHigh(const int dim);                ///< Get low rank 
    int getHigh() { return getHigh(0); }                ///< Get low rank 
    int getPadl(const int dim);                ///< Get padding on low side  
    int getPadl() { return getPadl(0); }                ///< Get padding on low side  
    int getPadh(const int dim);                ///< Get padding on high side  
    int getPadh() { return getPadh(0); }                ///< Get padding on high side  
    int getNd() { return nd; } ///< GetNd
    int getNd0() { return nd0; } ///< GetNd0
    int getNd1() { return nd1; } ///< GetNd1
    int getNd2() { return nd2; } ///< GetNd2
    int getD() { return d; } ///< Get d
    MPIdomaindecomp * getMpi() { return mpi; }
    int getIx0() { return ix0; } ///< Get ix0
    int getIy0() { return iy0; } ///< Get iy0
    int getIz0() { return iz0; } ///< Get iz0
    int getId(const int d, const int dim); ///< Get domain index in dimension dim 
    int getEdge(const int i) { if(i>=0 && i<12) return edge[i]; else {rs_error("Domain::getEdge: invalid corner index."); return -1; }}
    int getCorner(const int i) { if(i>=0 && i<8) return corner[i]; else {rs_error("Domain::getCorner: invalid corner index."); return -1; }}
    bool getStatus() { return status; }
    void setupDomain1D(const int nx, const int ny, const int nz, const int d, const int nd, const int order);
    void setupDomain3D(const int nx, const int ny, const int nz, const int d, const int nd0, const int nd1, const int nd2, const int order);

   private:
    int nx_orig, ny_orig, nz_orig;
    int nx_pad, ny_pad, nz_pad;
    int lpad; // Domain boundary padding
    int dim;
    Point3D<int> low,high;
    Point3D<int> padl, padh;
    int corner[8]; // 0 - lll; 1 - llh; 2 - lhl; 3 - lhh; 4 - hll; 5 - hlh; 6 - hhl; 7 - hhh; 
    int edge[12]; // 0-3 - xy; 4-7 - xz; 8-11 - yz; 
    int nd,nd0,nd1,nd2;
    int d;
    int ix0, iy0, iz0; // Global indexes for the origin of domain 
    bool status; // 1 for on, 0 for off. 
    bool allocated[7];
    T *wrk0; // Work array to be used to communicate boundaries
    T *wrk1; // Work array to be used to communicate boundaries
    T *wrk2; // Work array to be used to communicate boundaries
    T *wrkedg0; // Work array to be used to communicate edges
    T *wrkedg1; // Work array to be used to communicate edges
    T *wrkedg2; // Work array to be used to communicate edges
    T *wrkcrn; // Work array to be used to communicate corners
    size_t wrksize[3]; // Size of workarray
    size_t wrkedgsize[3]; // Size of edge workarray
    size_t wrkcrnsize; // Size of corner workarray
    MPIdomaindecomp *mpi; // Domain decomposition mpi class
    bool mpiset;
};

}
#endif //DOMAIN_H

