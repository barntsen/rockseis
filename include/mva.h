#ifndef MVA_H
#define MVA_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <time.h>
#include "geometry.h"
#include "utils.h"
#include "file.h"
#include "model.h"
#include "data.h"
#include "waves.h"
#include "der.h"
#include "snap.h"
#include "rtm.h"
#include "image.h"
#include "revolve.h"

#define MVA_OK 1
#define MVA_ERR 0

#define GMAP 1
#define SMAP 0

#define ki2D(i,j,k,l) ((l)*nhx*nz*nx + (k)*nx*nz + (j)*nx +(i))
#define km2D(i,j) ((j)*nx + (i))
#define ks2D(i,j) ((j)*nxs + (i))
#define kr2D(i,j) ((j)*nxr + (i))

#define ki3D(i,j,k,l,m,n) ((n)*nhy*nhx*nx*ny*nz + (m)*nhx*nx*ny*nz + (l)*nx*ny*nz + (k)*nx*ny + (j)*nx + (i))
#define km3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define ks3D(i,j,k) ((k)*nxs*nys + (j)*nxs + (i))
#define kr3D(i,j,k) ((k)*nxr*nyr + (j)*nxr + (i))

namespace rockseis {

// =============== ABSTRACT MVA CLASS =============== //
/** The abstract rtm class
 *
 */

template<typename T>
class Mva {
public:
    Mva();	///< Constructor
    Mva(int order, int snapinc);  ///< Constructor
    
    // Mva functions
    int getOrder() { return order; } ///< Get order of FD stencil
    int getSnapinc() { return snapinc; } ///< Get snap increment
    std::string getLogfile() { return logfile; } ///< Get name of logfile
    std::string getFwsnapfile() { return fwsnapfile; } ///< Sets checkpoint filename
    std::string getBwsnapfile() { return bwsnapfile; } ///< Sets checkpoint filename
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil
    void setSnapinc(int _snapinc) {snapinc = _snapinc;} ///< Set snap increment for recording snapshots
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setSnapmethod(rs_snapmethod val) { snapmethod = val; } ///< Sets choice of snapshot saving
    void setIncore(bool val) { incore = val; } ///< Sets optimal checkpoint incore flag
    void setNcheck(int val) { ncheck = val; } ///< Sets optimal checkpointing number of snaps 
    void setFwsnapfile(std::string file) { fwsnapfile = file; } ///< Sets checkpoint filename
    void setBwsnapfile(std::string file) { bwsnapfile = file; } ///< Sets checkpoint filename
    int getNcheck() { return ncheck; } ///< Gets the number of checkpoints for the optimal checkpointing scheme
    bool getIncore() { return incore; } ///< Gets the incore flag for the optimal checkpointing scheme

    void setMisfit(T val) { misfit = val; } ///< Sets data misfit value
    T getMisfit() { return misfit; }   ///< Gets misfit value
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

    ~Mva();	///< Destructor


private:
    int order;  ///< Order of the FD stencil 
    int snapinc;  ///< Snap interval
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    rs_snapmethod snapmethod; ///< Choice of checkpointing method
    T misfit; ///< Misfit value
    bool incore; ///< Incore flag for optimal checkpointing (No IO)
    int ncheck; ///< Number of checkpoints in optimal checkpointing
    std::string fwsnapfile;
    std::string bwsnapfile;
};

/** The 2D Acoustic Mva class
 *
 */
template<typename T>
class MvaAcoustic2D: public Mva<T> {
public:
    MvaAcoustic2D();					///< Constructor
    MvaAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Image2D<T>> pimage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T *vp, T* rho);
    void computeMisfit();
    void computeResiduals();

    ~MvaAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Image2D<T>> pimage;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataP;
    bool modelset;
    bool pimageset;
    bool vpgradset;
    bool sourceset;
    bool dataPset;
};


}
#endif //MVA_H
