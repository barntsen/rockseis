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
#define kw2D(i,j) ((j)*nxw + (i))
#define ks2D(i,j) ((j)*nxs + (i))
#define kr2D(i,j) ((j)*nxr + (i))

#define ki3D(i,j,k,l,m,n) ((n)*nhy*nhx*nx*ny*nz + (m)*nhx*nx*ny*nz + (l)*nx*ny*nz + (k)*nx*ny + (j)*nx + (i))
#define km3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define kw3D(i,j,k) ((k)*nxw*nyw + (j)*nxw + (i))
#define ks3D(i,j,k) ((k)*nxs*nys + (j)*nxs + (i))
#define kr3D(i,j,k) ((k)*nxr*nyr + (j)*nxr + (i))
  
#define GAUSS(x,y) expf(-1.0*((x)*(x))/(2.0*(y)*(y)));

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
    std::string getSnapfile() { return snapfile; } ///< Sets checkpoint filename
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil
    void setSnapinc(int _snapinc) {snapinc = _snapinc;} ///< Set snap increment for recording snapshots
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setSnapmethod(rs_snapmethod val) { snapmethod = val; } ///< Sets choice of snapshot saving
    void setIncore(bool val) { incore = val; } ///< Sets optimal checkpoint incore flag
    void setNcheck(int val) { ncheck = val; } ///< Sets optimal checkpointing number of snaps 
    void setSnapfile(std::string file) { snapfile = file; } ///< Sets checkpoint filename
    int getNcheck() { return ncheck; } ///< Gets the number of checkpoints for the optimal checkpointing scheme
    bool getIncore() { return incore; } ///< Gets the incore flag for the optimal checkpointing scheme

    void setMisfit(T val) { misfit = val; } ///< Sets data misfit value
    T getMisfit() { return misfit; }   ///< Gets misfit value

    rs_wemvamisfit getMisfit_type() { return misfit_type; }  ///< Gets misfit type
    void setMisfit_type(rs_wemvamisfit type) { misfit_type = type;} ///< Sets wemva misfit type
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
    rs_wemvamisfit misfit_type; ///< Misfit type can be either difference semblance or stack-power or a combination of both
    rs_snapmethod snapmethod; ///< Choice of checkpointing method
    T misfit; ///< Misfit value
    bool incore; ///< Incore flag for optimal checkpointing (No IO)
    int ncheck; ///< Number of checkpoints in optimal checkpointing
    std::string snapfile;
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

    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T *vp, T* rho, T* adjsrc);
    void calcAdjointsource(T *adjsrc_fw, T* wsp, int pads, T *adjsrc_bw, T* wrp, int padr);
    void insertAdjointsource(std::shared_ptr<WavesAcoustic2D<T>> waves_fw, T* adjsrc_fw, std::shared_ptr<WavesAcoustic2D<T>> waves_bw, T* adjsrc_bw, T *L);
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

/** The 2D Elastic PP Mva class
 *
 */
template<typename T>
class PPmvaElastic2D: public Mva<T> {
public:
    PPmvaElastic2D();					///< Constructor
    PPmvaElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Image2D<T>> pimage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataUx, std::shared_ptr<Data2D<T>> dataUz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setDataUx(std::shared_ptr<Data2D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
    void setDataUz(std::shared_ptr<Data2D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    void crossCorr(T *wsx, T *wsz, int pads, std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, T *adjsrc);
    void calcAdjointsource(T *adjsrc_fw, T *wsx, T *wsz, int pads, T *adjsrc_bw, T* wrx, T* wrz, int padr, std::shared_ptr<ModelElastic2D<T>> model);
    void insertAdjointsource(std::shared_ptr<WavesElastic2D_DS<T>> waves_fw, T* adjsrc_fw, std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, T* adjsrc_bw, std::shared_ptr<ModelElastic2D<T>> model);
    ~PPmvaElastic2D();	///< Destructor

private:
    std::shared_ptr<ModelElastic2D<T>> model;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Image2D<T>> pimage;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataUx;
    std::shared_ptr<Data2D<T>> dataUz;
    bool modelset;
    bool pimageset;
    bool vpgradset;
    bool sourceset;
    bool dataUxset, dataUzset;
};

/** The 2D Elastic PS Mva class
 *
 */
//template<typename T>
//class PSmvaElastic2D: public Mva<T> {
//public:
//    PSmvaElastic2D();					///< Constructor
//    PSmvaElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Image2D<T>> simage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataUx, std::shared_ptr<Data2D<T>> dataUz, int order, int snapinc);					///< Constructor 
//    int run(); ///< Runs rtm with full snapshoting
//    int run_optimal(); ///< Runs rtm with optimal checkpointing
//    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
//    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
//    void setVsgrad(std::shared_ptr<Image2D<T>> _vsgrad) { vsgrad = _vsgrad; vsgradset = true; }
//    void setDataUx(std::shared_ptr<Data2D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
//    void setDataUz(std::shared_ptr<Data2D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
//    bool checkStability(); ///< Check stability of finite difference modelling
//
//    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T *vp, T* rho, T* adjsrc);
//    void calcAdjointsource(T *adjsrc_fw, T* wsp, int pads, T *adjsrc_bw, T* wrp, int padr);
//    void insertAdjointsource(std::shared_ptr<WavesElastic2D<T>> waves_fw, T* adjsrc_fw, std::shared_ptr<WavesElastic2D<T>> waves_bw, T* adjsrc_bw, T *L);
//    ~PSmvaElastic2D();	///< Destructor
//
//private:
//    std::shared_ptr<ModelElastic2D<T>> model;
//    std::shared_ptr<Image2D<T>> vsgrad;
//    std::shared_ptr<Image2D<T>> simage;
//    std::shared_ptr<Data2D<T>> source;
//    std::shared_ptr<Data2D<T>> dataUx;
//    std::shared_ptr<Data2D<T>> dataUz;
//    bool modelset;
//    bool simageset;
//    bool vsgradset;
//    bool sourceset;
//    bool dataUxset, dataUzset;
//};

}
#endif //MVA_H
