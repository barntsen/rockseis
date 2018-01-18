#ifndef RTM_H
#define RTM_H

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
#include "image.h"
#include "revolve.h"

#define RTM_OK 1
#define RTM_ERR 0

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

// =============== ENUMS =============== //
typedef enum {FULL, OPTIMAL, EDGES} rs_snapmethod; ///< Snapshot saving method

// =============== ABSTRACT RTM CLASS =============== //
/** The abstract rtm class
 *
 */

template<typename T>
class Rtm {
public:
    Rtm();	///< Constructor
    Rtm(int order, int snapinc);  ///< Constructor
    
    // Rtm functions
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
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

    ~Rtm();	///< Destructor


private:
    int order;  ///< Order of the FD stencil 
    int snapinc;  ///< Snap interval
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    rs_snapmethod snapmethod; ///< Choice of checkpointing method
    bool incore; ///< Incore flag for optimal checkpointing (No IO)
    int ncheck; ///< Number of checkpoints in optimal checkpointing
    std::string snapfile;
};

/** The 2D Acoustic Rtm class
 *
 */
template<typename T>
class RtmAcoustic2D: public Rtm<T> {
public:
    RtmAcoustic2D();					///< Constructor
    RtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Image2D<T>> pimage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_edge(); ///< Runs rtm with edge boundary saving
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data2D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    void crossCorr(T *ws, int pads, T* wr, int padr);

    ~RtmAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Image2D<T>> pimage;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataP;
    std::shared_ptr<Data2D<T>> dataAz;
    bool modelset;
    bool pimageset;
    bool sourceset;
    bool dataPset, dataAxset, dataAzset;
};

/** The 3D Acoustic Rtm class
 *
 */
template<typename T>
class RtmAcoustic3D: public Rtm<T> {
public:
    RtmAcoustic3D();					///< Constructor
    RtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Image3D<T>> pimage, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setDataP(std::shared_ptr<Data3D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data3D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling
    void crossCorr(T *ws, int pads, T* wr, int padr);

    ~RtmAcoustic3D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic3D<T>> model;
    std::shared_ptr<Image3D<T>> pimage;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataP;
    std::shared_ptr<Data3D<T>> dataAz;
    bool modelset;
    bool pimageset;
    bool sourceset;
    bool dataPset, dataAxset, dataAyset, dataAzset;
};

/** The 2D Elastic Rtm class
 *
 */
template<typename T>
class RtmElastic2D: public Rtm<T> {
public:
    RtmElastic2D();					///< Constructor
    RtmElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataUx, std::shared_ptr<Data2D<T>> dataUz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setDataUx(std::shared_ptr<Data2D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
    void setDataUz(std::shared_ptr<Data2D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
    void setPimage(std::shared_ptr<Image2D<T>> _pimage) { pimage = _pimage; pimageset = true; }
    void setSimage(std::shared_ptr<Image2D<T>> _simage) { simage = _simage; simageset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling
    void crossCorr(T *wsx, T *wsz, int pads, T* wrx, T* wrz, int padr, T* Vp, T* Vs, T* Rho);

    ~RtmElastic2D();	///< Destructor

private:
    std::shared_ptr<ModelElastic2D<T>> model;
    std::shared_ptr<Image2D<T>> pimage;
    std::shared_ptr<Image2D<T>> simage;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataUx;
    std::shared_ptr<Data2D<T>> dataUz;
    bool modelset;
    bool pimageset;
    bool simageset;
    bool sourceset;
    bool dataUxset, dataUzset;
};

/** The 3D Elastic Rtm class
 *
 */
template<typename T>
class RtmElastic3D: public Rtm<T> {
public:
    RtmElastic3D();					///< Constructor
    RtmElastic3D(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataVx, std::shared_ptr<Data3D<T>> dataVy, std::shared_ptr<Data3D<T>> dataVz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setDataUx(std::shared_ptr<Data3D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
    void setDataUy(std::shared_ptr<Data3D<T>> _dataUy) { dataUy = _dataUy; dataUyset = true; }
    void setDataUz(std::shared_ptr<Data3D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
    void setPimage(std::shared_ptr<Image3D<T>> _pimage) { pimage = _pimage; pimageset = true; }
    void setSimage(std::shared_ptr<Image3D<T>> _simage) { simage = _simage; simageset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling
    void crossCorr(T *wsx, T *wsy, T *wsz, int pads, T* wrx, T *wry, T* wrz, int padr, T* Vp, T* Vs, T* Rho);


    ~RtmElastic3D();	///< Destructor

private:
    std::shared_ptr<ModelElastic3D<T>> model;
    std::shared_ptr<Image3D<T>> pimage;
    std::shared_ptr<Image3D<T>> simage;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataUx;
    std::shared_ptr<Data3D<T>> dataUy;
    std::shared_ptr<Data3D<T>> dataUz;
    bool modelset;
    bool pimageset;
    bool simageset;
    bool sourceset;
    bool dataUxset, dataUyset, dataUzset;
};


}
#endif //RTM_H
