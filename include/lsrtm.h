#ifndef LSRTM_H
#define LSRTM_H

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

#define LSRTM_OK 1
#define LSRTM_ERR 0

#define GMAP 1
#define SMAP 0

#define PNOISE_ACOUSTIC 0.01
#define PNOISE_ELASTIC 0.01
#define STDEV 0.01
#define HMAX 500
#define PCLIP 97

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define CLOSETO(x, y) ((ABS((x) - (y)) <= FLT_EPSILON*ABS(y)) ? 1 : 0)

#define kim2D(i,j,k,l) ((l)*nhx*nz*nx + (k)*nx*nz + (j)*nx +(i))
#define ki2D(i,j) ((j)*nx +(i))
#define kw2D(i,j) ((j)*nxw + (i))
#define km2D(i,j) ((j)*nx + (i))
#define ks2D(i,j) ((j)*nxs + (i))
#define kr2D(i,j) ((j)*nxr + (i))

#define kim3D(i,j,k,l,m,n) ((n)*nhy*nhx*nx*ny*nz + (m)*nhx*nx*ny*nz + (l)*nx*ny*nz + (k)*nx*ny + (j)*nx + (i))
#define ki3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define km3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define kw3D(i,j,k) ((k)*nxw*nyw + (j)*nxw + (i))
#define ks3D(i,j,k) ((k)*nxs*nys + (j)*nxs + (i))
#define kr3D(i,j,k) ((k)*nxr*nyr + (j)*nxr + (i))

#define kwav(i,j) ((j)*nt + (i))

namespace rockseis {

// =============== ENUMS =============== //
typedef enum {FULL, OPTIMAL, EDGES} rs_snapmethod; ///< Snapshot saving method


// =============== ABSTRACT LSRTM CLASS =============== //
/** The abstract fwi class
 *
 */

template<typename T>
class Lsrtm {
public:
    Lsrtm();	///< Constructor
    Lsrtm(int order, int snapinc);  ///< Constructor
    ~Lsrtm();	///< Destructor
    
    // Lsrtm functions
    int getOrder() { return order; } ///< Get order of FD stencil
    int getSnapinc() { return snapinc; } ///< Get snap increment
    bool getFilter() { return filter; } ///< Get filter status
    T getFreqs(int i) { if(i >=0 && i < 4) return freqs[i]; else return -1;} ///< Gets filter corner
    std::string getLogfile() { return logfile; } ///< Get name of logfile
    std::string getSnapfile() { return snapfile; } ///< Sets checkpoint filename
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil
    void setSnapinc(int _snapinc) {snapinc = _snapinc;} ///< Set snap increment for recording snapshots
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setSnapmethod(rs_snapmethod val) { snapmethod = val; } ///< Sets choice of snapshot saving
    void setIncore(bool val) { incore = val; } ///< Sets optimal checkpoint incore flag
    void setNcheck(int val) { ncheck = val; } ///< Sets optimal checkpointing number of snaps 
    void setSnapfile(std::string file) { snapfile = file; } ///< Sets checkpoint filename
    void setFilter(bool val) { filter = val; } ///< Sets filter on or off
    void setFreqs(T freq, int i) { if(i >=0 && i < 4) freqs[i] = freq; } ///< Sets filter corner
    void setAllfreqs(T *freqptr) { freqs[0] = freqptr[0]; freqs[1] = freqptr[1]; freqs[2] = freqptr[2]; freqs[3] = freqptr[3]; }///< Sets all filter corner
    void setMisfit(T val) { misfit = val; } ///< Sets data misfit value
    T getMisfit() { return misfit; }   ///< Gets misfit value
    int getNcheck() { return ncheck; } ///< Gets the number of checkpoints for the optimal checkpointing scheme
    bool getIncore() { return incore; } ///< Gets the incore flag for the optimal checkpointing scheme
    rs_fwimisfit getMisfit_type() { return misfit_type; }  ///< Gets misfit type
    void setMisfit_type(rs_fwimisfit type) { misfit_type = type;} ///< Sets fwi misfit type
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);
    void setNoreverse(bool val) { noreverse = val; }
    bool getNoreverse() { return noreverse; }
    void stoep(int n, T r[], T g[], T f[], T a[]);
    void convolve(int lx, int ifx, T *x, int ly, int ify, T *y, int lz, int ifz, T *z);
    void xcor(int lx, int ifx, T *x,int ly, int ify, T *y, int lz, int ifz, T *z);
    T gauss(int it, T stdev);
    T linear(int it, T stdev);
    //void apply_filter(T *data, int nt, T dt);

private:
    int order;  ///< Order of the FD stencil 
    int snapinc;  ///< Snap interval
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    rs_snapmethod snapmethod; ///< Choice of checkpointing method
    bool incore; ///< Incore flag for optimal checkpointing (No IO)
    rs_fwimisfit misfit_type; ///< Misfit type can be either difference or correlation
    int ncheck; ///< Number of checkpoints in optimal checkpointing
    std::string snapfile;
    T misfit; ///< Misfit value
    bool noreverse; ///< Do only forward loop flag
    bool filter; ///< Whether or not to filter residuals in Adaptive misfit 
    T freqs[4];  ///< Four corner filter  
};

/** The 2D Acoustic Lsrtm class
 *
 */
template<typename T>
class LsrtmAcoustic2D: public Lsrtm<T> {
public:
    LsrtmAcoustic2D();					///< Constructor
    LsrtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Image2D<T>> pimage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_edge(); ///< Runs fwi with edge boundary saving
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setSrcilum(std::shared_ptr<Image2D<T>> _srcilum) { srcilum = _srcilum; srcilumset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data2D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }
    void setDatamodP(std::shared_ptr<Data2D<T>> _datamodP) { datamodP = _datamodP; datamodPset = true; }
    void setDatamodAz(std::shared_ptr<Data2D<T>> _datamodAz) { datamodAz = _datamodAz; datamodAzset = true; }
    void setDataresP(std::shared_ptr<Data2D<T>> _dataresP) { dataresP = _dataresP; dataresPset = true; }
    void setDataresAz(std::shared_ptr<Data2D<T>> _dataresAz) { dataresAz = _dataresAz; dataresAzset = true; }
    void setDataweight(std::shared_ptr<Data2D<T>> _dataweight) { dataweight = _dataweight; dataweightset = true; }

    void calcAdjointsource(T *adjsrc_bw, T* wrp, int padr);
    void insertAdjointsource(std::shared_ptr<WavesAcoustic2D<T>> waves_fw, T* adjsrc_fw, T *L);
    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T *vp, T* rho);

    void computeMisfit();
    void computeResiduals();
    ~LsrtmAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Image2D<T>> srcilum;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Image2D<T>> pimage;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataP;
    std::shared_ptr<Data2D<T>> dataAz;
    std::shared_ptr<Data2D<T>> datamodP;
    std::shared_ptr<Data2D<T>> datamodAz;
    std::shared_ptr<Data2D<T>> dataresP;
    std::shared_ptr<Data2D<T>> dataresAz;
    std::shared_ptr<Data2D<T>> dataweight;
    bool modelset;
    bool srcilumset;
    bool vpgradset;
    bool pimageset;
    bool sourceset;
    bool dataPset, dataAzset;
    bool datamodPset, datamodAzset;
    bool dataresPset, dataresAzset;
    bool dataweightset;
};

/** The 3D Acoustic Lsrtm class
 *
 */
template<typename T>
class LsrtmAcoustic3D: public Lsrtm<T> {
public:
    LsrtmAcoustic3D();					///< Constructor
    LsrtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setWavgrad(std::shared_ptr<Data3D<T>> _wavgrad) { wavgrad = _wavgrad; wavgradset = true; }
    void setVpgrad(std::shared_ptr<Image3D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setRhograd(std::shared_ptr<Image3D<T>> _rhograd) { rhograd = _rhograd; rhogradset = true; }
    void setDataP(std::shared_ptr<Data3D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data3D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }
    void setDatamodP(std::shared_ptr<Data3D<T>> _datamodP) { datamodP = _datamodP; datamodPset = true; }
    void setDatamodAz(std::shared_ptr<Data3D<T>> _datamodAz) { datamodAz = _datamodAz; datamodAzset = true; }
    void setDataresP(std::shared_ptr<Data3D<T>> _dataresP) { dataresP = _dataresP; dataresPset = true; }
    void setDataresAz(std::shared_ptr<Data3D<T>> _dataresAz) { dataresAz = _dataresAz; dataresAzset = true; }
    void setDataweight(std::shared_ptr<Data3D<T>> _dataweight) { dataweight = _dataweight; dataweightset = true; }
    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wry, T* wrz, int padr, T *vp, T* rho);
    void computeMisfit();
    void computeResiduals();

    ~LsrtmAcoustic3D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic3D<T>> model;
    std::shared_ptr<Image3D<T>> vpgrad;
    std::shared_ptr<Image3D<T>> rhograd;
    std::shared_ptr<Data3D<T>> wavgrad;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataP;
    std::shared_ptr<Data3D<T>> dataAz;
    std::shared_ptr<Data3D<T>> datamodP;
    std::shared_ptr<Data3D<T>> datamodAz;
    std::shared_ptr<Data3D<T>> dataresP;
    std::shared_ptr<Data3D<T>> dataresAz;
    std::shared_ptr<Data3D<T>> dataweight;
    bool modelset;
    bool vpgradset;
    bool rhogradset;
    bool wavgradset;
    bool sourceset;
    bool dataPset, dataAzset;
    bool datamodPset, datamodAzset;
    bool dataresPset, dataresAzset;
    bool dataweightset;
};

/** The 2D Elastic Lsrtm class
 *
 */
template<typename T>
class LsrtmElastic2D: public Lsrtm<T> {
public:
    LsrtmElastic2D();					///< Constructor
    LsrtmElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataUx, std::shared_ptr<Data2D<T>> dataUz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setDataUx(std::shared_ptr<Data2D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
    void setDataUz(std::shared_ptr<Data2D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
    void setDatamodUx(std::shared_ptr<Data2D<T>> _datamodUx) { datamodUx = _datamodUx; datamodUxset = true; }
    void setDatamodUz(std::shared_ptr<Data2D<T>> _datamodUz) { datamodUz = _datamodUz; datamodUzset = true; }
    void setDataresUx(std::shared_ptr<Data2D<T>> _dataresUx) { dataresUx = _dataresUx; dataresUxset = true; }
    void setDataresUz(std::shared_ptr<Data2D<T>> _dataresUz) { dataresUz = _dataresUz; dataresUzset = true; }
    void setWavgrad(std::shared_ptr<Data2D<T>> _wavgrad) { wavgrad = _wavgrad; wavgradset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setVsgrad(std::shared_ptr<Image2D<T>> _vsgrad) { vsgrad = _vsgrad; vsgradset = true; }
    void setRhograd(std::shared_ptr<Image2D<T>> _rhograd) { rhograd = _rhograd; rhogradset = true; }
    void setDataweightx(std::shared_ptr<Data2D<T>> _dataweightx) { dataweightx = _dataweightx; dataweightxset = true; }
    void setDataweightz(std::shared_ptr<Data2D<T>> _dataweightz) { dataweightz = _dataweightz; dataweightzset = true; }
    void crossCorr(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it);
    void crossCorr2(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it);
    void scaleGrad(std::shared_ptr<ModelElastic2D<T>> model);
    void computeMisfit();
    void computeResiduals();
    void computeResiduals2();

    ~LsrtmElastic2D();	///< Destructor

private:
    std::shared_ptr<ModelElastic2D<T>> model;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Image2D<T>> vsgrad;
    std::shared_ptr<Image2D<T>> rhograd;
    std::shared_ptr<Data2D<T>> wavgrad;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataUx;
    std::shared_ptr<Data2D<T>> dataUz;
    std::shared_ptr<Data2D<T>> datamodUx;
    std::shared_ptr<Data2D<T>> datamodUz;
    std::shared_ptr<Data2D<T>> dataresUx;
    std::shared_ptr<Data2D<T>> dataresUz;
    std::shared_ptr<Data2D<T>> dataweightx;
    std::shared_ptr<Data2D<T>> dataweightz;
    bool dataweightxset;
    bool dataweightzset;
    bool modelset;
    bool vpgradset;
    bool vsgradset;
    bool rhogradset;
    bool wavgradset;
    bool sourceset;
    bool dataUxset, dataUzset;
    bool datamodUxset, datamodUzset;
    bool dataresUxset, dataresUzset;
};

/** The 3D Elastic Lsrtm class
 *
 */
template<typename T>
class LsrtmElastic3D: public Lsrtm<T> {
public:
    LsrtmElastic3D();					///< Constructor
    LsrtmElastic3D(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataVx, std::shared_ptr<Data3D<T>> dataVy, std::shared_ptr<Data3D<T>> dataVz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setDataUx(std::shared_ptr<Data3D<T>> _dataUx) { dataUx = _dataUx; dataUxset = true; }
    void setDataUy(std::shared_ptr<Data3D<T>> _dataUy) { dataUy = _dataUy; dataUyset = true; }
    void setDataUz(std::shared_ptr<Data3D<T>> _dataUz) { dataUz = _dataUz; dataUzset = true; }
    void setDatamodUx(std::shared_ptr<Data3D<T>> _datamodUx) { datamodUx = _datamodUx; datamodUxset = true; }
    void setDatamodUy(std::shared_ptr<Data3D<T>> _datamodUy) { datamodUy = _datamodUy; datamodUyset = true; }
    void setDatamodUz(std::shared_ptr<Data3D<T>> _datamodUz) { datamodUz = _datamodUz; datamodUzset = true; }
    void setDataresUx(std::shared_ptr<Data3D<T>> _dataresUx) { dataresUx = _dataresUx; dataresUxset = true; }
    void setDataresUy(std::shared_ptr<Data3D<T>> _dataresUy) { dataresUy = _dataresUy; dataresUyset = true; }
    void setDataresUz(std::shared_ptr<Data3D<T>> _dataresUz) { dataresUz = _dataresUz; dataresUzset = true; }

    void setWavgrad(std::shared_ptr<Data3D<T>> _wavgrad) { wavgrad = _wavgrad; wavgradset = true; }
    void setVpgrad(std::shared_ptr<Image3D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setVsgrad(std::shared_ptr<Image3D<T>> _vsgrad) { vsgrad = _vsgrad; vsgradset = true; }
    void setRhograd(std::shared_ptr<Image3D<T>> _rhograd) { rhograd = _rhograd; rhogradset = true; }
    void setDataweightx(std::shared_ptr<Data3D<T>> _dataweightx) { dataweightx = _dataweightx; dataweightxset = true; }
    void setDataweighty(std::shared_ptr<Data3D<T>> _dataweighty) { dataweighty = _dataweighty; dataweightyset = true; }
    void setDataweightz(std::shared_ptr<Data3D<T>> _dataweightz) { dataweightz = _dataweightz; dataweightzset = true; }
    void crossCorr(T *wsx, T*wsy, T *wsz, int pads, std::shared_ptr<WavesElastic3D_DS<T>> waves_bw, std::shared_ptr<ModelElastic3D<T>> model, int it);
    void scaleGrad(std::shared_ptr<ModelElastic3D<T>> model);
    void computeMisfit();
    void computeResiduals();

    ~LsrtmElastic3D();	///< Destructor

private:
    std::shared_ptr<ModelElastic3D<T>> model;
    std::shared_ptr<Image3D<T>> vpgrad;
    std::shared_ptr<Image3D<T>> vsgrad;
    std::shared_ptr<Image3D<T>> rhograd;
    std::shared_ptr<Data3D<T>> wavgrad;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataUx;
    std::shared_ptr<Data3D<T>> datamodUx;
    std::shared_ptr<Data3D<T>> dataresUx;
    std::shared_ptr<Data3D<T>> dataUy;
    std::shared_ptr<Data3D<T>> datamodUy;
    std::shared_ptr<Data3D<T>> dataresUy;
    std::shared_ptr<Data3D<T>> dataUz;
    std::shared_ptr<Data3D<T>> datamodUz;
    std::shared_ptr<Data3D<T>> dataresUz;
    std::shared_ptr<Data3D<T>> dataweightx;
    std::shared_ptr<Data3D<T>> dataweighty;
    std::shared_ptr<Data3D<T>> dataweightz;
    bool modelset;
    bool vpgradset;
    bool vsgradset;
    bool rhogradset;
    bool wavgradset;
    bool sourceset;
    bool dataUxset, dataUyset, dataUzset;
    bool datamodUxset, datamodUyset, datamodUzset;
    bool dataresUxset, dataresUyset, dataresUzset;
    bool dataweightxset;
    bool dataweightyset;
    bool dataweightzset;
};

}
#endif //LSRTM_H
