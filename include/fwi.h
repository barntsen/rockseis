#ifndef FWI_H
#define FWI_H

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

#define FWI_OK 1
#define FWI_ERR 0

#define GMAP 1
#define SMAP 0

#define ki2D(i,j) ((j)*nx +(i))
#define km2D(i,j) ((j)*nx + (i))
#define ks2D(i,j) ((j)*nxs + (i))
#define kr2D(i,j) ((j)*nxr + (i))

#define ki3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define km3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define ks3D(i,j,k) ((k)*nxs*nys + (j)*nxs + (i))
#define kr3D(i,j,k) ((k)*nxr*nyr + (j)*nxr + (i))

#define kwav(i,j) ((j)*nt + (i))

namespace rockseis {

// =============== ENUMS =============== //
typedef enum {FULL, OPTIMAL, EDGES} rs_snapmethod; ///< Snapshot saving method

/** The Progress struct
 *
 */
typedef struct{
	clock_t previous, current; ///< Time book keeping
	float persec; ///< Iterations per second 
	char speed[48]; ///< Iterations per second string
	char progress[128]; ///< Progress string
} Progress;


// =============== ABSTRACT FWI CLASS =============== //
/** The abstract fwi class
 *
 */

template<typename T>
class Fwi {
public:
    Fwi();	///< Constructor
    Fwi(int order, int snapinc);  ///< Constructor
    
    // Fwi functions
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

    ~Fwi();	///< Destructor


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
};

/** The 2D Acoustic Fwi class
 *
 */
template<typename T>
class FwiAcoustic2D: public Fwi<T> {
public:
    FwiAcoustic2D();					///< Constructor
    FwiAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_edge(); ///< Runs fwi with edge boundary saving
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setWavgrad(std::shared_ptr<Data2D<T>> _wavgrad) { wavgrad = _wavgrad; wavgradset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setRhograd(std::shared_ptr<Image2D<T>> _rhograd) { rhograd = _rhograd; rhogradset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data2D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }
    void setDatamodP(std::shared_ptr<Data2D<T>> _datamodP) { datamodP = _datamodP; datamodPset = true; }
    void setDatamodAz(std::shared_ptr<Data2D<T>> _datamodAz) { datamodAz = _datamodAz; datamodAzset = true; }
    void setDataresP(std::shared_ptr<Data2D<T>> _dataresP) { dataresP = _dataresP; dataresPset = true; }
    void setDataresAz(std::shared_ptr<Data2D<T>> _dataresAz) { dataresAz = _dataresAz; dataresAzset = true; }
    void setDataweight(std::shared_ptr<Data2D<T>> _dataweight) { dataweight = _dataweight; dataweightset = true; }

    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wrz, int padr, T *vp, T* rho);

    void computeMisfit();
    void computeResiduals();
    ~FwiAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Image2D<T>> rhograd;
    std::shared_ptr<Data2D<T>> wavgrad;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataP;
    std::shared_ptr<Data2D<T>> dataAz;
    std::shared_ptr<Data2D<T>> datamodP;
    std::shared_ptr<Data2D<T>> datamodAz;
    std::shared_ptr<Data2D<T>> dataresP;
    std::shared_ptr<Data2D<T>> dataresAz;
    std::shared_ptr<Data2D<T>> dataweight;
    bool modelset;
    bool wavgradset;
    bool vpgradset;
    bool rhogradset;
    bool sourceset;
    bool dataPset, dataAzset;
    bool datamodPset, datamodAzset;
    bool dataresPset, dataresAzset;
    bool dataweightset;
};

/** The 3D Acoustic Fwi class
 *
 */
template<typename T>
class FwiAcoustic3D: public Fwi<T> {
public:
    FwiAcoustic3D();					///< Constructor
    FwiAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataP, int order, int snapinc);					///< Constructor 
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
    void crossCorr(T* wsp, int pads, T* wrp, T* wrx, T* wry, T* wrz, int padr, T *vp, T* rho);
    void computeResiduals();

    ~FwiAcoustic3D();	///< Destructor

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

/** The 2D Elastic Fwi class
 *
 */
template<typename T>
class FwiElastic2D: public Fwi<T> {
public:
    FwiElastic2D();					///< Constructor
    FwiElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataUx, std::shared_ptr<Data2D<T>> dataUz, int order, int snapinc);					///< Constructor 
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
    void setDataweight(std::shared_ptr<Data2D<T>> _dataweight) { dataweight = _dataweight; dataweightset = true; }
    void crossCorr(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D_DS<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it);
    void crossCorr2(T *wsx, T *wsz, int pads,std::shared_ptr<WavesElastic2D<T>> waves_bw, std::shared_ptr<ModelElastic2D<T>> model, int it);
    void scaleGrad(std::shared_ptr<ModelElastic2D<T>> model);
    void computeMisfit();
    void computeResiduals();
    void computeResiduals2();

    ~FwiElastic2D();	///< Destructor

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
    std::shared_ptr<Data2D<T>> dataweight;
    bool dataweightset;
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

/** The 3D Elastic Fwi class
 *
 */
template<typename T>
class FwiElastic3D: public Fwi<T> {
public:
    FwiElastic3D();					///< Constructor
    FwiElastic3D(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataVx, std::shared_ptr<Data3D<T>> dataVy, std::shared_ptr<Data3D<T>> dataVz, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs fwi with full snapshoting
    int run_optimal(); ///< Runs fwi with optimal checkpointing
    void setModel(std::shared_ptr<ModelElastic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setDataVx(std::shared_ptr<Data3D<T>> _dataVx) { dataVx = _dataVx; dataVxset = true; }
    void setDataVy(std::shared_ptr<Data3D<T>> _dataVy) { dataVy = _dataVy; dataVyset = true; }
    void setDataVz(std::shared_ptr<Data3D<T>> _dataVz) { dataVz = _dataVz; dataVzset = true; }
    void setDatamodVx(std::shared_ptr<Data3D<T>> _datamodVx) { datamodVx = _datamodVx; datamodVxset = true; }
    void setDatamodVy(std::shared_ptr<Data3D<T>> _datamodVy) { datamodVy = _datamodVy; datamodVyset = true; }
    void setDatamodVz(std::shared_ptr<Data3D<T>> _datamodVz) { datamodVz = _datamodVz; datamodVzset = true; }
    void setDataresVx(std::shared_ptr<Data3D<T>> _dataresVx) { dataresVx = _dataresVx; dataresVxset = true; }
    void setDataresVy(std::shared_ptr<Data3D<T>> _dataresVy) { dataresVy = _dataresVy; dataresVyset = true; }
    void setDataresVz(std::shared_ptr<Data3D<T>> _dataresVz) { dataresVz = _dataresVz; dataresVzset = true; }

    void setWavgrad(std::shared_ptr<Data3D<T>> _wavgrad) { wavgrad = _wavgrad; wavgradset = true; }
    void setVpgrad(std::shared_ptr<Image3D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setVsgrad(std::shared_ptr<Image3D<T>> _vsgrad) { vsgrad = _vsgrad; vsgradset = true; }
    void setRhograd(std::shared_ptr<Image3D<T>> _rhograd) { rhograd = _rhograd; rhogradset = true; }
    void crossCorr(T *wsx, T*wsy, T *wsz, int pads, std::shared_ptr<WavesElastic3D<T>> waves_bw, std::shared_ptr<ModelElastic3D<T>> model, int it);
    void computeResiduals();

    ~FwiElastic3D();	///< Destructor

private:
    std::shared_ptr<ModelElastic3D<T>> model;
    std::shared_ptr<Image3D<T>> vpgrad;
    std::shared_ptr<Image3D<T>> vsgrad;
    std::shared_ptr<Image3D<T>> rhograd;
    std::shared_ptr<Data3D<T>> wavgrad;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataVx;
    std::shared_ptr<Data3D<T>> datamodVx;
    std::shared_ptr<Data3D<T>> dataresVx;
    std::shared_ptr<Data3D<T>> dataVy;
    std::shared_ptr<Data3D<T>> datamodVy;
    std::shared_ptr<Data3D<T>> dataresVy;
    std::shared_ptr<Data3D<T>> dataVz;
    std::shared_ptr<Data3D<T>> datamodVz;
    std::shared_ptr<Data3D<T>> dataresVz;
    std::shared_ptr<Data3D<T>> dataweight;
    bool modelset;
    bool vpgradset;
    bool vsgradset;
    bool rhogradset;
    bool wavgradset;
    bool sourceset;
    bool dataVxset, dataVyset, dataVzset;
    bool datamodVxset, datamodVyset, datamodVzset;
    bool dataresVxset, dataresVyset, dataresVzset;
    bool dataweightset;
};


}
#endif //FWI_H
