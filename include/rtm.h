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
    void writeLog(char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

    ~Rtm();	///< Destructor


private:
    int order;  ///< Order of the FD stencil 
    int snapinc;  ///< Snap interval
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    rs_snapmethod snapmethod; 
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
    RtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<ImageAcoustic2D<T>> pimage, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_edge(); ///< Runs rtm with edge boundary saving
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data2D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }


    ~RtmAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<ImageAcoustic2D<T>> pimage;
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
    RtmAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<ImageAcoustic3D<T>> pimage, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm with full snapshoting
    int run_edge(); ///< Runs rtm with edge boundary saving
    int run_optimal(); ///< Runs rtm with optimal checkpointing
    void setModel(std::shared_ptr<ModelAcoustic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setDataP(std::shared_ptr<Data3D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data3D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }

    ~RtmAcoustic3D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic3D<T>> model;
    std::shared_ptr<ImageAcoustic3D<T>> pimage;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> dataP;
    std::shared_ptr<Data3D<T>> dataAz;
    bool modelset;
    bool pimageset;
    bool sourceset;
    bool dataPset, dataAxset, dataAyset, dataAzset;
};


}
#endif //RTM_H
