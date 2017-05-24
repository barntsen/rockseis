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

#define RTM_OK 1
#define RTM_ERR 0

#define GMAP 1
#define SMAP 0

namespace rockseis {

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
    std::string getCpfile() { return cpfile; } ///< Sets checkpoint filename
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil
    void setSnapinc(int _snapinc) {snapinc = _snapinc;} ///< Set snap increment for recording snapshots
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setOptimal(bool val) { optimal = val; } ///< Sets optimal checkpointing flag
    void setIncore(bool val) { incore = val; } ///< Sets optimal checkpoint incore flag
    void setCpfile(std::string file) { cpfile = file; } ///< Sets checkpoint filename
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
    bool optimal; 
    bool incore;
    bool edges;
    std::string cpfile;
};

/** The 2D Acoustic Rtm class
 *
 */
template<typename T>
class RtmAcoustic2D: public Rtm<T> {
public:
    RtmAcoustic2D();					///< Constructor
    RtmAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> dataP, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs rtm
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setDataP(std::shared_ptr<Data2D<T>> _dataP) { dataP = _dataP; dataPset = true; }
    void setDataAz(std::shared_ptr<Data2D<T>> _dataAz) { dataAz = _dataAz; dataAzset = true; }


    ~RtmAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> dataP;
    std::shared_ptr<Data2D<T>> dataAz;
    bool modelset;
    bool sourceset;
    bool dataPset, dataAxset, dataAzset;
};

}
#endif //RTM_H
