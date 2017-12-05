#ifndef MODELLING_H
#define MODELLING_H

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

#define MOD_OK 1
#define MOD_ERR 0

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


// =============== ABSTRACT MODELLING CLASS =============== //
/** The abstract modelling class
 *
 */

template<typename T>
class Modelling {
public:
    Modelling();	///< Constructor
    Modelling(int order, int snapinc);  ///< Constructor
    
    // Modelling functions
    int getOrder() { return order; } ///< Get order of FD stencil
    int getSnapinc() { return snapinc; } ///< Get snap increment
    std::string getLogfile() { return logfile; } ///< Get name of logfile
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil
    void setSnapinc(int _snapinc) {snapinc = _snapinc;} ///< Set snap increment for recording snapshots
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

    ~Modelling();	///< Destructor


private:
    int order;  ///< Order of the FD stencil 
    int snapinc;  ///< Snap interval
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
};

/** The 2D Acoustic Modelling class
 *
 */
template<typename T>
class ModellingAcoustic2D: public Modelling<T> {
public:
    ModellingAcoustic2D();					///< Constructor
    ModellingAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Data2D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelAcoustic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data2D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecAx(std::shared_ptr<Data2D<T>> _recAx) { recAx = _recAx; recAxset = true; }
    void setRecAz(std::shared_ptr<Data2D<T>> _recAz) { recAz = _recAz; recAzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapAx(std::string _snapAx) { snapAx = _snapAx; snapAxset = true; }
    void setSnapAz(std::string _snapAz) { snapAz = _snapAz; snapAzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic2D<T>> model;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> recP;
    std::shared_ptr<Data2D<T>> recAx;
    std::shared_ptr<Data2D<T>> recAz;
    bool modelset;
    bool sourceset;
    bool recPset, recAxset, recAzset;
    std::string snapP, snapAx, snapAz;
    bool snapPset, snapAxset, snapAzset;
};

/** The 3D Acoustic Modelling class
 *
 */
template<typename T>
class ModellingAcoustic3D: public Modelling<T> {
public:
    ModellingAcoustic3D();					///< Constructor
    ModellingAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Data3D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelAcoustic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data3D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecAx(std::shared_ptr<Data3D<T>> _recAx) { recAx = _recAx; recAxset = true; }
    void setRecAy(std::shared_ptr<Data3D<T>> _recAy) { recAy = _recAy; recAyset = true; }
    void setRecAz(std::shared_ptr<Data3D<T>> _recAz) { recAz = _recAz; recAzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapAx(std::string _snapAx) { snapAx = _snapAx; snapAxset = true; }
    void setSnapAy(std::string _snapAy) { snapAy = _snapAy; snapAyset = true; }
    void setSnapAz(std::string _snapAz) { snapAz = _snapAz; snapAzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingAcoustic3D();	///< Destructor

private:
    std::shared_ptr<ModelAcoustic3D<T>> model;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> recP;
    std::shared_ptr<Data3D<T>> recAx;
    std::shared_ptr<Data3D<T>> recAy;
    std::shared_ptr<Data3D<T>> recAz;
    bool modelset;
    bool sourceset;
    bool recPset, recAxset, recAyset, recAzset;
    std::string snapP, snapAx, snapAy, snapAz;
    bool snapPset, snapAxset, snapAyset, snapAzset;
};


/** The 2D Elastic Modelling class
 *
 */
template<typename T>
class ModellingElastic2D: public Modelling<T> {
public:
    ModellingElastic2D();					///< Constructor
    ModellingElastic2D(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Data2D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data2D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecVx(std::shared_ptr<Data2D<T>> _recVx) { recVx = _recVx; recVxset = true; }
    void setRecVz(std::shared_ptr<Data2D<T>> _recVz) { recVz = _recVz; recVzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapSxx(std::string _snapSxx) { snapSxx = _snapSxx; snapSxxset = true; }
    void setSnapSzz(std::string _snapSzz) { snapSzz = _snapSzz; snapSzzset = true; }
    void setSnapSxz(std::string _snapSxz) { snapSxz = _snapSxz; snapSxzset = true; }
    void setSnapVx(std::string _snapVx) { snapVx = _snapVx; snapVxset = true; }
    void setSnapVz(std::string _snapVz) { snapVz = _snapVz; snapVzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingElastic2D();	///< Destructor

private:
    std::shared_ptr<ModelElastic2D<T>> model;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> recP;
    std::shared_ptr<Data2D<T>> recVx;
    std::shared_ptr<Data2D<T>> recVz;
    bool modelset;
    bool sourceset;
    bool recPset, recVxset, recVzset;
    std::string snapP, snapSxx, snapSzz, snapSxz, snapVx, snapVz;
    bool snapPset, snapSxxset, snapSzzset, snapSxzset, snapVxset, snapVzset;
};



/** The 3D Elastic Modelling class
 *
 */
template<typename T>
class ModellingElastic3D: public Modelling<T> {
public:
    ModellingElastic3D();					///< Constructor
    ModellingElastic3D(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Data3D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelElastic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data3D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecVx(std::shared_ptr<Data3D<T>> _recVx) { recVx = _recVx; recVxset = true; }
    void setRecVy(std::shared_ptr<Data3D<T>> _recVy) { recVy = _recVy; recVyset = true; }
    void setRecVz(std::shared_ptr<Data3D<T>> _recVz) { recVz = _recVz; recVzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapSxx(std::string _snapSxx) { snapSxx = _snapSxx; snapSxxset = true; }
    void setSnapSyy(std::string _snapSyy) { snapSyy = _snapSyy; snapSyyset = true; }
    void setSnapSzz(std::string _snapSzz) { snapSzz = _snapSzz; snapSzzset = true; }
    void setSnapSyz(std::string _snapSyz) { snapSyz = _snapSyz; snapSyzset = true; }
    void setSnapSxz(std::string _snapSxz) { snapSxz = _snapSxz; snapSxzset = true; }
    void setSnapSxy(std::string _snapSxy) { snapSxy = _snapSxy; snapSxyset = true; }
    void setSnapVx(std::string _snapVx) { snapVx = _snapVx; snapVxset = true; }
    void setSnapVy(std::string _snapVy) { snapVy = _snapVy; snapVyset = true; }
    void setSnapVz(std::string _snapVz) { snapVz = _snapVz; snapVzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingElastic3D();	///< Destructor

private:
    std::shared_ptr<ModelElastic3D<T>> model;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> recP;
    std::shared_ptr<Data3D<T>> recVx;
    std::shared_ptr<Data3D<T>> recVy;
    std::shared_ptr<Data3D<T>> recVz;
    bool modelset;
    bool sourceset;
    bool recPset, recVxset, recVyset, recVzset;
    std::string snapP, snapSxx, snapSyy, snapSzz, snapSyz, snapSxz, snapSxy, snapVx, snapVy, snapVz;
    bool snapPset, snapSxxset, snapSyyset, snapSzzset, snapSxzset, snapSyzset, snapSxyset, snapVxset, snapVyset, snapVzset;
};

/** The 2D Elastic Displancement-Stress Modelling class
 *
 */
template<typename T>
class ModellingElastic2D_DS: public Modelling<T> {
public:
    ModellingElastic2D_DS();					///< Constructor
    ModellingElastic2D_DS(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Data2D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelElastic2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data2D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecUx(std::shared_ptr<Data2D<T>> _recUx) { recUx = _recUx; recUxset = true; }
    void setRecUz(std::shared_ptr<Data2D<T>> _recUz) { recUz = _recUz; recUzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapSxx(std::string _snapSxx) { snapSxx = _snapSxx; snapSxxset = true; }
    void setSnapSzz(std::string _snapSzz) { snapSzz = _snapSzz; snapSzzset = true; }
    void setSnapSxz(std::string _snapSxz) { snapSxz = _snapSxz; snapSxzset = true; }
    void setSnapUx(std::string _snapUx) { snapUx = _snapUx; snapUxset = true; }
    void setSnapUz(std::string _snapUz) { snapUz = _snapUz; snapUzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingElastic2D_DS();	///< Destructor

private:
    std::shared_ptr<ModelElastic2D<T>> model;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> recP;
    std::shared_ptr<Data2D<T>> recUx;
    std::shared_ptr<Data2D<T>> recUz;
    bool modelset;
    bool sourceset;
    bool recPset, recUxset, recUzset;
    std::string snapP, snapSxx, snapSzz, snapSxz, snapUx, snapUz;
    bool snapPset, snapSxxset, snapSzzset, snapSxzset, snapUxset, snapUzset;
};

/** The 3D Displacement Stress Elastic Modelling class
 *
 */
template<typename T>
class ModellingElastic3D_DS: public Modelling<T> {
public:
    ModellingElastic3D_DS();					///< Constructor
    ModellingElastic3D_DS(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Data3D<T>> source, int order, int snapinc);					///< Constructor 
    int run(); ///< Runs modelling
    void setModel(std::shared_ptr<ModelElastic3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setRecP(std::shared_ptr<Data3D<T>> _recP) { recP = _recP; recPset = true; }
    void setRecUx(std::shared_ptr<Data3D<T>> _recUx) { recUx = _recUx; recUxset = true; }
    void setRecUy(std::shared_ptr<Data3D<T>> _recUy) { recUy = _recUy; recUyset = true; }
    void setRecUz(std::shared_ptr<Data3D<T>> _recUz) { recUz = _recUz; recUzset = true; }

    void setSnapP(std::string _snapP) { snapP = _snapP; snapPset = true; }
    void setSnapSxx(std::string _snapSxx) { snapSxx = _snapSxx; snapSxxset = true; }
    void setSnapSyy(std::string _snapSyy) { snapSyy = _snapSyy; snapSyyset = true; }
    void setSnapSzz(std::string _snapSzz) { snapSzz = _snapSzz; snapSzzset = true; }
    void setSnapSyz(std::string _snapSyz) { snapSyz = _snapSyz; snapSyzset = true; }
    void setSnapSxz(std::string _snapSxz) { snapSxz = _snapSxz; snapSxzset = true; }
    void setSnapSxy(std::string _snapSxy) { snapSxy = _snapSxy; snapSxyset = true; }
    void setSnapUx(std::string _snapUx) { snapUx = _snapUx; snapUxset = true; }
    void setSnapUy(std::string _snapUy) { snapUy = _snapUy; snapUyset = true; }
    void setSnapUz(std::string _snapUz) { snapUz = _snapUz; snapUzset = true; }
    bool checkStability(); ///< Check stability of finite difference modelling

    ~ModellingElastic3D_DS();	///< Destructor

private:
    std::shared_ptr<ModelElastic3D<T>> model;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> recP;
    std::shared_ptr<Data3D<T>> recUx;
    std::shared_ptr<Data3D<T>> recUy;
    std::shared_ptr<Data3D<T>> recUz;
    bool modelset;
    bool sourceset;
    bool recPset, recUxset, recUyset, recUzset;
    std::string snapP, snapSxx, snapSyy, snapSzz, snapSyz, snapSxz, snapSxy, snapUx, snapUy, snapUz;
    bool snapPset, snapSxxset, snapSyyset, snapSzzset, snapSxzset, snapSyzset, snapSxyset, snapUxset, snapUyset, snapUzset;
};

}
#endif //MODELLING_H
