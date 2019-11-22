#ifndef FAT_H
#define FAT_H

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
#include "rays.h"
#include "der.h"
#include "snap.h"
#include "image.h"
#include "revolve.h"

#define FAT_OK 1
#define FAT_ERR 0

#define GMAP 1
#define SMAP 0

#define PCLIP 97

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define CLOSETO(x, y) ((ABS((x) - (y)) <= FLT_EPSILON*ABS(y)) ? 1 : 0)
#define CUB(x) ((x)*(x)*(x) + 1e-2)

namespace rockseis {


// =============== ABSTRACT FAT CLASS =============== //
/** The abstract fwi class
 *
 */

template<typename T>
class Fat {
public:
    Fat();	///< Constructor
    ~Fat();	///< Destructor
    
    // Fat functions
    std::string getLogfile() { return logfile; } ///< Get name of logfile
    std::string getSnapfile() { return snapfile; } ///< Sets checkpoint filename
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setSnapfile(std::string file) { snapfile = file; } ///< Sets checkpoint filename
    void setMisfit(T val) { misfit = val; } ///< Sets data misfit value
    T getMisfit() { return misfit; }   ///< Gets misfit value
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

private:
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    std::string snapfile;
    T misfit; ///< Misfit value
};

/** The 2D Acoustic Fat class
 *
 */
template<typename T>
class FatAcoustic2D: public Fat<T> {
public:
    FatAcoustic2D();					///< Constructor
    FatAcoustic2D(std::shared_ptr<ModelEikonal2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> Tdata);					///< Constructor 
    int solve(); ///< Runs forward eikonal solver
    int solve_adj(); ///< Runs adjoint eikonal solver
    void setModel(std::shared_ptr<ModelEikonal2D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data2D<T>> _source) { source = _source; sourceset = true; }
    void setVpgrad(std::shared_ptr<Image2D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setTdata(std::shared_ptr<Data2D<T>> _Tdata) { Tdata = _Tdata; Tdataset = true; }
    void setTmod(std::shared_ptr<Data2D<T>> _Tmod) { Tmod = _Tmod; Tmodset = true; }
    void setTres(std::shared_ptr<Data2D<T>> _Tres) { Tres = _Tres; Tresset = true; }
    void setTweight(std::shared_ptr<Data2D<T>> _Tweight) { Tweight = _Tweight; Tweightset = true; }
    void scaleGrad(std::shared_ptr<rockseis::ModelEikonal2D<T>> model, T *lam, T *grad);
    void computeMisfit();
    int run();

    ~FatAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelEikonal2D<T>> model;
    std::shared_ptr<Image2D<T>> vpgrad;
    std::shared_ptr<Data2D<T>> source;
    std::shared_ptr<Data2D<T>> Tdata;
    std::shared_ptr<Data2D<T>> Tmod;
    std::shared_ptr<Data2D<T>> Tres;
    std::shared_ptr<Data2D<T>> Tweight;
    bool modelset;
    bool vpgradset;
    bool sourceset;
    bool Tdataset;
    bool Tmodset;
    bool Tresset;
    bool Tweightset;
};

/** The 3D Acoustic Fat class
 *
 */
template<typename T>
class FatAcoustic3D: public Fat<T> {
public:
    FatAcoustic3D();					///< Constructor
    FatAcoustic3D(std::shared_ptr<ModelEikonal3D<T>> model, std::shared_ptr<Data3D<T>> source, std::shared_ptr<Data3D<T>> Tdata);					///< Constructor 
    int solve(); ///< Runs forward eikonal solver
    int solve_adj(); ///< Runs adjoint eikonal solver
    void setModel(std::shared_ptr<ModelEikonal3D<T>> _model) { model = _model; modelset = true; }
    void setSource(std::shared_ptr<Data3D<T>> _source) { source = _source; sourceset = true; }
    void setVpgrad(std::shared_ptr<Image3D<T>> _vpgrad) { vpgrad = _vpgrad; vpgradset = true; }
    void setTdata(std::shared_ptr<Data3D<T>> _Tdata) { Tdata = _Tdata; Tdataset = true; }
    void setTmod(std::shared_ptr<Data3D<T>> _Tmod) { Tmod = _Tmod; Tmodset = true; }
    void setTres(std::shared_ptr<Data3D<T>> _Tres) { Tres = _Tres; Tresset = true; }
    void setTweight(std::shared_ptr<Data3D<T>> _Tweight) { Tweight = _Tweight; Tweightset = true; }
    void scaleGrad(std::shared_ptr<rockseis::ModelEikonal3D<T>> model, T *lam, T *grad);
    void computeMisfit();
    int run();

    ~FatAcoustic3D();	///< Destructor

private:
    std::shared_ptr<ModelEikonal3D<T>> model;
    std::shared_ptr<Image3D<T>> vpgrad;
    std::shared_ptr<Data3D<T>> source;
    std::shared_ptr<Data3D<T>> Tdata;
    std::shared_ptr<Data3D<T>> Tmod;
    std::shared_ptr<Data3D<T>> Tres;
    std::shared_ptr<Data3D<T>> Tweight;
    bool modelset;
    bool vpgradset;
    bool sourceset;
    bool Tdataset;
    bool Tmodset;
    bool Tresset;
    bool Tweightset;
};


}
#endif //FAT_H
