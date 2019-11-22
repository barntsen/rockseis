#ifndef KDMIG_H
#define KDMIG_H

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

#define KDMIG_OK 1
#define KDMIG_ERR 0

#define GMAP 1
#define SMAP 0

#define ki2D(i,j,k,l) ((l)*nhx*nz*nx + (k)*nx*nz + (j)*nx +(i))
#define km2D(i,j) ((j)*nx + (i))
#define kt2D(i,j) ((j)*nxt + (i))

#define ki3D(i,j,k,l,m,n) ((n)*nhy*nhx*nx*ny*nz + (m)*nhx*nx*ny*nz + (l)*nx*ny*nz + (k)*nx*ny + (j)*nx + (i))
#define km3D(i,j,k) ((k)*nx*ny + (j)*nx + (i))
#define kt3D(i,j,k) ((k)*nxt*nyt + (j)*nxt + (i))

#define PCLIP 97

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define CLOSETO(x, y) ((ABS((x) - (y)) <= FLT_EPSILON*ABS(y)) ? 1 : 0)
#define CUB(x) ((x)*(x)*(x) + 1e-2)
#define CMULR(ra,ia,rb,ib) (((ra)*(rb)) - ((ia)*(ib)))
#define CMULI(ra,ia,rb,ib) (((ia)*(rb)) + ((ra)*(ib)))

namespace rockseis {


// =============== ABSTRACT KDMIG CLASS =============== //
/** The abstract kdmig class
 *
 */

template<typename T>
class Kdmig {
public:
    Kdmig();	///< Constructor
    ~Kdmig();	///< Destructor
    
    // Kdmig functions
    std::string getLogfile() { return logfile; } ///< Get name of logfile
    std::string getSnapfile() { return snapfile; } ///< Sets checkpoint filename
    void setLogfile(std::string name) { logfile = name; } ///< Set name of logfile
    void setSnapfile(std::string file) { snapfile = file; } ///< Sets checkpoint filename
    bool createLog(std::string name); ///< Set name of logfile and open for writing
    void writeLog(std::string text);  ///< Write string to log file
    void writeLog(const char * text); ///< Write c_string to log file
    int getFreqinc() { return freqinc; } ///< Get freq increment
    void setFreqinc(int _freqinc) {freqinc = _freqinc;} ///< Set freq increment for recording freqshots
    void writeProgressbar(int x, int n, int r, int w);
    void writeProgress(int x, int n, int r, int w);

private:
    std::string logfile; ///< Log file name
    std::ofstream Flog; ///< Logfile
	Progress prog; ///< Progress counter
    std::string snapfile;
    int freqinc;  ///< Frequency decimation interval
};

/** The 2D Acoustic Kdmig class
 *
 */
template<typename T>
class KdmigAcoustic2D: public Kdmig<T> {
public:
    KdmigAcoustic2D();					///< Constructor
    KdmigAcoustic2D(std::shared_ptr<ModelEikonal2D<T>> model, std::shared_ptr<Data2D<T>> data, std::shared_ptr<Image2D<T>> pimage);					///< Constructor 
    int solve(); ///< Runs forward eikonal solver
    int solve_adj(); ///< Runs adjoint eikonal solver
    void setModel(std::shared_ptr<ModelEikonal2D<T>> _model) { model = _model; modelset = true; }
    void setData(std::shared_ptr<Data2D<T>> _data) { data = _data; dataset = true; }
    void crossCorr(std::shared_ptr<RaysAcoustic2D<T>> rays_sou, std::shared_ptr<RaysAcoustic2D<T>> rays_rec, T* cdata, unsigned long nfs, T df, T ot, int pad);
    int run();

    ~KdmigAcoustic2D();	///< Destructor

private:
    std::shared_ptr<ModelEikonal2D<T>> model;
    std::shared_ptr<Data2D<T>> data;
    std::shared_ptr<Image2D<T>> pimage;
    bool modelset;
    bool dataset;
    bool pimageset;
};

/** The 2D Elastic Kdmig class
 *
 */
template<typename T>
class KdmigElastic2D: public Kdmig<T> {
public:
    KdmigElastic2D();					///< Constructor
    KdmigElastic2D(std::shared_ptr<ModelEikonal2D<T>> vpmodel, std::shared_ptr<ModelEikonal2D<T>> vsmodel, std::shared_ptr<Data2D<T>> data, std::shared_ptr<Image2D<T>> pimage);					///< Constructor 
    int solve(); ///< Runs forward eikonal solver
    int solve_adj(); ///< Runs adjoint eikonal solver
    void setVpmodel(std::shared_ptr<ModelEikonal2D<T>> _vpmodel) { vpmodel = _vpmodel; vpmodelset = true; }
    void setVsmodel(std::shared_ptr<ModelEikonal2D<T>> _vsmodel) { vsmodel = _vsmodel; vsmodelset = true; }
    void setData(std::shared_ptr<Data2D<T>> _data) { data = _data; dataset = true; }
    void crossCorr(std::shared_ptr<RaysAcoustic2D<T>> rays_sou, std::shared_ptr<RaysAcoustic2D<T>> rays_rec, T* cdata, unsigned long nfs, T df, T ot, int pad);
    int run();

    ~KdmigElastic2D();	///< Destructor

private:
    std::shared_ptr<ModelEikonal2D<T>> vpmodel;
    std::shared_ptr<ModelEikonal2D<T>> vsmodel;
    std::shared_ptr<Data2D<T>> data;
    std::shared_ptr<Image2D<T>> simage;
    bool vpmodelset;
    bool vsmodelset;
    bool dataset;
    bool simageset;
};

}
#endif //KDMIG_H