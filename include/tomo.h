#ifndef TOMO_H
#define TOMO_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "geometry.h"
#include "utils.h"
#include "file.h"
#include "model.h"
#include "data.h"
#include "rays.h"
#include "waves.h"
#include "der.h"
#include "snap.h"
#include "interp.h"
#include "parallel.h"
#include "fat.h"
#include "sort.h"
#include "bspl.h"

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

#define TOM_ERR 0
#define TOM_OK 1

#define LOGFILE "tomo.log"

#define PROGLOGFILE "progress.log"

#define VPLSFILE "vp_ls.rss"

#define VP0FILE "vp_0.rss"

#define VPGRADFILE "vp_grad.rss"

#define VPREGGRADFILE "vp_reg_grad.rss"

#define VPGRADCOMBFILE "vp_grad_comb.rss"

#define VPGRADMUTEFILE "vp_grad_muted.rss"

#define VPPROJGRADFILE "vp_proj_grad.rss"

#define MISFITFILE "misfit.rss"

#define VPREGMISFITFILE "vpreg_misfit.rss"

#define TMODFILE "tmod.rss"
#define TRESFILE "tres.rss"

#define VP_UP "vp_up.rss"

#define RESULTDIR "Results"

#define ABS(x) ((x) < 0 ? -(x) : (x))


namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PAR_GRID, PAR_BSPLINE, PAR_1D} rs_paramtype; ///< Type of parameterisation

// ##### TOMO CLASS
template<typename T>
class Tomo {
public:
    Tomo(); ///<Constructor
    Tomo(MPImodeling *_mpi); ///<Constructor
    ~Tomo(); ///<Destructor
    void createResult();
    double vector_norm(double *v, const int type, const int n);
    void setMpi(MPImodeling *_mpi) { mpi = _mpi; }
    void writeLog(std::string msg);
    void writeProgress(std::string msg);
    MPImodeling * getMpi() { return mpi; }
    void setLpml(int val) { lpml = val; }
    int getLpml() { return lpml; }
    bool getFilter() { return filter; }
    void setIncore(bool val) { incore = val; } 
    bool getIncore() { return incore; }
    int getOrder() { return order; }
    void setOrder(int val) { order = val; }
    int getSnapinc() { return snapinc; }
    void setSnapinc(int val) { snapinc = val; }
    int getNsnaps() { return nsnaps; }
    void setNsnaps(int val) { nsnaps = val; }
    rs_paramtype getParamtype() { return paramtype; } 
    void setParamtype(rs_paramtype val) { paramtype = val; }
    void setMisfit(T val) { misfit = val; } ///< Sets data misfit value
    T getMisfit() { return misfit; }   ///< Gets misfit value

    T getTmax() { return tmax; }   ///< Gets tmax
    void setTmax(T val) { tmax = val; }
    T getDtx() { return dtx; }
    T getDty() { return dty; }
    T getDtz() { return dtz; }
    void setDtx(T val) { dtx = val; }
    void setDty(T val) { dty = val; }
    void setDtz(T val) { dtz = val; }
    void setFnorm(double norm) { fnorm = norm; }
    double getFnorm() { return fnorm; }
    void normalize(double *v, double *f, int n);
    void setNoreverse(bool val) { noreverse = val; }
    bool getNoreverse() { return noreverse; }


private:
	int lpml;
    bool incore;
	int order;
	int snapinc;
	int nsnaps;
    bool filter;
    double fnorm;
    std::string logfile; ///< Log file
    std::string progressfile; ///< Log file
    std::string progresslogfile; ///< Log file
    std::ofstream Flog; ///< Log stream
    bool createLog(); ///< Set name of logfile and open for writing
    bool createProglog(); ///< Set name of progress logfile and open for writing
    rs_paramtype paramtype;
    MPImodeling *mpi;
    T dtx, dty, dtz;
    bool noreverse;
    T misfit; ///< Misfit value
    T tmax; ///< Maximum traveltime in domain
};

// ##### ACOUSTIC 2D TOMO CLASS
template<typename T>
class TomoAcoustic2D: public Tomo<T> {
public:
    TomoAcoustic2D(); ///<Constructor
    TomoAcoustic2D(MPImodeling *_mpi); ///<Constructor
    ~TomoAcoustic2D(); ///<Destructor

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setTrecordfile(std::string file) { Trecordfile = file; }
    std::string getTrecordfile() { return Trecordfile; }

    void setTmodelledfile(std::string file) { Tmodelledfile = file; }
    std::string getTmodelledfile() { return Tmodelledfile; }

    void setTresidualfile(std::string file) { Tresidualfile = file; }
    std::string getTresidualfile() { return Tresidualfile; }

    void setMutefile(std::string file) { Mutefile = file; }
    std::string getMutefile() { return Mutefile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setVpregalpha(T val) { reg_alpha[0] = val; }
    T getVpregalpha() { return reg_alpha[0]; }

    void setVpregeps(T val) { reg_eps[0] = val; }
    T getVpregeps() { return reg_eps[0]; }

    void setDataweight(bool val) {dataweight = val; }

    void setDataweightfile(std::string file) {Dataweightfile = file; }

    // Run gradient
    void runGrad();

    void clipGrad(std::shared_ptr<rockseis::Image2D<T>> grad);
   
    // Run BSProjection
    void runBsproj();

    // Mute gradient
    void applyMute();

    // Regularisation computation
    void computeRegularisation(double *x);

    // Combine gradients
    void combineGradients();

    // Set initial
    int setInitial(double *x, std::string vpfile);

    // Save line search models
    void saveLinesearch(double *x);

    // Save Results
    void saveResults(int i);

    // Read gradient
    void readGrad(double *g);

    // Read misfit
    void readMisfit(double *f);


private:
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vpgradfile;
    std::string Misfitfile;
    std::string Snapfile;
    std::string Trecordfile;
    std::string Tmodelledfile;
    std::string Tresidualfile;
    std::string Dataweightfile;
    bool dataweight;
    std::string Mutefile;
    T apertx;
    T kvp;
    T reg_eps[2];
    T reg_alpha[2];
};

}
#endif //TOMO_H
