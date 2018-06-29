#ifndef WEMVA_H
#define WEMVA_H

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
#include "waves.h"
#include "der.h"
#include "snap.h"
#include "interp.h"
#include "parallel.h"
#include "rtm.h"
#include "sort.h"
#include "bspl.h"

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

#define WVA_ERR 0
#define WVA_OK 1

#define LOGFILE "wemva.log"

#define PROGLOGFILE "progress.log"

#define VPLSFILE "vp_ls.rss"
#define RHOLSFILE "rho_ls.rss"
#define VSLSFILE "vs_ls.rss"
#define SOURCELSFILE "source_ls.rss"

#define VP0FILE "vp_0.rss"
#define RHO0FILE "rho_0.rss"
#define VS0FILE "vs_0.rss"
#define SOURCE0FILE "source_0.rss"

#define PIMAGEFILE "pimage.rss"
#define SIMAGEFILE "simage.rss"

#define VPGRADFILE "vp_grad.rss"
#define VSGRADFILE "vs_grad.rss"

#define VPREGGRADFILE "vp_reg_grad.rss"
#define VSREGGRADFILE "vs_reg_grad.rss"

#define SRCILUMFILE "src_ilum.rss"

#define VPGRADCOMBFILE "vp_grad_comb.rss"
#define VSGRADCOMBFILE "vs_grad_comb.rss"

#define VPGRADMUTEFILE "vp_grad_muted.rss"
#define VSGRADMUTEFILE "vs_grad_muted.rss"

#define VPPROJGRADFILE "vp_proj_grad.rss"
#define VSPROJGRADFILE "vs_proj_grad.rss"

#define MISFITFILE "misfit.rss"
#define VPREGMISFITFILE "vpreg_misfit.rss"
#define VSREGMISFITFILE "vsreg_misfit.rss"

#define VP_UP "vp_up.rss"
#define VS_UP "vs_up.rss"

#define RESULTDIR "Results"



namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PAR_GRID, PAR_BSPLINE, PAR_1D} rs_paramtype; ///< Type of parameterisation

// ##### WEMVA CLASS
template<typename T>
class Wemva {
public:
    Wemva(); ///<Constructor
    Wemva(MPImodeling *_mpi); ///<Constructor
    ~Wemva(); ///<Destructor
    void createResult();
    double vector_norm(double *v, const int type, const int n);
    void setMpi(MPImodeling *_mpi) { mpi = _mpi; }
    void writeLog(std::string msg);
    void writeProgress(std::string msg);
    MPImodeling * getMpi() { return mpi; }
    void setLpml(int val) { lpml = val; }
    int getLpml() { return lpml; }
    void setFs(bool val) { fs = val; } 
    bool getFs() { return fs; }
    bool getFilter() { return filter; }
    void setIncore(bool val) { incore = val; } 
    bool getIncore() { return incore; }
    int getOrder() { return order; }
    void setOrder(int val) { order = val; }
    int getSnapinc() { return snapinc; }
    void setSnapinc(int val) { snapinc = val; }
    int getNsnaps() { return nsnaps; }
    void setNsnaps(int val) { nsnaps = val; }
    rs_snapmethod getSnapmethod() { return snapmethod; } 
    void setSnapmethod(rs_snapmethod val) { snapmethod = val; }
    rs_paramtype getParamtype() { return paramtype; } 
    void setParamtype(rs_paramtype val) { paramtype = val; }
    T getDtx() { return dtx; }
    T getDty() { return dty; }
    T getDtz() { return dtz; }
    void setDtx(T val) { dtx = val; }
    void setDty(T val) { dty = val; }
    void setDtz(T val) { dtz = val; }
    int getNhx() { return nhx; }
    int getNhy() { return nhy; }
    int getNhz() { return nhz; }
    void setNhx(int val) { nhx = val; }
    void setNhy(int val) { nhy = val; }
    void setNhz(int val) { nhz = val; }
    void setFnorm(double norm) { fnorm = norm; }
    double getFnorm() { return fnorm; }
    void normalize(double *v, double *f, int n);
    void setNoreverse(bool val) { noreverse = val; }
    bool getNoreverse() { return noreverse; }


private:
	int lpml;
	bool fs;
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
    rs_snapmethod snapmethod;
    rs_paramtype paramtype;
    MPImodeling *mpi;
    T dtx, dty, dtz;
    int nhx,nhy,nhz;
    bool noreverse;
};

// ##### ACOUSTIC 2D WEMVA CLASS
template<typename T>
class WemvaAcoustic2D: public Wemva<T> {
public:
    WemvaAcoustic2D(); ///<Constructor
    WemvaAcoustic2D(MPImodeling *_mpi); ///<Constructor
    ~WemvaAcoustic2D(); ///<Destructor

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setRhofile(std::string file) { Rhofile = file; }
    std::string getRhofile() { return Rhofile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setPimagefile(std::string file) { Pimagefile = file; }
    std::string getPimagefile() { return Pimagefile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setPsnapfile(std::string file) { Psnapfile = file; }
    std::string getPsnapfile() { return Psnapfile; }

    void setPrecordfile(std::string file) { Precordfile = file; }
    std::string getPrecordfile() { return Precordfile; }

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

    // Run gradient
    void runGrad();
   
    // Run BSProjection
    void runBsproj();

    // Mute gradient
    void applyMute();

    // Regularisation computation
    void computeRegularisation(double *x);

    // Combine gradients
    void combineGradients();

    // Set initial
    int setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile);

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
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Pimagefile;
    std::string Misfitfile;
    std::string Psnapfile;
    std::string Precordfile;
    std::string Mutefile;
    T apertx;
    T kvp;
    T reg_eps[2];
    T reg_alpha[2];
};

//// ##### ACOUSTIC 3D WEMVA CLASS
//template<typename T>
//class WemvaAcoustic3D: public Wemva<T> {
//public:
//    WemvaAcoustic3D(); ///<Constructor
//    WemvaAcoustic3D(MPImodeling *_mpi); ///<Constructor
//    ~WemvaAcoustic3D(); ///<Destructor
//
//    void setDataweight(bool val) { dataweight = val; }
//    bool getDataweight() { return dataweight; }
//
//    void setDataweightfile(std::string file) { Dataweightfile = file; }
//    std::string getDataweightfile() { return Dataweightfile; }
//
//    void setWaveletfile(std::string file) { Waveletfile = file; }
//    std::string getWaveletfile() { return Waveletfile; }
//
//    void setVpfile(std::string file) { Vpfile = file; }
//    std::string getVpfile() { return Vpfile; }
//
//    void setRhofile(std::string file) { Rhofile = file; }
//    std::string getRhofile() { return Rhofile; }
//
//    void setVpgradfile(std::string file) { Vpgradfile = file; }
//    std::string getVpgradfile() { return Vpgradfile; }
//
//    void setSrcilumfile(std::string file) { Srcilumfile = file; }
//    std::string getSrcilumfile() { return Srcilumfile; }
//
//    void setMisfitfile(std::string file) { Misfitfile = file; }
//    std::string getMisfitfile() { return Misfitfile; }
//
//    void setPsnapfile(std::string file) { Psnapfile = file; }
//    std::string getPsnapfile() { return Psnapfile; }
//
//    void setPrecordfile(std::string file) { Precordfile = file; }
//    std::string getPrecordfile() { return Precordfile; }
//
//    void setPmodelledfile(std::string file) { Pmodelledfile = file; }
//    std::string getPmodelledfile() { return Pmodelledfile; }
//
//    void setMutefile(std::string file) { Mutefile = file; }
//    std::string getMutefile() { return Mutefile; }
//
//    void setPresidualfile(std::string file) { Presidualfile = file; }
//    std::string getPresidualfile() { return Presidualfile; }
//
//    void setApertx(T val) { apertx = val; }
//    T getApertx() { return apertx; }
//
//    void setAperty(T val) { aperty = val; }
//    T getAperty() { return aperty; }
//
//    void setKvp(T val) { kvp = val; }
//    T getKvp() { return kvp; }
//
//    void setKrho(T val) { krho = val; }
//    T getKrho() { return krho; }
//
//    void setKsource(T val) { ksource = val; }
//    T getKsource() { return ksource; }
//
//    void setVpregalpha(T val) { reg_alpha[0] = val; }
//    T getVpregalpha() { return reg_alpha[0]; }
//
//    void setVpregeps(T val) { reg_eps[0] = val; }
//    T getVpregeps() { return reg_eps[0]; }
//
//    // Run gradient
//    void runGrad();
//   
//    // Run BSProjection
//    void runBsproj();
//
//    // Mute gradient
//    void applyMute();
//
//    // Correct for source ilumination
//    void applySrcilum();
//
//    // Regularisation computation
//    void computeRegularisation(double *x);
//
//    // Combine gradients
//    void combineGradients();
//
//    // Set initial
//    int setInitial(double *x, std::string vpfile, std::string rhofile, std::string sourcefile);
//
//    // Save line search models
//    void saveLinesearch(double *x);
//
//    // Save Results
//    void saveResults(int i);
//
//    // Read gradient
//    void readGrad(double *g);
//
//    // Read misfit
//    void readMisfit(double *f);
//
//    void setUpdates(bool vp, bool rho, bool source) { update_vp = vp; update_rho = rho; update_source = source; }
//    void setUpdate_vp(bool vp) { update_vp = vp;}
//    void setUpdate_rho(bool rho) { update_rho = rho;}
//    void setUpdate_source(bool source) { update_source = source;}
//    void setSrcilum(bool val) { srcilumset = val;}
//    bool getSrcilum() { return srcilumset;}
//
//private:
//    bool update_vp;
//    bool update_rho;
//    bool update_source;
//    bool dataweight;
//    bool srcilumset;
//    std::string Waveletfile;
//    std::string Vpfile;
//    std::string Rhofile;
//    std::string Vpgradfile;
//    std::string Dataweightfile;
//    std::string Misfitfile;
//    std::string Psnapfile;
//    std::string Precordfile;
//    std::string Pmodelledfile;
//    std::string Presidualfile;
//    std::string Mutefile;
//    std::string Srcilumfile;
//    T apertx, aperty;
//    T kvp, krho, ksource;
//    T reg_eps[2];
//    T reg_alpha[2];
//};
//
//// ##### ELASTIC 2D WEMVA CLASS
//template<typename T>
//class WemvaElastic2D: public Wemva<T> {
//public:
//    WemvaElastic2D(); ///<Constructor
//    WemvaElastic2D(MPImodeling *_mpi); ///<Constructor
//    ~WemvaElastic2D(); ///<Destructor
//
//    void setDataweightx(bool val) { dataweightx = val; }
//    bool getDataweightx() { return dataweightx; }
//
//    void setDataweightz(bool val) { dataweightz = val; }
//    bool getDataweightz() { return dataweightz; }
//
//    void setDataweightxfile(std::string file) { Dataweightxfile = file; }
//    std::string getDataweightxfile() { return Dataweightxfile; }
//
//    void setDataweightzfile(std::string file) { Dataweightzfile = file; }
//    std::string getDataweightzfile() { return Dataweightzfile; }
//
//    void setWaveletfile(std::string file) { Waveletfile = file; }
//    std::string getWaveletfile() { return Waveletfile; }
//
//    void setVpfile(std::string file) { Vpfile = file; }
//    std::string getVpfile() { return Vpfile; }
//
//    void setSrcilumfile(std::string file) { Srcilumfile = file; }
//    std::string getSrcilumfile() { return Srcilumfile; }
//
//    void setVsfile(std::string file) { Vsfile = file; }
//    std::string getVsfile() { return Vsfile; }
//
//    void setRhofile(std::string file) { Rhofile = file; }
//    std::string getRhofile() { return Rhofile; }
//
//    void setVpgradfile(std::string file) { Vpgradfile = file; }
//    std::string getVpgradfile() { return Vpgradfile; }
//
//    void setVsgradfile(std::string file) { Vsgradfile = file; }
//    std::string getVsgradfile() { return Vsgradfile; }
//
//    void setMisfitfile(std::string file) { Misfitfile = file; }
//    std::string getMisfitfile() { return Misfitfile; }
//
//    void setUxrecordfile(std::string file) { Uxrecordfile = file; }
//    std::string getUxrecordfile() { return Uxrecordfile; }
//
//    void setUxmodelledfile(std::string file) { Uxmodelledfile = file; }
//    std::string getUxmodelledfile() { return Uxmodelledfile; }
//
//    void setUxresidualfile(std::string file) { Uxresidualfile = file; }
//    std::string getUxresidualfile() { return Uxresidualfile; }
//
//    void setSnapfile(std::string file) { Snapfile = file; }
//    std::string getSnapfile() { return Snapfile; }
//
//    void setUzrecordfile(std::string file) { Uzrecordfile = file; }
//    std::string getUzrecordfile() { return Uzrecordfile; }
//
//    void setUzmodelledfile(std::string file) { Uzmodelledfile = file; }
//    std::string getUzmodelledfile() { return Uzmodelledfile; }
//
//    void setUzresidualfile(std::string file) { Uzresidualfile = file; }
//    std::string getUzresidualfile() { return Uzresidualfile; }
//
//    void setMutefile(std::string file) { Mutefile = file; }
//    std::string getMutefile() { return Mutefile; }
//
//    void setApertx(T val) { apertx = val; }
//    T getApertx() { return apertx; }
//
//    void setKvp(T val) { kvp = val; }
//    T getKvp() { return kvp; }
//
//    void setKvs(T val) { kvs = val; }
//    T getKvs() { return kvs; }
//
//    void setKrho(T val) { krho = val; }
//    T getKrho() { return krho; }
//
//    void setKsource(T val) { ksource = val; }
//    T getKsource() { return ksource; }
//
//    void setVpregalpha(T val) { reg_alpha[0] = val; }
//    T getVpregalpha() { return reg_alpha[0]; }
//
//    void setVsregalpha(T val) { reg_alpha[1] = val; }
//    T getVsregalpha() { return reg_alpha[1]; }
//
//    void setVpregeps(T val) { reg_eps[0] = val; }
//    T getVpregeps() { return reg_eps[0]; }
//
//    void setVsregeps(T val) { reg_eps[1] = val; }
//    T getVsregeps() { return reg_eps[1]; }
//
//    void setSourcetype(int type) { sourcetype = type; }
//    int getSourcetype() { return sourcetype; }
//
//    // Run gradient
//    void runGrad();
//   
//    // Run BSProjection
//    void runBsproj();
//
//    // Mute gradient
//    void applyMute();
//
//    // Regularisation computation
//    void computeRegularisation(double *x);
//
//    // Combine gradients
//    void combineGradients();
//
//    // Set initial
//    int setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile);
//
//    // Save line search models
//    void saveLinesearch(double *x);
//
//    // Save Results
//    void saveResults(int i);
//
//    // Read gradient
//    void readGrad(double *g);
//
//    // Read misfit
//    void readMisfit(double *f);
//
//    void setUpdates(bool vp, bool vs, bool rho, bool source) { update_vp = vp; update_vs= vs; update_rho = rho; update_source = source; }
//    void setUpdate_vp(bool vp) { update_vp = vp;}
//    void setUpdate_vs(bool vs) { update_vs = vs;}
//    void setUpdate_rho(bool rho) { update_rho = rho;}
//    void setUpdate_source(bool source) { update_source = source;}
//
//private:
//    bool update_vp;
//    bool update_vs;
//    bool update_rho;
//    bool update_source;
//    bool srcilum;
//    bool recilum;
//    bool dataweightx;
//    bool dataweightz;
//    std::string Waveletfile;
//    std::string Vpfile;
//    std::string Vsfile;
//    std::string Rhofile;
//    std::string Vpgradfile;
//    std::string Vsgradfile;
//    std::string Dataweightxfile;
//    std::string Dataweightzfile;
//    std::string Misfitfile;
//    std::string Uxrecordfile;
//    std::string Uxmodelledfile;
//    std::string Uxresidualfile;
//    std::string Uzrecordfile;
//    std::string Uzmodelledfile;
//    std::string Uzresidualfile;
//    std::string Snapfile;
//    std::string Mutefile;
//    std::string Srcilumfile;
//    T apertx;
//    T kvp, kvs, krho, ksource;
//    T reg_eps[3];
//    T reg_alpha[3];
//    int sourcetype;
//};
//
//// ##### ELASTIC 3D WEMVA CLASS
//template<typename T>
//class WemvaElastic3D: public Wemva<T> {
//public:
//    WemvaElastic3D(); ///<Constructor
//    WemvaElastic3D(MPImodeling *_mpi); ///<Constructor
//    ~WemvaElastic3D(); ///<Destructor
//
//    void setDataweightx(bool val) { dataweightx = val; }
//    bool getDataweightx() { return dataweightx; }
//
//    void setDataweightxfile(std::string file) { Dataweightxfile = file; }
//    std::string getDataweightxfile() { return Dataweightxfile; }
//
//    void setDataweighty(bool val) { dataweighty = val; }
//    bool getDataweighty() { return dataweighty; }
//
//    void setDataweightyfile(std::string file) { Dataweightyfile = file; }
//    std::string getDataweightyfile() { return Dataweightyfile; }
//
//    void setDataweightz(bool val) { dataweightz = val; }
//    bool getDataweightz() { return dataweightz; }
//
//    void setDataweightzfile(std::string file) { Dataweightzfile = file; }
//    std::string getDataweightzfile() { return Dataweightzfile; }
//
//    void setWaveletfile(std::string file) { Waveletfile = file; }
//    std::string getWaveletfile() { return Waveletfile; }
//
//    void setVpfile(std::string file) { Vpfile = file; }
//    std::string getVpfile() { return Vpfile; }
//
//    void setVsfile(std::string file) { Vsfile = file; }
//    std::string getVsfile() { return Vsfile; }
//
//    void setRhofile(std::string file) { Rhofile = file; }
//    std::string getRhofile() { return Rhofile; }
//
//    void setVpgradfile(std::string file) { Vpgradfile = file; }
//    std::string getVpgradfile() { return Vpgradfile; }
//
//    void setVsgradfile(std::string file) { Vsgradfile = file; }
//    std::string getVsgradfile() { return Vsgradfile; }
//
//    void setMisfitfile(std::string file) { Misfitfile = file; }
//    std::string getMisfitfile() { return Misfitfile; }
//
//    void setUxrecordfile(std::string file) { Uxrecordfile = file; }
//    std::string getUxrecordfile() { return Uxrecordfile; }
//
//    void setUxmodelledfile(std::string file) { Uxmodelledfile = file; }
//    std::string getUxmodelledfile() { return Uxmodelledfile; }
//
//    void setUxresidualfile(std::string file) { Uxresidualfile = file; }
//    std::string getUxresidualfile() { return Uxresidualfile; }
//
//    void setUyrecordfile(std::string file) { Uyrecordfile = file; }
//    std::string getUyrecordfile() { return Uyrecordfile; }
//
//    void setUymodelledfile(std::string file) { Uymodelledfile = file; }
//    std::string getUymodelledfile() { return Uymodelledfile; }
//
//    void setUyresidualfile(std::string file) { Uyresidualfile = file; }
//    std::string getUyresidualfile() { return Uyresidualfile; }
//
//    void setUzrecordfile(std::string file) { Uzrecordfile = file; }
//    std::string getUzrecordfile() { return Uzrecordfile; }
//
//    void setUzmodelledfile(std::string file) { Uzmodelledfile = file; }
//    std::string getUzmodelledfile() { return Uzmodelledfile; }
//
//    void setUzresidualfile(std::string file) { Uzresidualfile = file; }
//    std::string getUzresidualfile() { return Uzresidualfile; }
//
//    void setSnapfile(std::string file) { Snapfile = file; }
//    std::string getSnapfile() { return Snapfile; }
//
//    void setMutefile(std::string file) { Mutefile = file; }
//    std::string getMutefile() { return Mutefile; }
//
//    void setApertx(T val) { apertx = val; }
//    T getApertx() { return apertx; }
//
//    void setAperty(T val) { aperty = val; }
//    T getAperty() { return aperty; }
//
//    void setKvp(T val) { kvp = val; }
//    T getKvp() { return kvp; }
//
//    void setKvs(T val) { kvs = val; }
//    T getKvs() { return kvs; }
//
//    void setKrho(T val) { krho = val; }
//    T getKrho() { return krho; }
//
//    void setKsource(T val) { ksource = val; }
//    T getKsource() { return ksource; }
//
//    void setVpregalpha(T val) { reg_alpha[0] = val; }
//    T getVpregalpha() { return reg_alpha[0]; }
//
//    void setVsregalpha(T val) { reg_alpha[1] = val; }
//    T getVsregalpha() { return reg_alpha[1]; }
//
//    void setVpregeps(T val) { reg_eps[0] = val; }
//    T getVpregeps() { return reg_eps[0]; }
//
//    void setVsregeps(T val) { reg_eps[1] = val; }
//    T getVsregeps() { return reg_eps[1]; }
//
//    void setSourcetype(int type) { sourcetype = type; }
//    int getSourcetype() { return sourcetype; }
//
//    // Run gradient
//    void runGrad();
//   
//    // Run BSProjection
//    void runBsproj();
//
//    // Mute gradient
//    void applyMute();
//
//    // Regularisation computation
//    void computeRegularisation(double *x);
//
//    // Combine gradients
//    void combineGradients();
//
//    // Set initial
//    int setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile);
//
//    // Save line search models
//    void saveLinesearch(double *x);
//
//    // Save Results
//    void saveResults(int i);
//
//    // Read gradient
//    void readGrad(double *g);
//
//    // Read misfit
//    void readMisfit(double *f);
//
//    void setUpdates(bool vp, bool vs, bool rho, bool source) { update_vp = vp; update_vs= vs; update_rho = rho; update_source = source; }
//    void setUpdate_vp(bool vp) { update_vp = vp;}
//    void setUpdate_vs(bool vs) { update_vs = vs;}
//    void setUpdate_rho(bool rho) { update_rho = rho;}
//    void setUpdate_source(bool source) { update_source = source;}
//
//private:
//    bool update_vp;
//    bool update_vs;
//    bool update_rho;
//    bool update_source;
//    bool dataweightx;
//    bool dataweighty;
//    bool dataweightz;
//    std::string Waveletfile;
//    std::string Vpfile;
//    std::string Vsfile;
//    std::string Rhofile;
//    std::string Vpgradfile;
//    std::string Vsgradfile;
//    std::string Dataweightxfile;
//    std::string Dataweightyfile;
//    std::string Dataweightzfile;
//    std::string Misfitfile;
//    std::string Uxrecordfile;
//    std::string Uxmodelledfile;
//    std::string Uxresidualfile;
//    std::string Uyrecordfile;
//    std::string Uymodelledfile;
//    std::string Uyresidualfile;
//    std::string Uzrecordfile;
//    std::string Uzmodelledfile;
//    std::string Uzresidualfile;
//    std::string Snapfile;
//    std::string Mutefile;
//    T apertx, aperty;
//    T kvp, kvs, krho, ksource;
//    T reg_eps[3];
//    T reg_alpha[3];
//    int sourcetype;
//};
//
//

}
#endif //WEMVA_H
