#ifndef LSMIGINV_H
#define LSMIGINV_H

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
#include "lsrtm.h"
#include "sort.h"
#include "bspl.h"

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

#define INV_ERR 0
#define INV_OK 1

#define LOGFILE "inversion.log"

#define PROGLOGFILE "progress.log"

#define VPLSFILE "vp_ls.rss"
#define VSLSFILE "vs_ls.rss"
#define RHOLSFILE "rho_ls.rss"
#define SOURCELSFILE "source_ls.rss"

#define VP0FILE "vp_0.rss"
#define VS0FILE "vs_0.rss"
#define RHO0FILE "rho_0.rss"
#define SOURCE0FILE "source_0.rss"

#define VPGRADFILE "vp_grad.rss"
#define VSGRADFILE "vs_grad.rss"
#define RHOGRADFILE "rho_grad.rss"
#define SOURCEGRADFILE "source_grad.rss"

#define VPREGGRADFILE "vp_reg_grad.rss"
#define VSREGGRADFILE "vs_reg_grad.rss"
#define RHOREGGRADFILE "rho_reg_grad.rss"

#define SRCILUMFILE "src_ilum.rss"

#define VPGRADCOMBFILE "vp_grad_comb.rss"
#define VSGRADCOMBFILE "vs_grad_comb.rss"
#define RHOGRADCOMBFILE "rho_grad_comb.rss"

#define VPGRADMUTEFILE "vp_grad_muted.rss"
#define VSGRADMUTEFILE "vs_grad_muted.rss"
#define RHOGRADMUTEFILE "rho_grad_muted.rss"

#define VPPROJGRADFILE "vp_proj_grad.rss"
#define VSPROJGRADFILE "vs_proj_grad.rss"
#define RHOPROJGRADFILE "rho_proj_grad.rss"

#define MISFITFILE "misfit.rss"
#define VPREGMISFITFILE "vpreg_misfit.rss"
#define VSREGMISFITFILE "vsreg_misfit.rss"
#define RHOREGMISFITFILE "rhoreg_misfit.rss"

#define PMODFILE "pmod.rss"
#define PRESFILE "pres.rss"

#define UXMODFILE "uxmod.rss"
#define UXRESFILE "uxres.rss"

#define UYMODFILE "uymod.rss"
#define UYRESFILE "uyres.rss"

#define UZMODFILE "uzmod.rss"
#define UZRESFILE "uzres.rss"

#define VP_UP "vp_up.rss"
#define VS_UP "vs_up.rss"
#define RHO_UP "rho_up.rss"
#define SOURCE_UP "source_up.rss"
#define PIMAGE_UP "pimage_up.rss"
#define SIMAGE_UP "simage_up.rss"

#define RESULTDIR "Results"



namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PAR_GRID, PAR_BSPLINE, PAR_1D} rs_paramtype; ///< Type of parameterisation

// ##### LSMIGINV CLASS
template<typename T>
class Lsmiginv {
public:
    Lsmiginv(); ///<Constructor
    Lsmiginv(MPImodeling *_mpi); ///<Constructor
    ~Lsmiginv(); ///<Destructor
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
    void setFilter(bool val) { filter = val; }
    T* getFreqs() { return freqs; }
    void setFreqs(T *pointer) { freqs = pointer; }
    rs_fwimisfit getMisfit_type() { return misfit_type; } 
    void setMisfit_type(rs_fwimisfit val) { misfit_type = val; }
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
    void setFnorm(double norm) { fnorm = norm; }
    double getFnorm() { return fnorm; }
    void normalize(double *v, double *f, int n);
    void setNoreverse(bool val) { noreverse = val; }
    bool getNoreverse() { return noreverse; }
    void setZder(bool val) { zder = val; }
    bool getZder() { return zder; }


private:
	int lpml;
	bool fs;
    bool incore;
	int order;
	int snapinc;
	int nsnaps;
    bool filter;
    T *freqs;
    double fnorm;
    std::string logfile; ///< Log file
    std::string progressfile; ///< Log file
    std::string progresslogfile; ///< Log file
    std::ofstream Flog; ///< Log stream
    bool createLog(); ///< Set name of logfile and open for writing
    bool createProglog(); ///< Set name of progress logfile and open for writing
    rs_fwimisfit misfit_type;
    rs_snapmethod snapmethod;
    rs_paramtype paramtype;
    MPImodeling *mpi;
    T dtx, dty, dtz;
    bool noreverse;
    bool zder;
};

// ##### ACOUSTIC 2D LSMIGINV CLASS
template<typename T>
class LsmiginvAcoustic2D: public Lsmiginv<T> {
public:
    LsmiginvAcoustic2D(); ///<Constructor
    LsmiginvAcoustic2D(MPImodeling *_mpi); ///<Constructor
    ~LsmiginvAcoustic2D(); ///<Destructor

    void setDataweight(bool val) { dataweight = val; }
    bool getDataweight() { return dataweight; }

    void setDataweightfile(std::string file) { Dataweightfile = file; }
    std::string getDataweightfile() { return Dataweightfile; }

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setRhofile(std::string file) { Rhofile = file; }
    std::string getRhofile() { return Rhofile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setSrcilumfile(std::string file) { Srcilumfile = file; }
    std::string getSrcilumfile() { return Srcilumfile; }

    void setPimagefile(std::string file) { Pimagefile = file; }
    std::string getPimagefile() { return Pimagefile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setPsnapfile(std::string file) { Psnapfile = file; }
    std::string getPsnapfile() { return Psnapfile; }

    void setPrecordfile(std::string file) { Precordfile = file; }
    std::string getPrecordfile() { return Precordfile; }

    void setMultiplefile(std::string file) { Multiplefile = file; }
    std::string getMultiplefile() { return Multiplefile; }

    void setPmodelledfile(std::string file) { Pmodelledfile = file; }
    std::string getPmodelledfile() { return Pmodelledfile; }

    void setMutefile(std::string file) { Mutefile = file; }
    std::string getMutefile() { return Mutefile; }

    void setPresidualfile(std::string file) { Presidualfile = file; }
    std::string getPresidualfile() { return Presidualfile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setVpregalpha(T val) { reg_alpha = val; }
    T getVpregalpha() { return reg_alpha; }

    void setVpregeps(T val) { reg_eps = val; }
    T getVpregeps() { return reg_eps; }

    // Run gradient
    void runGrad();
    void runGrad_Multiples();
   
    // Run BSProjection
    void runBsproj();

    // Mute gradient
    void applyMute();

    // Correct for source ilumination
    void applySrcilum();

    // Regularisation computation
    void computeRegularisation(double *x);
    void computeTikhonovRegularisation(double *x);

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

    void setUpdate_vp(bool vp) { update_vp = vp;}
    void setSrcilum(bool val) { srcilumset = val;}
    bool getSrcilum() { return srcilumset;}


private:
    bool update_vp;
    bool dataweight;
    bool srcilumset;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Pimagefile;
    std::string Misfitfile;
    std::string Psnapfile;
    std::string Precordfile;
    std::string Multiplefile;
    std::string Pmodelledfile;
    std::string Presidualfile;
    std::string Dataweightfile;
    std::string Mutefile;
    std::string Srcilumfile;
    T apertx;
    T kvp;
    T reg_eps;
    T reg_alpha;
};

// ##### ACOUSTIC 3D LSMIGINV CLASS
template<typename T>
class LsmiginvAcoustic3D: public Lsmiginv<T> {
public:
    LsmiginvAcoustic3D(); ///<Constructor
    LsmiginvAcoustic3D(MPImodeling *_mpi); ///<Constructor
    ~LsmiginvAcoustic3D(); ///<Destructor

    void setDataweight(bool val) { dataweight = val; }
    bool getDataweight() { return dataweight; }

    void setDataweightfile(std::string file) { Dataweightfile = file; }
    std::string getDataweightfile() { return Dataweightfile; }

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setRhofile(std::string file) { Rhofile = file; }
    std::string getRhofile() { return Rhofile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setRhogradfile(std::string file) { Rhogradfile = file; }
    std::string getRhogradfile() { return Rhogradfile; }

    void setWavgradfile(std::string file) { Wavgradfile = file; }
    std::string getWavgradfile() { return Wavgradfile; }

    void setSrcilumfile(std::string file) { Srcilumfile = file; }
    std::string getSrcilumfile() { return Srcilumfile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setPsnapfile(std::string file) { Psnapfile = file; }
    std::string getPsnapfile() { return Psnapfile; }

    void setPrecordfile(std::string file) { Precordfile = file; }
    std::string getPrecordfile() { return Precordfile; }

    void setPmodelledfile(std::string file) { Pmodelledfile = file; }
    std::string getPmodelledfile() { return Pmodelledfile; }

    void setMutefile(std::string file) { Mutefile = file; }
    std::string getMutefile() { return Mutefile; }

    void setPresidualfile(std::string file) { Presidualfile = file; }
    std::string getPresidualfile() { return Presidualfile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setAperty(T val) { aperty = val; }
    T getAperty() { return aperty; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setKrho(T val) { krho = val; }
    T getKrho() { return krho; }

    void setKsource(T val) { ksource = val; }
    T getKsource() { return ksource; }

    void setVpregalpha(T val) { reg_alpha[0] = val; }
    T getVpregalpha() { return reg_alpha[0]; }

    void setRhoregalpha(T val) { reg_alpha[1] = val; }
    T getRhoregalpha() { return reg_alpha[1]; }

    void setVpregeps(T val) { reg_eps[0] = val; }
    T getVpregeps() { return reg_eps[0]; }

    void setRhoregeps(T val) { reg_eps[1] = val; }
    T getRhoregeps() { return reg_eps[1]; }

    // Run gradient
    void runGrad();
   
    // Run BSProjection
    void runBsproj();

    // Mute gradient
    void applyMute();

    // Correct for source ilumination
    void applySrcilum();

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

    void setUpdates(bool vp, bool rho, bool source) { update_vp = vp; update_rho = rho; update_source = source; }
    void setUpdate_vp(bool vp) { update_vp = vp;}
    void setUpdate_rho(bool rho) { update_rho = rho;}
    void setUpdate_source(bool source) { update_source = source;}
    void setSrcilum(bool val) { srcilumset = val;}
    bool getSrcilum() { return srcilumset;}

private:
    bool update_vp;
    bool update_rho;
    bool update_source;
    bool dataweight;
    bool srcilumset;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightfile;
    std::string Misfitfile;
    std::string Psnapfile;
    std::string Precordfile;
    std::string Pmodelledfile;
    std::string Presidualfile;
    std::string Mutefile;
    std::string Srcilumfile;
    T apertx, aperty;
    T kvp, krho, ksource;
    T reg_eps[2];
    T reg_alpha[2];
};

// ##### ELASTIC 2D LSMIGINV CLASS
template<typename T>
class LsmiginvElastic2D: public Lsmiginv<T> {
public:
    LsmiginvElastic2D(); ///<Constructor
    LsmiginvElastic2D(MPImodeling *_mpi); ///<Constructor
    ~LsmiginvElastic2D(); ///<Destructor

    void setDataweightx(bool val) { dataweightx = val; }
    bool getDataweightx() { return dataweightx; }

    void setDataweightz(bool val) { dataweightz = val; }
    bool getDataweightz() { return dataweightz; }

    void setDataweightxfile(std::string file) { Dataweightxfile = file; }
    std::string getDataweightxfile() { return Dataweightxfile; }

    void setDataweightzfile(std::string file) { Dataweightzfile = file; }
    std::string getDataweightzfile() { return Dataweightzfile; }

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setSrcilumfile(std::string file) { Srcilumfile = file; }
    std::string getSrcilumfile() { return Srcilumfile; }

    void setVsfile(std::string file) { Vsfile = file; }
    std::string getVsfile() { return Vsfile; }

    void setRhofile(std::string file) { Rhofile = file; }
    std::string getRhofile() { return Rhofile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setVsgradfile(std::string file) { Vsgradfile = file; }
    std::string getVsgradfile() { return Vsgradfile; }

    void setRhogradfile(std::string file) { Rhogradfile = file; }
    std::string getRhogradfile() { return Rhogradfile; }

    void setWavgradfile(std::string file) { Wavgradfile = file; }
    std::string getWavgradfile() { return Wavgradfile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setUxrecordfile(std::string file) { Uxrecordfile = file; }
    std::string getUxrecordfile() { return Uxrecordfile; }

    void setUxmodelledfile(std::string file) { Uxmodelledfile = file; }
    std::string getUxmodelledfile() { return Uxmodelledfile; }

    void setUxresidualfile(std::string file) { Uxresidualfile = file; }
    std::string getUxresidualfile() { return Uxresidualfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setUzrecordfile(std::string file) { Uzrecordfile = file; }
    std::string getUzrecordfile() { return Uzrecordfile; }

    void setUzmodelledfile(std::string file) { Uzmodelledfile = file; }
    std::string getUzmodelledfile() { return Uzmodelledfile; }

    void setUzresidualfile(std::string file) { Uzresidualfile = file; }
    std::string getUzresidualfile() { return Uzresidualfile; }

    void setMutefile(std::string file) { Mutefile = file; }
    std::string getMutefile() { return Mutefile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setKvs(T val) { kvs = val; }
    T getKvs() { return kvs; }

    void setKrho(T val) { krho = val; }
    T getKrho() { return krho; }

    void setKsource(T val) { ksource = val; }
    T getKsource() { return ksource; }

    void setVpregalpha(T val) { reg_alpha[0] = val; }
    T getVpregalpha() { return reg_alpha[0]; }

    void setVsregalpha(T val) { reg_alpha[1] = val; }
    T getVsregalpha() { return reg_alpha[1]; }

    void setRhoregalpha(T val) { reg_alpha[2] = val; }
    T getRhoregalpha() { return reg_alpha[2]; }

    void setVpregeps(T val) { reg_eps[0] = val; }
    T getVpregeps() { return reg_eps[0]; }

    void setVsregeps(T val) { reg_eps[1] = val; }
    T getVsregeps() { return reg_eps[1]; }

    void setRhoregeps(T val) { reg_eps[2] = val; }
    T getRhoregeps() { return reg_eps[2]; }

    void setSourcetype(int type) { sourcetype = type; }
    int getSourcetype() { return sourcetype; }

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
    int setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile);

    // Save line search models
    void saveLinesearch(double *x);

    // Save Results
    void saveResults(int i);

    // Read gradient
    void readGrad(double *g);

    // Read misfit
    void readMisfit(double *f);

    void setUpdates(bool vp, bool vs, bool rho, bool source) { update_vp = vp; update_vs= vs; update_rho = rho; update_source = source; }
    void setUpdate_vp(bool vp) { update_vp = vp;}
    void setUpdate_vs(bool vs) { update_vs = vs;}
    void setUpdate_rho(bool rho) { update_rho = rho;}
    void setUpdate_source(bool source) { update_source = source;}

private:
    bool update_vp;
    bool update_vs;
    bool update_rho;
    bool update_source;
    bool srcilum;
    bool recilum;
    bool dataweightx;
    bool dataweightz;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Vsgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightxfile;
    std::string Dataweightzfile;
    std::string Misfitfile;
    std::string Uxrecordfile;
    std::string Uxmodelledfile;
    std::string Uxresidualfile;
    std::string Uzrecordfile;
    std::string Uzmodelledfile;
    std::string Uzresidualfile;
    std::string Snapfile;
    std::string Mutefile;
    std::string Srcilumfile;
    T apertx;
    T kvp, kvs, krho, ksource;
    T reg_eps[3];
    T reg_alpha[3];
    int sourcetype;
};

// ##### ELASTIC 3D LSMIGINV CLASS
template<typename T>
class LsmiginvElastic3D: public Lsmiginv<T> {
public:
    LsmiginvElastic3D(); ///<Constructor
    LsmiginvElastic3D(MPImodeling *_mpi); ///<Constructor
    ~LsmiginvElastic3D(); ///<Destructor

    void setDataweightx(bool val) { dataweightx = val; }
    bool getDataweightx() { return dataweightx; }

    void setDataweightxfile(std::string file) { Dataweightxfile = file; }
    std::string getDataweightxfile() { return Dataweightxfile; }

    void setDataweighty(bool val) { dataweighty = val; }
    bool getDataweighty() { return dataweighty; }

    void setDataweightyfile(std::string file) { Dataweightyfile = file; }
    std::string getDataweightyfile() { return Dataweightyfile; }

    void setDataweightz(bool val) { dataweightz = val; }
    bool getDataweightz() { return dataweightz; }

    void setDataweightzfile(std::string file) { Dataweightzfile = file; }
    std::string getDataweightzfile() { return Dataweightzfile; }

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setVsfile(std::string file) { Vsfile = file; }
    std::string getVsfile() { return Vsfile; }

    void setRhofile(std::string file) { Rhofile = file; }
    std::string getRhofile() { return Rhofile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setVsgradfile(std::string file) { Vsgradfile = file; }
    std::string getVsgradfile() { return Vsgradfile; }

    void setRhogradfile(std::string file) { Rhogradfile = file; }
    std::string getRhogradfile() { return Rhogradfile; }

    void setWavgradfile(std::string file) { Wavgradfile = file; }
    std::string getWavgradfile() { return Wavgradfile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setUxrecordfile(std::string file) { Uxrecordfile = file; }
    std::string getUxrecordfile() { return Uxrecordfile; }

    void setUxmodelledfile(std::string file) { Uxmodelledfile = file; }
    std::string getUxmodelledfile() { return Uxmodelledfile; }

    void setUxresidualfile(std::string file) { Uxresidualfile = file; }
    std::string getUxresidualfile() { return Uxresidualfile; }

    void setUyrecordfile(std::string file) { Uyrecordfile = file; }
    std::string getUyrecordfile() { return Uyrecordfile; }

    void setUymodelledfile(std::string file) { Uymodelledfile = file; }
    std::string getUymodelledfile() { return Uymodelledfile; }

    void setUyresidualfile(std::string file) { Uyresidualfile = file; }
    std::string getUyresidualfile() { return Uyresidualfile; }

    void setUzrecordfile(std::string file) { Uzrecordfile = file; }
    std::string getUzrecordfile() { return Uzrecordfile; }

    void setUzmodelledfile(std::string file) { Uzmodelledfile = file; }
    std::string getUzmodelledfile() { return Uzmodelledfile; }

    void setUzresidualfile(std::string file) { Uzresidualfile = file; }
    std::string getUzresidualfile() { return Uzresidualfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setMutefile(std::string file) { Mutefile = file; }
    std::string getMutefile() { return Mutefile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setAperty(T val) { aperty = val; }
    T getAperty() { return aperty; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setKvs(T val) { kvs = val; }
    T getKvs() { return kvs; }

    void setKrho(T val) { krho = val; }
    T getKrho() { return krho; }

    void setKsource(T val) { ksource = val; }
    T getKsource() { return ksource; }

    void setVpregalpha(T val) { reg_alpha[0] = val; }
    T getVpregalpha() { return reg_alpha[0]; }

    void setVsregalpha(T val) { reg_alpha[1] = val; }
    T getVsregalpha() { return reg_alpha[1]; }

    void setRhoregalpha(T val) { reg_alpha[2] = val; }
    T getRhoregalpha() { return reg_alpha[2]; }

    void setVpregeps(T val) { reg_eps[0] = val; }
    T getVpregeps() { return reg_eps[0]; }

    void setVsregeps(T val) { reg_eps[1] = val; }
    T getVsregeps() { return reg_eps[1]; }

    void setRhoregeps(T val) { reg_eps[2] = val; }
    T getRhoregeps() { return reg_eps[2]; }

    void setSourcetype(int type) { sourcetype = type; }
    int getSourcetype() { return sourcetype; }

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
    int setInitial(double *x, std::string vpfile, std::string vsfile, std::string rhofile, std::string sourcefile);

    // Save line search models
    void saveLinesearch(double *x);

    // Save Results
    void saveResults(int i);

    // Read gradient
    void readGrad(double *g);

    // Read misfit
    void readMisfit(double *f);

    void setUpdates(bool vp, bool vs, bool rho, bool source) { update_vp = vp; update_vs= vs; update_rho = rho; update_source = source; }
    void setUpdate_vp(bool vp) { update_vp = vp;}
    void setUpdate_vs(bool vs) { update_vs = vs;}
    void setUpdate_rho(bool rho) { update_rho = rho;}
    void setUpdate_source(bool source) { update_source = source;}

private:
    bool update_vp;
    bool update_vs;
    bool update_rho;
    bool update_source;
    bool dataweightx;
    bool dataweighty;
    bool dataweightz;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Vsgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightxfile;
    std::string Dataweightyfile;
    std::string Dataweightzfile;
    std::string Misfitfile;
    std::string Uxrecordfile;
    std::string Uxmodelledfile;
    std::string Uxresidualfile;
    std::string Uyrecordfile;
    std::string Uymodelledfile;
    std::string Uyresidualfile;
    std::string Uzrecordfile;
    std::string Uzmodelledfile;
    std::string Uzresidualfile;
    std::string Snapfile;
    std::string Mutefile;
    T apertx, aperty;
    T kvp, kvs, krho, ksource;
    T reg_eps[3];
    T reg_alpha[3];
    int sourcetype;
};



}
#endif //LSMIGINV_H
