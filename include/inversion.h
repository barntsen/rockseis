#ifndef INVERSION_H
#define INVERSION_H

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
#include "fwi.h"
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

#define VXMODFILE "vxmod.rss"
#define VXRESFILE "vxres.rss"

#define VYMODFILE "vymod.rss"
#define VYRESFILE "vyres.rss"

#define VZMODFILE "vzmod.rss"
#define VZRESFILE "vzres.rss"

#define VP_UP "vp_up.rss"
#define VS_UP "vs_up.rss"
#define RHO_UP "rho_up.rss"
#define SOURCE_UP "source_up.rss"

#define RESULTDIR "Results"


namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PAR_GRID, PAR_BSPLINE, PAR_1D} rs_paramtype; ///< Type of parameterisation

// ##### INVERSION CLASS
template<typename T>
class Inversion {
public:
    Inversion(); ///<Constructor
    Inversion(MPImodeling *_mpi); ///<Constructor
    ~Inversion(); ///<Destructor
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
    void setIncore(bool val) { incore = val; } 
    bool getIncore() { return incore; }
    int getOrder() { return order; }
    void setOrder(int val) { order = val; }
    int getSnapinc() { return snapinc; }
    void setSnapinc(int val) { snapinc = val; }
    int getNsnaps() { return nsnaps; }
    void setNsnaps(int val) { nsnaps = val; }
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

private:
	int lpml;
	bool fs;
    bool incore;
	int order;
	int snapinc;
	int nsnaps;
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
};

// ##### ACOUSTIC 2D INVERSION CLASS
template<typename T>
class InversionAcoustic2D: public Inversion<T> {
public:
    InversionAcoustic2D(); ///<Constructor
    InversionAcoustic2D(MPImodeling *_mpi); ///<Constructor
    ~InversionAcoustic2D(); ///<Destructor

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

private:
    bool update_vp;
    bool update_rho;
    bool update_source;
    bool dataweight;
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
    T apertx;
    T kvp, krho, ksource;
    T reg_eps[2];
    T reg_alpha[2];
};

// ##### ELASTIC 2D INVERSION CLASS
template<typename T>
class InversionElastic2D: public Inversion<T> {
public:
    InversionElastic2D(); ///<Constructor
    InversionElastic2D(MPImodeling *_mpi); ///<Constructor
    ~InversionElastic2D(); ///<Destructor

    void setDataweight(bool val) { dataweight = val; }
    bool getDataweight() { return dataweight; }

    void setDataweightfile(std::string file) { Dataweightfile = file; }
    std::string getDataweightfile() { return Dataweightfile; }

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

    void setVxrecordfile(std::string file) { Vxrecordfile = file; }
    std::string getVxrecordfile() { return Vxrecordfile; }

    void setVxmodelledfile(std::string file) { Vxmodelledfile = file; }
    std::string getVxmodelledfile() { return Vxmodelledfile; }

    void setVxresidualfile(std::string file) { Vxresidualfile = file; }
    std::string getVxresidualfile() { return Vxresidualfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setVzrecordfile(std::string file) { Vzrecordfile = file; }
    std::string getVzrecordfile() { return Vzrecordfile; }

    void setVzmodelledfile(std::string file) { Vzmodelledfile = file; }
    std::string getVzmodelledfile() { return Vzmodelledfile; }

    void setVzresidualfile(std::string file) { Vzresidualfile = file; }
    std::string getVzresidualfile() { return Vzresidualfile; }

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
    bool dataweight;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Vsgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightfile;
    std::string Misfitfile;
    std::string Vxrecordfile;
    std::string Vxmodelledfile;
    std::string Vxresidualfile;
    std::string Vzrecordfile;
    std::string Vzmodelledfile;
    std::string Vzresidualfile;
    std::string Snapfile;
    std::string Mutefile;
    T apertx;
    T kvp, kvs, krho, ksource;
    T reg_eps[3];
    T reg_alpha[3];
    int sourcetype;
};

}
#endif //INVERSION_H
