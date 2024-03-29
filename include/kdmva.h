#ifndef KDMVA_H
#define KDMVA_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "geometry.h"
#include "utils.h"
#include "file.h"
#include "model.h"
#include "data.h"
#include "der.h"
#include "interp.h"
#include "parallel.h"
#include "kdmig.h"
#include "ttable.h"
#include "sort.h"
#include "bspl.h"
#include "hilbert.h"

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

#define KVA_ERR 0
#define KVA_OK 1

#define LOGFILE "kdmva.log"

#define PROGLOGFILE "progress.log"

#define VPLSFILE "vp_ls.rss"
#define VSLSFILE "vs_ls.rss"
#define SOURCELSFILE "source_ls.rss"

#define TTABLELSFILE "ttable_ls.rss"

#define VP0FILE "vp_0.rss"
#define VS0FILE "vs_0.rss"
#define SOURCE0FILE "source_0.rss"

#define PIMAGEFILE "pimage.rss"
#define SIMAGEFILE "simage.rss"

#define PIMAGERESFILE "pimage_res.rss"
#define SIMAGERESFILE "simage_res.rss"

#define VPGRADFILE "vp_grad.rss"
#define VSGRADFILE "vs_grad.rss"

#define VPREGGRADFILE "vp_reg_grad.rss"
#define VSREGGRADFILE "vs_reg_grad.rss"

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

#define PIMAGE_UP "pimage_up.rss"
#define SIMAGE_UP "simage_up.rss"

#define RESULTDIR "Results"

#define THRES 70

namespace rockseis {

// =============== ENUMS =============== //
typedef enum {PAR_GRID, PAR_BSPLINE, PAR_1D, PAR_AVG} rs_paramtype; ///< Type of parameterisation

// ##### KDMVA CLASS
template<typename T>
class Kdmva {
public:
    Kdmva(); ///<Constructor
    Kdmva(MPImodeling *_mpi); ///<Constructor
    ~Kdmva(); ///<Destructor
    void createResult();
    double vector_norm(double *v, const int type, const int n);
    void find_max(T *v, T *max, int *imax, const int n);
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
    int getSouinc() { return souinc; }
    void setSouinc(int val) { souinc = val; }
    int getRecinc() { return recinc; }
    void setRecinc(int val) { recinc = val; }
    int getSnapinc() { return snapinc; }
    void setSnapinc(int val) { snapinc = val; }
    int getNsnaps() { return nsnaps; }
    void setNsnaps(int val) { nsnaps = val; }
    rs_paramtype getParamtype() { return paramtype; } 
    void setParamtype(rs_paramtype val) { paramtype = val; }
    rs_wemvamisfit getMisfit_type() { return misfit_type; } 
    void setMisfit_type(rs_wemvamisfit val) { misfit_type = val; }
    T getRadius() { return radius; }
    void setRadius(T val) { radius = val; }
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
    bool getConstrain() { return constrain; }
    void setConstrain(bool val) { constrain = val; }

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
    rs_paramtype paramtype;
    rs_wemvamisfit misfit_type;
    MPImodeling *mpi;
    T dtx, dty, dtz;
    int nhx,nhy,nhz;
    int souinc;
    int recinc;
    T radius;
    bool constrain;
};

// ##### ACOUSTIC 2D KDMVA CLASS
template<typename T>
class KdmvaAcoustic2D: public Kdmva<T> {
public:
    KdmvaAcoustic2D(); ///<Constructor
    KdmvaAcoustic2D(MPImodeling *_mpi); ///<Constructor
    ~KdmvaAcoustic2D(); ///<Destructor

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setTtablefile(std::string file) { Ttablefile = file; }
    std::string getTtablefile() { return Ttablefile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setPimagefile(std::string file) { Pimagefile = file; }
    std::string getPimagefile() { return Pimagefile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setPrecordfile(std::string file) { Precordfile = file; }
    std::string getPrecordfile() { return Precordfile; }

    void setModelmutefile(std::string file) { Modelmutefile = file; }
    std::string getModelmutefile() { return Modelmutefile; }

    void setResidualmutefile(std::string file) { Residualmutefile = file; }
    std::string getResidualmutefile() { return Residualmutefile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setLboundfile(std::string file) { Lboundfile = file; }
    std::string getLboundfile() { return Lboundfile; }

    void setUboundfile(std::string file) { Uboundfile = file; }
    std::string getUboundfile() { return Uboundfile; }

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

    // Apply Chain rule for constrained optimisation using logistic model
    void applyChainrule(double *x);

    // Regularisation computation
    void computeRegularisation(double *x);

    // Compute Misfit
    void computeMisfit(std::shared_ptr<rockseis::Image2D<T>> pimage);

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
    std::string Ttablefile;
    std::string Vpgradfile;
    std::string Pimagefile;
    std::string Misfitfile;
    std::string Snapfile;
    std::string Precordfile;
    std::string Modelmutefile;
    std::string Residualmutefile;
    std::string Lboundfile;
    std::string Uboundfile;
    T apertx;
    T kvp;
    T reg_eps[2];
    T reg_alpha[2];
};

// ##### ELASTIC 2D KDMVA CLASS
template<typename T>
class KdmvaElastic2D: public Kdmva<T> {
public:
    KdmvaElastic2D(); ///<Constructor
    KdmvaElastic2D(MPImodeling *_mpi); ///<Constructor
    ~KdmvaElastic2D(); ///<Destructor

    void setWaveletfile(std::string file) { Waveletfile = file; }
    std::string getWaveletfile() { return Waveletfile; }

    void setVpfile(std::string file) { Vpfile = file; }
    std::string getVpfile() { return Vpfile; }

    void setVsfile(std::string file) { Vsfile = file; }
    std::string getVsfile() { return Vsfile; }

    void setVpgradfile(std::string file) { Vpgradfile = file; }
    std::string getVpgradfile() { return Vpgradfile; }

    void setVsgradfile(std::string file) { Vsgradfile = file; }
    std::string getVsgradfile() { return Vsgradfile; }

    void setPimagefile(std::string file) { Pimagefile = file; }
    std::string getPimagefile() { return Pimagefile; }

    void setSimagefile(std::string file) { Simagefile = file; }
    std::string getSimagefile() { return Simagefile; }

    void setMisfitfile(std::string file) { Misfitfile = file; }
    std::string getMisfitfile() { return Misfitfile; }

    void setUxrecordfile(std::string file) { Uxrecordfile = file; }
    std::string getUxrecordfile() { return Uxrecordfile; }

    void setSnapfile(std::string file) { Snapfile = file; }
    std::string getSnapfile() { return Snapfile; }

    void setUzrecordfile(std::string file) { Uzrecordfile = file; }
    std::string getUzrecordfile() { return Uzrecordfile; }

    void setModelmutefile(std::string file) { Modelmutefile = file; }
    std::string getModelmutefile() { return Modelmutefile; }

    void setResidualmutefile(std::string file) { Residualmutefile = file; }
    std::string getResidualmutefile() { return Residualmutefile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    void setKvp(T val) { kvp = val; }
    T getKvp() { return kvp; }

    void setKvs(T val) { kvs = val; }
    T getKvs() { return kvs; }

    void setVpregalpha(T val) { reg_alpha[0] = val; }
    T getVpregalpha() { return reg_alpha[0]; }

    void setVsregalpha(T val) { reg_alpha[1] = val; }
    T getVsregalpha() { return reg_alpha[1]; }

    void setVpregeps(T val) { reg_eps[0] = val; }
    T getVpregeps() { return reg_eps[0]; }

    void setVsregeps(T val) { reg_eps[1] = val; }
    T getVsregeps() { return reg_eps[1]; }

    void setSourcetype(int type) { sourcetype = type; }
    int getSourcetype() { return sourcetype; }

    void setWavemode(int mode);
    int getWavemode() { return wavemode; }

    // Run PP gradient
    void runPPgrad();

    // Run PS gradient
    void runPSgrad();

    // Run gradient
    void runGrad() { if(wavemode == 1) runPSgrad(); else runPPgrad(); } 
   
    // Run BSProjection
    void runBsproj();

    // Mute gradient
    void applyMute();

    // Regularisation computation
    void computeRegularisation(double *x);

    // Combine gradients
    void combineGradients();

    // Compute Misfit
    void computeMisfit(std::shared_ptr<rockseis::Image2D<T>> pimage, std::string imageresfile);

    // Set initial
    int setInitial(double *x, std::string vpfile, std::string vsfile);

    // Save line search models
    void saveLinesearch(double *x);

    // Save Results
    void saveResults(int i);

    // Read gradient
    void readGrad(double *g);

    // Read misfit
    void readMisfit(double *f);

    void setUpdate_vp(bool vp) { update_vp = vp;}
    void setUpdate_vs(bool vs) { update_vs = vs;}

private:
    bool update_vp;
    bool update_vs;
    bool srcilum;
    bool recilum;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Vpgradfile;
    std::string Vsgradfile;
    std::string Misfitfile;
    std::string Uxrecordfile;
    std::string Uzrecordfile;
    std::string Snapfile;
    std::string Modelmutefile;
    std::string Residualmutefile;
    std::string Pimagefile;
    std::string Simagefile;
    T apertx;
    T kvp, kvs;
    T reg_eps[2];
    T reg_alpha[2];
    int sourcetype;
    int wavemode;
};

}
#endif //KDMVA_H
