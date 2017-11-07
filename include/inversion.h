#ifndef INVERSION_H
#define INVERSION_H

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
#include "bspl.h"
#include "interp.h"
#include "parallel.h"
#include "fwi.h"
#include "sort.h"
#include <inparse.h>

#define RUN_F_GRAD 0
#define RUN_BS_PROJ 1
#define BREAK_LOOP 2

namespace rockseis {

// ##### INVERSION CLASS
template<typename T>
class Inversion {
public:
    Inversion(); ///<Constructor
    Inversion(MPImodeling *_mpi); ///<Constructor
    ~Inversion(); ///<Destructor
    void setMpi(MPImodeling *_mpi) { mpi = _mpi; }
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




private:
	int lpml;
	bool fs;
    bool incore;
	int order;
	int snapinc;
	int nsnaps;
    rs_fwimisfit misfit_type;
    rs_snapmethod snapmethod;
    MPImodeling *mpi;
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

    void setPresidualfile(std::string file) { Presidualfile = file; }
    std::string getPresidualfile() { return Presidualfile; }

    void setApertx(T val) { apertx = val; }
    T getApertx() { return apertx; }

    // Run gradient
    void runAcousticfwigrad2d();


private:
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
    T apertx;
};

}
#endif //INVERSION_H
