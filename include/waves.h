#ifndef WAVES_H
#define WAVES_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "model.h"
#include "pml.h"
#include "der.h"
#include "geometry.h"
#include "utils.h"
#include "data.h"
#include "file.h"
#include "interp.h"

#define WAVES_OK 1;
#define WAVES_ERR 0;

#define LANC_SIZE 3
#define LANC(x,a) (this->sinc(x)*this->sinc((x)/a))

namespace rockseis {

/** The abstract waves class
 *
 */
template<typename T>
class Waves {
public:
    Waves();		///< Constructor
    Waves(const int _dim, const int _nx, const int _ny, const int _nz, const int _nt, const int _lpml, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    virtual ~Waves();	///< Destructor
    
    // Get functions
    int getDim() { return dim; }		///< Get dimension
    int getNx() { return geometry->getN(1); }	///< Get Nx
    int getNy() { return geometry->getN(2); }	///< Get Ny
    int getNz() { return geometry->getN(3); }	///< Get Nz
    int getNx_pml();	///< Get Nx padded
    int getNy_pml();	///< Get Ny padded
    int getNz_pml();	///< Get Nz padded
    int getNt() { return geometry->getN(4); }	///< Get Nt
    int getLpml() { return lpml; }		///< Get lpml
    T getDx() { return geometry->getD(1); }	///< Get Dx
    T getDy() { return geometry->getD(2); }	///< Get Dy
    T getDz() { return geometry->getD(3); }	///< Get Dz
    T getDt() { return geometry->getD(4); }	///< Get Dt
    T getOx() { return geometry->getO(1); }	///< Get Ox
    T getOy() { return geometry->getO(2); }	///< Get Oy
    T getOz() { return geometry->getO(3); }	///< Get Oz
    T getOt() { return geometry->getO(4); }	///< Get Ot

    bool getDomdec() { return domdec; } ///< Get domain decomposition flag

    //Interpolation function
    T sinc(T x) { if(x == 0) return 1; else return std::sin(PI*x)/(PI*x); } 

        // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }///< Set Nz
    void setNt(const int _nt) { geometry->setN(4, _nt); }///< Set Nt
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set lpml
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setDt(const T _dt) { geometry->setD(4, _dt); }	///< Set Dt
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setOt(const T _ot) { geometry->setO(4, _ot); }	///< Set Ot
    void setDim(const int _dim) { dim = _dim; }		///< Set the dimension
    void setDomain(std::shared_ptr<Domain<T>> dom) { domain = dom; if(dom->getStatus()) domdec = true;} ///<Set domain decomposition obj

private:
    int dim;
    std::shared_ptr<Geometry<T>> geometry; // regular geometry
    std::shared_ptr<Domain<T>> domain; // Domain pointer 
    int lpml; // PML boundary size
    bool domdec; // Domain decomposition flag
};


/** The 2D Acoustic WAVES class
 *
 */
template<typename T>
class WavesAcoustic2D: public Waves<T> {
public:
    WavesAcoustic2D();					///< Constructor
    WavesAcoustic2D(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesAcoustic2D(std::shared_ptr<rockseis::ModelAcoustic2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesAcoustic2D();					///< Destructor
    
    // Time stepping
    void forwardstepVelocity(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with velocity
    void forwardstepStress(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with stress

    // Get functions
    std::shared_ptr<PmlAcoustic2D<T>> getPml() { return Pml; } ///< Get pml
    T * getP() { return P; }  ///< Get advanced pressure
    T * getVx() { return Vx; }  ///< Get x-component of velocity 
    T * getVz() { return Vz; }  ///< Get z-component of velocity

    void computeABC() { Pml->callcompABC(); } ///< Compute the PML decay constants (needed if changes are made to the Amax, Kmax and Smax)

    // Insert source functions
    void insertSource(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Velocity type or Pressure )

    // Record data at receivers functions
    void recordData(std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Velocity type or Pressure )

private:
    T *P; // Pressure at time t+1
    T *Vx; // Velocity component at time t
    T *Vz; // Velocity component at time t
    std::shared_ptr<PmlAcoustic2D<T>> Pml; // Associated PML class
    
};

/** The 3D Acoustic WAVES class
 *
 */
template<typename T>
class WavesAcoustic3D: public Waves<T> {
public:
    WavesAcoustic3D();	///< Constructor
    WavesAcoustic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _L, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    WavesAcoustic3D(std::shared_ptr<rockseis::ModelAcoustic3D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesAcoustic3D();	///< Destructor

    // Get functions
    std::shared_ptr<PmlAcoustic3D<T>> getPml() { return Pml; }  ///< Get pml
    T * getP() { return P; }  ///< Get advanced pressure
    T * getVx() { return Vx; }  ///< Get current x-velocity
    T * getVy() { return Vy; }  ///< Get current y-velocity
    T * getVz() { return Vz; }  ///< Get current z-velocity

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with velocity
    void forwardstepStress(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with stress
    
    // Insert source functions
    void insertSource(std::shared_ptr<ModelAcoustic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Velocity type or Pressure )

    // Record data at receivers functions
    void recordData(std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Velocity type or Pressure )

    void computeABC() { Pml->callcompABC(); } ///< Compute the PML decay constants (needed if changes are made to the Amax, Kmax and Smax)

private:
    T *P; // Pressure at time t+1
    T *Vx; // Velocity component at time t
    T *Vy; // Velocity component at time t
    T *Vz; // Velocity component at time t
    std::shared_ptr<PmlAcoustic3D<T>> Pml; // Associated Pml class

};

/** The 2D Elastic WAVES class
 *
 */
template<typename T>
class WavesElastic2D: public Waves<T> {
public:
    WavesElastic2D();	///< Constructor
    ~WavesElastic2D();	///< Destructor
    WavesElastic2D(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesElastic2D(std::shared_ptr<rockseis::ModelElastic2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor

    // Get functions
    std::shared_ptr<PmlElastic2D<T>> getPml() { return Pml; } ///< Get Pml 
    T * getSxx() { return Sxx; }  ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }  ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }  ///< Get Stress component at time t+1

    T * getVx() { return Vx; }  ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }  ///< Get Velocity component at time t+1/2

    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with Stress
    void backwardstepVelocity(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step backward with particle velocity 
    void backwardstepStress(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with Stress

    // Insert source functions
    void insertSource(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling 

    // Record data at receivers functions
    void recordData(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling 



private:
    T *Sxx;  // Stress component at time t+1
    T *Szz; // Stress component at time t+1
    T *Sxz; // Stress component at time t+1
    T *Vx; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic2D<T>> Pml; // Associated Pml class
    bool adjoint; 

};

/** The 3D Elastic WAVES class
 *
 */
template<typename T>
class WavesElastic3D: public Waves<T> {
public:
    WavesElastic3D();	///< Constructor
    WavesElastic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _L, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    WavesElastic3D(std::shared_ptr<rockseis::ModelElastic3D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesElastic3D();	///< Destructor

    // Get functions
    std::shared_ptr<PmlElastic3D<T>> getPml() { return Pml; }
    T * getSxx() { return Sxx; }    ///< Get Stress component at time t+1
    T * getSyy() { return Syy; }    ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }    ///< Get Stress component at time t+1
    T * getSyz() { return Syz; }    ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }    ///< Get Stress component at time t+1
    T * getSxy() { return Sxy; }    ///< Get Stress component at time t+1
    T * getVx() { return Vx; }    ///< Get Velocity component at time t+1/2
    T * getVy() { return Vy; }    ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }    ///< Get Velocity component at time t+1/2
    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity
    void forwardstepStress(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity
    void backwardstepVelocity(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity
    void backwardstepStress(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity

   // Insert source functions
    void insertSource(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling 
    void recordData(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it); ///< Record data from modeling 


private:
    T *Sxx;  // Stress component at time t+1
    T *Syy;  // Stress component at time t+1
    T *Szz;  // Stress component at time t+1
    T *Sxz;  // Stress component at time t+1
    T *Syz;  // Stress component at time t+1
    T *Sxy;  // Stress component at time t+1

    T *Vx; // Velocity component at time t+1/2
    T *Vy; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic3D<T>> Pml; // Associated Pml class
    bool adjoint; 
};

/** The 2D Elastic displacement-stress WAVES class
 *
*/

template<typename T>
class WavesElastic2D_DS: public Waves<T> {
public:
    WavesElastic2D_DS();	///< Constructor
    ~WavesElastic2D_DS();	///< Destructor
    WavesElastic2D_DS(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesElastic2D_DS(std::shared_ptr<rockseis::ModelElastic2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor

    // Get functions
    std::shared_ptr<PmlElastic2D<T>> getPml() { return Pml; } ///< Get Pml 
    T * getSxx() { return Sxx; }  ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }  ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }  ///< Get Stress component at time t+1

    T * getUx1() { return Ux1; }  ///< Get Displacement component at time t
    T * getUz1() { return Uz1; }  ///< Get Displacement component at time t

    T * getUx2() { return Ux2; }  ///< Get Displacement component at time t+1
    T * getUz2() { return Uz2; }  ///< Get Displacement component at time t+1

    // Time stepping functions
    void forwardstepDisplacement(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with Stress

    // Insert source functions
    void insertPressuresource(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Pressure type)
    void insertForcesource(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Force type )

    // Record data at receivers functions
    void recordData(std::shared_ptr<ModelElastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Displacement type or Pressure )

    void roll();  // Roll the displacement pointers

private:
    T *Sxx;  // Stress component at time t
    T *Szz; // Stress component at time t
    T *Sxz; // Stress component at time t
    T *Ux1; // Displacement component at time t
    T *Ux2; // Displacement component at time t+1
    T *Uz1; // Displacement component at time t
    T *Uz2; // Displacement component at time t+1
    std::shared_ptr<PmlElastic2D<T>> Pml; // Associated Pml class
};

/** The 3D Elastic Displacement-stress WAVES class
 *
 */
template<typename T>
class WavesElastic3D_DS: public Waves<T> {
public:
    WavesElastic3D_DS();	///< Constructor
    WavesElastic3D_DS(const int _nx, const int _ny, const int _nz, const int _nt, const int _L, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    WavesElastic3D_DS(std::shared_ptr<rockseis::ModelElastic3D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesElastic3D_DS();	///< Destructor

    // Get functions
    std::shared_ptr<PmlElastic3D<T>> getPml() { return Pml; }
    T * getSxx() { return Sxx; }    ///< Get Stress component at time t
    T * getSyy() { return Syy; }    ///< Get Stress component at time t
    T * getSzz() { return Szz; }    ///< Get Stress component at time t
    T * getSyz() { return Syz; }    ///< Get Stress component at time t
    T * getSxz() { return Sxz; }    ///< Get Stress component at time t
    T * getSxy() { return Sxy; }    ///< Get Stress component at time t
    T * getUx1() { return Ux1; }    ///< Get Displacement component at time t
    T * getUx2() { return Ux2; }    ///< Get Displacement component at time t+1
    T * getUy1() { return Uy1; }    ///< Get Displacement component at time t
    T * getUy2() { return Uy2; }    ///< Get Displacement component at time t+1
    T * getUz1() { return Uz1; }    ///< Get Displacement component at time t
    T * getUz2() { return Uz2; }    ///< Get Displacement component at time t+1

    // Time stepping functions
    void forwardstepDisplacement(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity 

    // Insert source functions
    void insertForcesource(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Displacement type or Pressure )
    void insertPressuresource(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Displacement type or Pressure )
    void recordData(std::shared_ptr<ModelElastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Displacement type or Pressure )

    void roll();  // Roll the displacement pointers
    
private:
    T *Sxx;  // Stress component at time t
    T *Syy;  // Stress component at time t
    T *Szz;  // Stress component at time t
    T *Sxz;  // Stress component at time t
    T *Syz;  // Stress component at time t
    T *Sxy;  // Stress component at time t
    
    T *Ux1; // Displacement component at time t
    T *Ux2; // Displacement component at time t+1
    T *Uy1; // Displacement component at time t
    T *Uy2; // Displacement component at time t+1
    T *Uz1; // Displacement component at time t
    T *Uz2; // Displacement component at time t+1
    std::shared_ptr<PmlElastic3D<T>> Pml; // Associated Pml class
};

/** The 2D Viscoelastic WAVES class
 *
 */
template<typename T>
class WavesViscoelastic2D: public Waves<T> {
public:
    WavesViscoelastic2D();	///< Constructor
    ~WavesViscoelastic2D();	///< Destructor
    WavesViscoelastic2D(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesViscoelastic2D(std::shared_ptr<rockseis::ModelViscoelastic2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor

    // Get functions
    std::shared_ptr<PmlElastic2D<T>> getPml() { return Pml; } ///< Get Pml 
    T * getSxx() { return Sxx; }  ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }  ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }  ///< Get Stress component at time t+1

    T * getMxx() { return Mxx; }  ///< Get Stress component at time t+1
    T * getMzz() { return Mzz; }  ///< Get Stress component at time t+1
    T * getMxz() { return Mxz; }  ///< Get Stress component at time t+1

    T * getVx() { return Vx; }  ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }  ///< Get Velocity component at time t+1/2
    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with Stress

    void backwardstepVelocity(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step backward with particle velocity 
    void backwardstepStress(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with Stress

    // Insert source functions
    void insertSource(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Velocity type or Pressure )

    // Record data at receivers functions
    void recordData(std::shared_ptr<ModelViscoelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Velocity type or Pressure )

private:
    T *Sxx;  // Stress component at time t+1
    T *Szz; // Stress component at time t+1
    T *Sxz; // Stress component at time t+1
    T *Mxx;  // Memory component at time t+1
    T *Mzz; // Memory component at time t+1
    T *Mxz; // Memory component at time t+1
    T *Vx; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic2D<T>> Pml; // Associated Pml class
    bool adjoint; 
};

/** The 3D Viscoelastic WAVES class
 *
 */
template<typename T>
class WavesViscoelastic3D: public Waves<T> {
public:
    WavesViscoelastic3D();	///< Constructor
    WavesViscoelastic3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _L, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    WavesViscoelastic3D(std::shared_ptr<rockseis::ModelViscoelastic3D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesViscoelastic3D();	///< Destructor

    // Get functions
    std::shared_ptr<PmlElastic3D<T>> getPml() { return Pml; }
    T * getSxx() { return Sxx; }    ///< Get Stress component at time t+1
    T * getSyy() { return Syy; }    ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }    ///< Get Stress component at time t+1
    T * getSyz() { return Syz; }    ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }    ///< Get Stress component at time t+1
    T * getSxy() { return Sxy; }    ///< Get Stress component at time t+1
    T * getMxx() { return Mxx; }    ///< Get Stress component at time t+1
    T * getMyy() { return Myy; }    ///< Get Stress component at time t+1
    T * getMzz() { return Mzz; }    ///< Get Stress component at time t+1
    T * getMyz() { return Myz; }    ///< Get Stress component at time t+1
    T * getMxz() { return Mxz; }    ///< Get Stress component at time t+1
    T * getMxy() { return Mxy; }    ///< Get Stress component at time t+1
    T * getVx() { return Vx; }    ///< Get Velocity component at time t+1/2
    T * getVy() { return Vy; }    ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }    ///< Get Velocity component at time t+1/2

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity
    void forwardstepStress(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity

    void backwardstepVelocity(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity
    void backwardstepStress(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity

    // Insert source functions
    void insertSource(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Velocity type or Pressure )
    void recordData(std::shared_ptr<ModelViscoelastic3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Velocity type or Pressure )


private:
    T *Sxx;  // Stress component at time t+1
    T *Syy;  // Stress component at time t+1
    T *Szz;  // Stress component at time t+1
    T *Sxz;  // Stress component at time t+1
    T *Syz;  // Stress component at time t+1
    T *Sxy;  // Stress component at time t+1

    T *Mxx;  // Memory component at time t+1
    T *Myy;  // Memory component at time t+1
    T *Mzz;  // Memory component at time t+1
    T *Mxz;  // Memory component at time t+1
    T *Myz;  // Memory component at time t+1
    T *Mxy;  // Memory component at time t+1

    T *Vx; // Velocity component at time t+1/2
    T *Vy; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic3D<T>> Pml; // Associated Pml class
};

/** The 2D Vti WAVES class
 *
 */
template<typename T>
class WavesVti2D: public Waves<T> {
public:
    WavesVti2D();	///< Constructor
    ~WavesVti2D();	///< Destructor
    WavesVti2D(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesVti2D(std::shared_ptr<rockseis::ModelVti2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor

    // Get functions
    std::shared_ptr<PmlElastic2D<T>> getPml() { return Pml; } ///< Get Pml 
    T * getSxx() { return Sxx; }  ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }  ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }  ///< Get Stress component at time t+1

    T * getVx() { return Vx; }  ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }  ///< Get Velocity component at time t+1/2

    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with Stress
    void backwardstepVelocity(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step backward with particle velocity 
    void backwardstepStress(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with Stress

    // Insert source functions
    void insertSource(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling 

    // Record data at receivers functions
    void recordData(std::shared_ptr<ModelVti2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling 



private:
    T *Sxx;  // Stress component at time t+1
    T *Szz; // Stress component at time t+1
    T *Sxz; // Stress component at time t+1
    T *Vx; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic2D<T>> Pml; // Associated Pml class
    bool adjoint; 

};

/** The 3D Ortho WAVES class
 *
 */
template<typename T>
class WavesOrtho3D: public Waves<T> {
public:
    WavesOrtho3D();	///< Constructor
    WavesOrtho3D(const int _nx, const int _ny, const int _nz, const int _nt, const int _L, const T _dx, const T _dy, const T _dz, const T _dt, const T _ox, const T _oy, const T _oz, const T _ot);	///< Constructor
    WavesOrtho3D(std::shared_ptr<rockseis::ModelOrtho3D<T>> model, int _nt, T _dt, T _ot);	///< Constructor
    ~WavesOrtho3D();	///< Destructor

    // Get functions
    std::shared_ptr<PmlElastic3D<T>> getPml() { return Pml; }
    T * getSxx() { return Sxx; }    ///< Get Stress component at time t+1
    T * getSyy() { return Syy; }    ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }    ///< Get Stress component at time t+1
    T * getSyz() { return Syz; }    ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }    ///< Get Stress component at time t+1
    T * getSxy() { return Sxy; }    ///< Get Stress component at time t+1
    T * getVx() { return Vx; }    ///< Get Velocity component at time t+1/2
    T * getVy() { return Vy; }    ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }    ///< Get Velocity component at time t+1/2
    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity
    void forwardstepStress(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with particle velocity
    void backwardstepVelocity(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity
    void backwardstepStress(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step backward with particle velocity

   // Insert source functions
    void insertSource(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int it); ///< Insert source for modeling 
    void recordData(std::shared_ptr<ModelOrtho3D<T>> model, std::shared_ptr<rockseis::Data3D<T>> data, bool maptype, int it); ///< Record data from modeling 


private:
    T *Sxx;  // Stress component at time t+1
    T *Syy;  // Stress component at time t+1
    T *Szz;  // Stress component at time t+1
    T *Sxz;  // Stress component at time t+1
    T *Syz;  // Stress component at time t+1
    T *Sxy;  // Stress component at time t+1

    T *Vx; // Velocity component at time t+1/2
    T *Vy; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlElastic3D<T>> Pml; // Associated Pml class
    bool adjoint; 
};

/** The 2D Poroelastic WAVES class
 *
 */
template<typename T>
class WavesPoroelastic2D: public Waves<T> {
public:
    WavesPoroelastic2D();	///< Constructor
    ~WavesPoroelastic2D();	///< Destructor
    WavesPoroelastic2D(const int _nx, const int _nz, const int _nt, const int _L, const T _dx, const T _dz, const T _dt, const T _ox, const T _oz, const T _ot);	///< Constructor
    WavesPoroelastic2D(std::shared_ptr<rockseis::ModelPoroelastic2D<T>> model, int _nt, T _dt, T _ot);	///< Constructor

    // Get functions
    std::shared_ptr<PmlPoroelastic2D<T>> getPml() { return Pml; } ///< Get Pml 
    T * getSxx() { return Sxx; }  ///< Get Stress component at time t+1
    T * getSzz() { return Szz; }  ///< Get Stress component at time t+1
    T * getSxz() { return Sxz; }  ///< Get Stress component at time t+1

    T * getP() { return P; }  ///< Get Pressure

    T * getVx() { return Vx; }  ///< Get Velocity component at time t+1/2
    T * getVz() { return Vz; }  ///< Get Velocity component at time t+1/2
    T * getQx() { return Qx; }  ///< Get Velocity component at time t+1/2
    T * getQz() { return Qz; }  ///< Get Velocity component at time t+1/2
    bool getAdjoint() {return adjoint;} ///< Return adjoint flag
    void setAdjoint(); ///< Set adjoint flag and allocate work array for adjoint computation

    // Time stepping functions
    void forwardstepVelocity(std::shared_ptr<ModelPoroelastic2D<T>> model, std::shared_ptr<Der<T>> der); ///< Advance one time step forward with particle velocity 
    void forwardstepStress(std::shared_ptr<ModelPoroelastic2D<T>> model, std::shared_ptr<Der<T>> der);  ///< Advance one time step forward with Stress

    // Insert source functions
    void insertSource(std::shared_ptr<ModelPoroelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int it); ///< Insert source for modeling ( Source types can be of Velocity type or Pressure )

    // Record data at receivers functions
    void recordData(std::shared_ptr<ModelPoroelastic2D<T>> model, std::shared_ptr<rockseis::Data2D<T>> data, bool maptype, int it); ///< Record data from modeling ( Data types can be of Velocity type or Pressure )

private:
    T *Sxx;  // Stress component at time t+1
    T *Szz; // Stress component at time t+1
    T *Sxz; // Stress component at time t+1
    T *P;  // Pressure component
    T *Vx; // Velocity component at time t+1/2
    T *Vz; // Velocity component at time t+1/2
    T *Qx; // Velocity component at time t+1/2
    T *Qz; // Velocity component at time t+1/2
    T *wrk; // work array used in adjoint state modelling
    std::shared_ptr<PmlPoroelastic2D<T>> Pml; // Associated Pml class
    bool adjoint; 
};


}
#endif //WAVES_H
