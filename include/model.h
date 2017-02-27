#ifndef MODEL_H
#define MODEL_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "file.h"

namespace rockseis {

// =============== ABSTRACT MODEL CLASS =============== //
/** The abstract model class
 *
 */
template<typename T>
class Model {
public:
    Model();		///< Constructor
    Model(const int _dim);	///< Constructor with dimension
    Model(const int _dim, const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs); ///< Constructer with all variables
    virtual ~Model();	///< Destructor
    
    // Get functions
    int getDim() { return dim; }		///< Get dimension
    int getNx() { return geometry->getN(1); }		///< Get Nx
    int getNy() { return geometry->getN(2); }		///< Get Ny
    int getNz() { return geometry->getN(3); }		///< Get Nz
    int getLpml() { return lpml; }		///< Get Lpml
    int getNx_pml() { return geometry->getN(1) + 2 * lpml; }	///< Nx_pml = Nx + 2*lpml 
    int getNy_pml() { return geometry->getN(2) + 2 * lpml; }	///< Ny_pml = Ny + 2*lpml 
    int getNz_pml() { return geometry->getN(3) + 2 * lpml; }	///< Nz_pml = Nz + 2*lpml 
    bool getFs() { return fs; }		///< Get fs
    T getDx() { return geometry->getD(1); }		///< Get Dx
    T getDy() { return geometry->getD(2); }		///< Get Dy
    T getDz() { return geometry->getD(3); }		///< Get Dz
    T getOx() { return geometry->getO(1); }		///< Get Ox
    T getOy() { return geometry->getO(2); }		///< Get Oy
    T getOz() { return geometry->getO(3); }		///< Get Oz
    std::shared_ptr<Geometry<T>> getGeom() { return geometry; } ///< Get geometry
    
    // Set functions
    void setNx(const int _nx) { geometry->setN(1, _nx); }	///< Set Nx
    void setNy(const int _ny) { geometry->setN(2, _ny); }	///< Set Ny
    void setNz(const int _nz) { geometry->setN(3, _nz); }	///< Set Nz
    void setLpml(const int _lpml) { lpml = _lpml; }	///< Set size of padding boundary for PML
    void setNx_pml() { nx_pml = geometry->getN(1) + 2 * lpml; }	///< Nx_pml = Nx + 2*lpml 
    void setNy_pml() { ny_pml = geometry->getN(2) + 2 * lpml; }	///< Ny_pml = Ny + 2*lpml 
    void setNz_pml() { nz_pml = geometry->getN(3) + 2 * lpml; }	///< Nz_pml = Nz + 2*lpml 
    void setFs(const bool _fs) { fs = _fs; }	///< Set free surface
    void setDx(const T _dx) { geometry->setD(1, _dx); }	///< Set Dx
    void setDy(const T _dy) { geometry->setD(2, _dy); }	///< Set Dy
    void setDz(const T _dz) { geometry->setD(3, _dz); }	///< Set Dz
    void setOx(const T _ox) { geometry->setO(1, _ox); }	///< Set Ox
    void setOy(const T _oy) { geometry->setO(2, _oy); }	///< Set Oy
    void setOz(const T _oz) { geometry->setO(3, _oz); }	///< Set Oz
    void setDim(const int _dim) { dim = _dim; } 	///< Set the dimension
    
    // I/O functions
    virtual void readModel() = 0;	///< Read model (virtual)
    virtual void writeModel() = 0;	///< write model (virtual)
    // PADDING AND STAGGERING FUNCTIONS 
    /** Pads 2-D model.
     * Pads the model by copying the edges over the padded area.  
     * The size of the padded arrays must be (nx+2*pad) * (ny+2*pad).
     * */
    void padmodel2d(T *padded, T *model, const int nx, const int ny, const int pad);

    /** Pads 3-D model.
     * Pads the model by copying the edges over the padded volume.  
     * The size of the padded arrays must be (nx+2*pad) * (ny+2*pad) * (nz+2*pad).
     * */
    void padmodel3d(T *padded, T *model, const int nx, const int ny, const int nz, const int pad); 

    /** Stagger the model forward in x direction.
     *  Uses an arithmetic running average
     * */
    void staggermodel_x(T *model, const int nx, const int ny, const int nz);
    /** Stagger the model forward in y direction.
     *  Uses an arithmetic running average
     * */
    void staggermodel_y(T *model, const int nx, const int ny, const int nz); 
    /** Stagger the model forward in z direction.
     *  Uses an arithmetic running average
     * */
    void staggermodel_z(T *model, const int nx, const int ny, const int nz); 

private:
    int dim; 
    int lpml;
    std::shared_ptr<Geometry<T>> geometry;
    int nx_pml;
    int ny_pml;
    int nz_pml;
    bool fs;
};

// =============== 2D ACOUSTIC MODEL CLASS =============== //
/** The 2D acoustic model class
 *
 */
template<typename T>
class ModelAcoustic2D: public Model<T> {
public:
    ModelAcoustic2D();	///< Constructor
    ModelAcoustic2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs);	///< Constructor
    ModelAcoustic2D(std::string _Vpfile, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelAcoustic2D();	///< Destructor
    
    // I/O functions
    void readModel();	 ///< Read a model from file
    void writeModel(); ///< Write a model to file
    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getR() { return R; }		///< Get R
    T *getRx() { return Rx; }		///< Get Rx
    T *getRz() { return Rz; }		///< Get Rz
    T *getL() { return L; }		///< Get L
    std::string getVpfile() { return Vpfile; }
    std::string getRfile() { return Rfile; }
    
    /** Stagger model functions. 
    It creates the padded Rx, Rz and L from the non-padded models R and Vp. 
    */
    void staggerModels(); 
    
private:
    T *Vp;  ///< P-wave velocity
    T *R;   ///< Density
    T *L;   ///< Bulk modulus
    T *Rx;  ///< Staggered inverse of density in x
    T *Rz;  ///< Staggered inverse of density in z
    std::string Vpfile; ///< Filename to vp model
    std::string Rfile; ///< Filename to density model

};

// =============== 3D ACOUSTIC MODEL CLASS =============== //
/** The 3D acoustic model class
 *
 */
template<typename T>
class ModelAcoustic3D: public Model<T> {
public:
    ModelAcoustic3D();	///< Default constructor (Not to be used)
    ModelAcoustic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs); ///< Constructor
    ModelAcoustic3D(std::string _Vpfile, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelAcoustic3D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from file
    void writeModel();	///< NOT IMPLEMENTED
    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getR() { return R; }		///< Get R
    T *getRx() { return Rx; }		///< Get Rx
    T *getRy() { return Ry; }		///< Get Ry
    T *getRz() { return Rz; }		///< Get Rz
    T *getL() { return L; }		///< Get L
    std::string getVpfile() { return Vpfile; }
    std::string getRfile() { return Rfile; }

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz and L from the non-padded models R and Vp. 
    */
    void staggerModels();

private:
    T *Vp;  ///< P-wave velocity
    T *R;   ///< Density
    T *L;   ///< Bulk modulus
    T *Rx;  ///< Staggered inverse of density in x
    T *Ry;  ///< Staggered inverse of density in y
    T *Rz;  ///< Staggered inverse of density in z
    std::string Vpfile; ///< Filename to vp model
    std::string Rfile; ///< Filename to density model
};

// =============== 2D ELASTIC MODEL CLASS =============== //
/** The 2D elastic model class
 *
 */
template<typename T>
class ModelElastic2D: public Model<T> {
public:
    ModelElastic2D();	///< Constructor
    ModelElastic2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs);	///< Constructor
    ModelElastic2D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelElastic2D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from file
    void writeModel();	///< NOT IMPLEMENTED

    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getVs() { return Vs; }	///< Get Vs
    T *getR() { return R; }		///< Get R
    T *getL() { return L; }		///< Get L
    T *getL2M() { return L2M; }		///< Get L
    T *getM() { return M; }		///< Get M
    T *getRx() { return Rx; }		///< Get Rx
    T *getRz() { return Rz; }		///< Get Rz
    std::string getVpfile() { return Vpfile; }
    std::string getVsfile() { return Vsfile; }
    std::string getRfile() { return Rfile; }

    /** Stagger model functions. 
    It creates the padded Rx, Rz, L, L2M and M from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();

private:
    T *Vp;  // P-wave velocity
    T *Vs;  // S-wave velocity
    T *R;   // Density
    T *L;   // Lame lambda  (padded)
    T *L2M; // Lame lambda + 2 Mu  (padded)
    T *M;   // Lame Mu  (padded)
    T *Rx;  // Staggered inverse of density in x (padded)
    T *Rz;  // Staggered inverse of density in z (padded)
    std::string Vpfile; ///< Filename to vp model
    std::string Vsfile; ///< Filename to vs model
    std::string Rfile; ///< Filename to density model
};

// =============== 3D ELASTIC MODEL CLASS =============== //
/** The 3D elastic model class
 *
 */
template<typename T>
class ModelElastic3D: public Model<T> {
public:
    ModelElastic3D();	///< Constructor
    ModelElastic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs); ///< Constructor
    ModelElastic3D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelElastic3D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from files
    void writeModel();	///< Write a to files

    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getVs() { return Vs; }	///< Get Vs
    T *getR() { return R; }		///< Get R
    T *getL() { return L; }		///< Get L
    T *getL2M() { return L2M; }		///< Get L2M
    T *getM_xz() { return M_xz; }	///< Get M_xz
    T *getM_yz() { return M_yz; }	///< Get M_yz
    T *getM_xy() { return M_xy; }	///< Get M_xy
    T *getRx() { return Rx; }		///< Get Rx
    T *getRy() { return Ry; }		///< Get Ry
    T *getRz() { return Rz; }		///< Get Rz
    std::string getVpfile() { return Vpfile; }
    std::string getVsfile() { return Vsfile; }
    std::string getRfile() { return Rfile; }

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz, L, L2M, M_xz, M_yz, M_xy from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();

private:
    T *Vp;  ///< P-wave velocity
    T *Vs;  ///< S-wave velocity
    T *R;   ///< Density 
    T *L;   ///< Lame lambda  (padded)
    T *L2M;   ///< Lame lambda  (padded)
    T *M_xz;   ///< Lame Mu  (padded)
    T *M_yz;   ///< Lame Mu  (padded)
    T *M_xy;   ///< Lame Mu  (padded)
    T *Rx;  ///< Staggered inverse of density in x (padded)
    T *Ry;  ///< Staggered inverse of density in y (padded)
    T *Rz;  ///< Staggered inverse of density in z (padded)
    std::string Vpfile; ///< Filename to vp model
    std::string Vsfile; ///< Filename to vs model
    std::string Rfile; ///< Filename to density model
};


}
#endif //MODEL_H
