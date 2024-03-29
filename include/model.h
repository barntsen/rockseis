#ifndef MODEL_H
#define MODEL_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "utils.h"
#include "geometry.h"
#include "domain.h"
#include "file.h"
#include "data.h"

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
    ~Model();	///< Destructor
    
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
    std::shared_ptr<Domain<T>> getDomain() { return domain; } ///< Get domain
    bool getRealized() { return realized; } ///< Check if model is allocated
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
    void setRealized(const bool val) { realized = val; } ///< Set if model is allocated


    // Calculate local model sizes
    void  getLocalsize2d(std::shared_ptr<Data2D<T>> data, T aperture, bool map, off_t *start, size_t *size);
    void  getLocalsize3d(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, off_t *start_x, size_t *size_x, off_t *start_y, size_t *size_y);

    // PADDING AND STAGGERING FUNCTIONS 
    /** Pads 1-D model.
     * Pads the model by copying the edges over the padded area.  
     * The size of the padded arrays must be (nx+2*pad).
     * */
    void padmodel1d(T *padded, T *model, const int nx, const int pad);

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

protected:
    T getMin(T* Model); ///< Finds minimum of a model (to be used with model arrays only)
    T getMax(T* Model); ///< Finds maximum of a model (to be used with model arrays only)

private:
    int dim; 
    int lpml;
    std::shared_ptr<Geometry<T>> geometry;
    std::shared_ptr<Domain<T>> domain;
    int nx_pml;
    int ny_pml;
    int nz_pml;
    bool fs;
    bool realized;  // Set to 1 when the model is allocated to the correct size
};

// =============== 2D EIKONAL MODEL CLASS =============== //
/** The 2D Eikonal model class
 *
 */
template<typename T>
class ModelEikonal2D: public Model<T> {
public:
    ModelEikonal2D();	///< Constructor
    ModelEikonal2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz);	///< Constructor
    ModelEikonal2D(std::string _Velocityfile, const int lpml);	///< Constructor
    ~ModelEikonal2D();	///< Destructor
    
    // I/O functions
    void readVelocity();	 ///< Read a model from file
    void writeVelocity(); ///< Write the model to file
    // Get functions
    T *getVelocity() { return Velocity; }	///< Get Velocity
    T *getL() { return L; }		///< Get L
    std::string getVelocityfile() { return Velocityfile; }
    void setVelocityfile(std::string name) { Velocityfile = name; }
    std::shared_ptr<ModelEikonal2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    T getMinvel() {return this->getMin(Velocity); } ///< Returns min of model
    T getMaxvel() {return this->getMax(Velocity); } ///< Returns max of model

    /** Stagger model functions. 
    It creates the padded model from the non-padded model. 
    */
    void Expand(); 

    /** Create model
    It creates an empty model 
    */
    void createModel();
    
private:
    T *Velocity;  ///< Velocity
    T *L;   ///< Expanded model
    std::string Velocityfile; ///< Filename to model
};

// =============== 3D EIKONAL MODEL CLASS =============== //
/** The 3D Eikonal model class
 *
 */
template<typename T>
class ModelEikonal3D: public Model<T> {
public:
    ModelEikonal3D();	///< Default constructor (Not to be used)
    ModelEikonal3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz); ///< Constructor
    ModelEikonal3D(std::string _Velocityfile, const int lpml);	///< Constructor
    ~ModelEikonal3D();	///< Destructor
    
    // I/O functions
    void readVelocity();	 ///< Read a model from file
    void writeVelocity(); ///< Write the model to file
    // Get functions
    T *getVelocity() { return Velocity; }	///< Get Velocity
    T *getL() { return L; }		///< Get L
    std::string getVelocityfile() { return Velocityfile; }
    void setVelocityfile(std::string name) { Velocityfile = name; }
    std::shared_ptr<ModelEikonal3D<T>> getLocal(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map);
    T getMinvel() {return this->getMin(Velocity); } ///< Returns min of model
    T getMaxvel() {return this->getMax(Velocity); } ///< Returns max of model

    /** Stagger model functions. 
    It creates the padded model from the non-padded model. 
    */
    void Expand(); 

    /** Create model
    It creates an empty model 
    */
    void createModel();
    
private:
    T *Velocity;  ///< Velocity
    T *L;   ///< Expanded model
    std::string Velocityfile; ///< Filename to model
};



// =============== 1D ACOUSTIC MODEL CLASS =============== //
/** The 1D acoustic model class
 *
 */
template<typename T>
class ModelAcoustic1D: public Model<T> {
public:
    ModelAcoustic1D();	///< Constructor
    ModelAcoustic1D(const int _nz, const int lpml, const T _dz, const T _oz, const bool _fs);	///< Constructor
    ModelAcoustic1D(std::string _Vpfile, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelAcoustic1D();	///< Destructor
    
    // I/O functions
    void readModel();	 ///< Read a model from file
    void writeModel(); ///< Write a model to file
    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getR() { return R; }		///< Get R
    T *getRz() { return Rz; }		///< Get Rz
    T *getL() { return L; }		///< Get L
    std::string getVpfile() { return Vpfile; }
    std::string getRfile() { return Rfile; }
    void setVpfile(std::string name) { Vpfile = name; }
    void setRfile(std::string name) { Rfile = name; }

    /** Stagger model functions. 
    It creates the padded Rx, Rz and L from the non-padded models R and Vp. 
    */
    void staggerModels(); 

    /** Create model
    It creates an empty model of Vp and R
    */
    void createModel();
    
private:
    T *Vp;  ///< P-wave velocity
    T *R;   ///< Density
    T *L;   ///< Bulk modulus
    T *Rz;  ///< Staggered inverse of density in z
    std::string Vpfile; ///< Filename to vp model
    std::string Rfile; ///< Filename to density model

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
    void writeVp(); ///< Write only the Vp model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeR(); } ///< Write a model to file
    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getR() { return R; }		///< Get R
    T *getRx() { return Rx; }		///< Get Rx
    T *getRz() { return Rz; }		///< Get Rz
    T *getL() { return L; }		///< Get L
    std::string getVpfile() { return Vpfile; }
    std::string getRfile() { return Rfile; }
    void setVpfile(std::string name) { Vpfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    std::shared_ptr<ModelAcoustic2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    std::shared_ptr<ModelAcoustic2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd, const int order); ///< Returns a model of a domain
    std::shared_ptr<ModelAcoustic2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd0, const int nd1, const int order); ///< Returns a model of a domain
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxR() {return this->getMax(R); } ///< Returns max R

    /** Stagger model functions. 
    It creates the padded Rx, Rz and L from the non-padded models R and Vp. 
    */
    void staggerModels(); 

    /** Create model
    It creates an empty model of Vp and R
    */
    void createModel();
    void createPaddedmodel();
    
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
    void writeVp(); ///< Write only the Vp model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeR(); } ///< Write a model to file
    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getR() { return R; }		///< Get R
    T *getRx() { return Rx; }		///< Get Rx
    T *getRy() { return Ry; }		///< Get Ry
    T *getRz() { return Rz; }		///< Get Rz
    T *getL() { return L; }		///< Get L
    std::string getVpfile() { return Vpfile; }
    std::string getRfile() { return Rfile; }
    void setVpfile(std::string name) { Vpfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxR() {return this->getMax(R); } ///< Returns max R

    std::shared_ptr<ModelAcoustic3D<T>> getLocal(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map);
    std::shared_ptr<ModelAcoustic3D<T>> getDomainmodel(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map, const int d, const int nd, const int order); ///< Returns a model of a domain
    std::shared_ptr<ModelAcoustic3D<T>> getDomainmodel(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order); ///< Returns a model of a domain

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz and L from the non-padded models R and Vp. 
    */
    void staggerModels();

    /** Create model
    It creates an empty model of Vp and R
    */
    void createModel();
    void createPaddedmodel();


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

    void writeVp(); ///< Write only the Vp model to file
    void writeVs(); ///< Write only the Vs model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeVs(); writeR(); } ///< Write a model to file

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
    void setVpfile(std::string name) { Vpfile = name; }
    void setVsfile(std::string name) { Vsfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinVs() {return this->getMin(Vs); } ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxVs() {return this->getMax(Vs); } ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R
    /** Stagger model functions. 
    It creates the padded Rx, Rz, L, L2M and M from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();
    std::shared_ptr<ModelElastic2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    std::shared_ptr<ModelElastic2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd0, const int nd1, const int order); ///< Returns a model of a domain

    /** Create model
    It creates an empty model of Vp, Vs and R
    */
    void createModel();
    void createPaddedmodel();

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
    void writeVp(); ///< Write only the Vp model to file
    void writeVs(); ///< Write only the Vs model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeVs(); writeR(); } ///< Write a model to file

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
    void setVpfile(std::string name) { Vpfile = name; }
    void setVsfile(std::string name) { Vsfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinVs() {return this->getMin(Vs); } ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxVs() {return this->getMax(Vs); } ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R

    std::shared_ptr<ModelElastic3D<T>> getLocal(std::shared_ptr<Data3D<T>>, T aperture_x, T aperture_y, bool map);
    std::shared_ptr<ModelElastic3D<T>> getDomainmodel(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order); ///< Returns a model of a domain

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz, L, L2M, M_xz, M_yz, M_xy from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();

    /** Create model
    It creates an empty model of Vp, Vs and R
    */
    void createModel();
    void createPaddedmodel();

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

// =============== 2D VISCO ELASTIC MODEL CLASS =============== //
/** The 2D visco-elastic model class
 *
 */
template<typename T>
class ModelViscoelastic2D: public Model<T> {
public:
    ModelViscoelastic2D();	///< Constructor
    ModelViscoelastic2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz, const T f0, const bool _fs);	///< Constructor
    ModelViscoelastic2D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, std::string _Qpfile, std::string _Qsfile, const int lpml, const T f0, const bool _fs);	///< Constructor
    ~ModelViscoelastic2D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from file

    void writeVp(); ///< Write only the Vp model to file
    void writeVs(); ///< Write only the Vs model to file
    void writeQp(); ///< Write only the Qp model to file
    void writeQs(); ///< Write only the Qs model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeVs(); writeQp(); writeQs(); writeR(); } ///< Write a model to file

    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getQp() { return Qp; }	///< Get Qp
    T *getVs() { return Vs; }	///< Get Vs
    T *getQs() { return Qs; }	///< Get Qs
    T *getR() { return R; }		///< Get R
    T *getM() { return M; }		///< Get M
    T *getL2M() { return L2M; }		///< Get L2M
    T *getM_xz() { return M_xz; }	///< Get M_xz
    T *getRx() { return Rx; }		///< Get Rx
    T *getRz() { return Rz; }		///< Get Rz
    T *getTp() { return Tp; }		///< Get Tp
    T *getTs() { return Ts; }		///< Get Ts
    T *getTs_xz() { return Ts_xz; }		///< Get Ts_xz
    T getF0() { return f0; } ///< Get the dominant frequency
    std::string getVpfile() { return Vpfile; }
    std::string getVsfile() { return Vsfile; }
    std::string getRfile() { return Rfile; }
    std::string getQpfile() { return Qpfile; }
    std::string getQsfile() { return Qsfile; }
    void setVpfile(std::string name) { Vpfile = name; }
    void setVsfile(std::string name) { Vsfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    void setQpfile(std::string name) { Qpfile = name; }
    void setQsfile(std::string name) { Qsfile = name; }
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinVs() {return this->getMin(Vs); } ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxVs() {return this->getMax(Vs); } ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R
    /** Stagger model functions. 
    It creates the padded Rx, Rz, M_xz, L2M and M from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();
    std::shared_ptr<ModelViscoelastic2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    std::shared_ptr<ModelViscoelastic2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd0, const int nd1, const int order); ///< Returns a model of a domain

    /** Create model
    It creates an empty model of Vp, Vs and R
    */
    void createModel();
    void createPaddedmodel();

private:
    T *Vp;  // P-wave velocity
    T *Vs;  // S-wave velocity
    T *Qp;  // Quality factor for P-wave
    T *Qs;  // Quality factor for S-wave
    T *R;   // Density
    T *L2M; // Lame lambda + 2 Mu  (padded)
    T *M;   // Lame Mu  (padded)
    T *M_xz; // Lame Mu  (staggered)
    T *Rx;  // Staggered inverse of density in x (padded)
    T *Rz;  // Staggered inverse of density in z (padded)
    T *Tp;  // Relaxation time for P-wave
    T *Ts;  // Relaxation time for S-wave
    T *Ts_xz;  // Relaxation time for S-wave staggered
    T f0; // Center frequency
    std::string Vpfile; ///< Filename to vp model
    std::string Vsfile; ///< Filename to vs model
    std::string Rfile; ///< Filename to density model
    std::string Qpfile; ///< Filename to P-wave Q model
    std::string Qsfile; ///< Filename to S-wave Q model
};

// =============== 3D VISCOELASTIC MODEL CLASS =============== //
/** The 3D elastic model class
 *
 */
template<typename T>
class ModelViscoelastic3D: public Model<T> {
public:
    ModelViscoelastic3D();	///< Constructor
    ModelViscoelastic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const T f0, const bool _fs); ///< Constructor
    ModelViscoelastic3D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, std::string _Qpfile, std::string _Qsfile, const T f0, const int lpml, const bool _fs);	///< Constructor
    ~ModelViscoelastic3D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from files
    void writeVp(); ///< Write only the Vp model to file
    void writeVs(); ///< Write only the Vs model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeVp(); writeVs(); writeR(); } ///< Write a model to file

    // Get functions
    T *getVp() { return Vp; }	///< Get Vp
    T *getQp() { return Qp; }	///< Get Qp
    T *getVs() { return Vs; }	///< Get Vs
    T *getQs() { return Qs; }	///< Get Qs
    T *getR() { return R; }		///< Get R
    T *getM() { return M; }		///< Get L
    T *getL2M() { return L2M; }		///< Get L2M
    T *getM_xz() { return M_xz; }	///< Get M_xz
    T *getM_yz() { return M_yz; }	///< Get M_yz
    T *getM_xy() { return M_xy; }	///< Get M_xy
    T *getTp() { return Tp; }		///< Get Tp
    T *getTs() { return Ts; }		///< Get Ts
    T *getTs_xz() { return Ts_xz; }	///< Get Ts_xz
    T *getTs_yz() { return Ts_yz; }	///< Get Ts_yz
    T *getTs_xy() { return Ts_xy; }	///< Get Ts_xy
    T *getRx() { return Rx; }		///< Get Rx
    T *getRy() { return Ry; }		///< Get Ry
    T *getRz() { return Rz; }		///< Get Rz
    T getF0() { return f0; } ///< Get the dominant frequency
    std::string getVpfile() { return Vpfile; }
    std::string getVsfile() { return Vsfile; }
    std::string getRfile() { return Rfile; }
    std::string getQpfile() { return Qpfile; }
    std::string getQsfile() { return Qsfile; }
    void setVpfile(std::string name) { Vpfile = name; }
    void setVsfile(std::string name) { Vsfile = name; }
    void setRfile(std::string name) { Rfile = name; }
    void setQpfile(std::string name) { Qpfile = name; }
    void setQsfile(std::string name) { Qsfile = name; }
    T getMinVp() {return this->getMin(Vp); } ///< Returns min Vp
    T getMinVs() {return this->getMin(Vs); } ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp() {return this->getMax(Vp); } ///< Returns max Vp
    T getMaxVs() {return this->getMax(Vs); } ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R

    std::shared_ptr<ModelViscoelastic3D<T>> getLocal(std::shared_ptr<Data3D<T>>, T aperture_x, T aperture_y, bool map);
    std::shared_ptr<ModelViscoelastic3D<T>> getDomainmodel(std::shared_ptr<Data3D<T>>,T aperture_x, T aperture_y, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order); ///< Returns a model of a domain

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz, L, L2M, M_xz, M_yz, M_xy from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();

    /** Create model
    It creates an empty model of Vp, Vs and R
    */
    void createModel();
    void createPaddedmodel();

private:
    T *Vp;  ///< P-wave velocity
    T *Vs;  ///< S-wave velocity
    T *Qp;  // Quality factor for P-wave
    T *Qs;  // Quality factor for S-wave
    T *R;   ///< Density 
    T *M;   // Lame Mu  (padded)
    T *L2M;   ///< Lame lambda  (padded)
    T *M_xz;   ///< Lame Mu  (padded)
    T *M_yz;   ///< Lame Mu  (padded)
    T *M_xy;   ///< Lame Mu  (padded)
    T *Tp;  ///< P-wave relax time
    T *Ts;  ///< S-wave relax time
    T *Ts_xz;   ///< Relax times S-wave  (staggered)
    T *Ts_yz;   ///< Relax times S-wave  (staggered)
    T *Ts_xy;   ///< Relax times S-wave  (staggered)
    T *Rx;  ///< Staggered inverse of density in x (padded)
    T *Ry;  ///< Staggered inverse of density in y (padded)
    T *Rz;  ///< Staggered inverse of density in z (padded)
    T f0; // Center frequency
    std::string Vpfile; ///< Filename to vp model
    std::string Vsfile; ///< Filename to vs model
    std::string Rfile; ///< Filename to density model
    std::string Qpfile; ///< Filename to P-wave Q model
    std::string Qsfile; ///< Filename to S-wave Q model
};

// =============== 2D VTI MODEL CLASS =============== //
/** The 2D vti model class
 *
 */
template<typename T>
class ModelVti2D: public Model<T> {
public:
    ModelVti2D();	///< Constructor
    ModelVti2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs);	///< Constructor
    ModelVti2D(std::string _C11file, std::string _C13file, std::string _C33file, std::string _C55file, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelVti2D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from file

    void writeC11(); ///< Write only the C11 model to file
    void writeC13(); ///< Write only the C13 model to file
    void writeC33(); ///< Write only the C33 model to file
    void writeC55(); ///< Write only the C55 model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeC11(); writeC13(); writeC33(); writeC55();} ///< Write a model to file

    // Get functions
    T *getC11() { return C11; }	///< Get C11
    T *getC13() { return C13; }	///< Get C13
    T *getC33() { return C33; }	///< Get C33
    T *getC55() { return C55; }	///< Get C55
    T *getC11p() { return C11p; }	///< Get C11 padded
    T *getC13p() { return C13p; }	///< Get C13 padded
    T *getC33p() { return C33p; }	///< Get C33 padded
    T *getC55p() { return C55p; }	///< Get C55 padded
    T *getR() { return R; }		///< Get R
    T *getRx() { return Rx; }		///< Get Rx
    T *getRz() { return Rz; }		///< Get Rz
    std::string getC11file() { return C11file; }
    std::string getC13file() { return C13file; }
    std::string getC33file() { return C33file; }
    std::string getC55file() { return C55file; }
    std::string getRfile() { return Rfile; }
    void setC11file(std::string name) { C11file = name; }
    void setC13file(std::string name) { C13file = name; }
    void setC33file(std::string name) { C33file = name; }
    void setC55file(std::string name) { C55file = name; }
    void setRfile(std::string name) { Rfile = name; }
    T getMinVp();  ///< Returns min Vp
    T getMinVs();  ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp();  ///< Returns max Vp
    T getMaxVs(); ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R
    /** Stagger model functions. 
    It creates the padded Rx, Rz, L, L2M and M from the non-padded models R, Vp and Vs. 
    */
    void staggerModels();
    std::shared_ptr<ModelVti2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    std::shared_ptr<ModelVti2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd0, const int nd1, const int order); ///< Returns a model of a domain

    /** Create model
    It creates an empty model of Vp, Vs and R
    */
    void createModel();
    void createPaddedmodel();

private:
    T *C11;  // C11 stiffness
    T *C13;  // C13 stiffness
    T *C33;   // C33 stiffness
    T *C55;  // C55 stiffness
    T *C11p;   // C11 padded
    T *C13p; // C13 padded
    T *C33p; // C33 padded
    T *C55p; // C55 padded
    T *R;  // Density model
    T *Rx;  // Staggered inverse of density in x (padded)
    T *Rz;  // Staggered inverse of density in z (padded)
    std::string C11file; ///< Filename to C11 model
    std::string C13file; ///< Filename to C13 model
    std::string C33file; ///< Filename to C33 model
    std::string C55file; ///< Filename to C55 model
    std::string Rfile; ///< Filename to density model
};

// =============== 3D ORTHOROMBIC MODEL CLASS =============== //
/** The 3D orthorombic model class
 *
 */
template<typename T>
class ModelOrtho3D: public Model<T> {
public:
    ModelOrtho3D();	///< Constructor
    ModelOrtho3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs); ///< Constructor
    ModelOrtho3D(std::string _C11file, std::string _C12file, std::string _C13file, std::string _C22file, std::string _C23file, std::string _C33file, std::string _C44file, std::string _C55file, std::string _C66file, std::string _Rfile, const int lpml, const bool _fs);	///< Constructor
    ~ModelOrtho3D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from files
    void writeC11(); ///< Write only the C11 model to file
    void writeC12(); ///< Write only the C12 model to file
    void writeC13(); ///< Write only the C13 model to file
    void writeC22(); ///< Write only the C22 model to file
    void writeC23(); ///< Write only the C23 model to file
    void writeC33(); ///< Write only the C33 model to file
    void writeC44(); ///< Write only the C44 model to file
    void writeC55(); ///< Write only the C55 model to file
    void writeC66(); ///< Write only the C66 model to file
    void writeR(); ///< Write only the Density model to file
    void writeModel() { writeC11(); writeC12(); writeC13(); writeC22(); writeC23(); writeC33(); writeC44(); writeC55(); writeC66(); writeR(); } ///< Write a model to file

    // Get functions
    T *getC11() { return C11; }	///< Get C11
    T *getC12() { return C12; }	///< Get C12
    T *getC13() { return C13; }	///< Get C13
    T *getC22() { return C22; }	///< Get C22
    T *getC23() { return C23; }	///< Get C23
    T *getC33() { return C33; }	///< Get C33
    T *getC44() { return C44; }	///< Get C44
    T *getC55() { return C55; }	///< Get C55
    T *getC66() { return C66; }	///< Get C66
    T *getR() { return R; }		///< Get R
    T *getC11p() { return C11p; }	///< Get C11 padded
    T *getC12p() { return C12p; }	///< Get C12 padded
    T *getC13p() { return C13p; }	///< Get C13 padded
    T *getC22p() { return C22p; }	///< Get C22 padded
    T *getC23p() { return C23p; }	///< Get C23 padded
    T *getC33p() { return C33p; }	///< Get C33 padded
    T *getC44p() { return C44p; }	///< Get C44 padded
    T *getC55p() { return C55p; }	///< Get C55 padded
    T *getC66p() { return C66p; }	///< Get C66 padded
    T *getRx() { return Rx; }		///< Get Rx
    T *getRy() { return Ry; }		///< Get Ry
    T *getRz() { return Rz; }		///< Get Rz
    std::string getC11file() { return C11file; }
    std::string getC12file() { return C12file; }
    std::string getC13file() { return C13file; }
    std::string getC22file() { return C22file; }
    std::string getC23file() { return C23file; }
    std::string getC33file() { return C33file; }
    std::string getC44file() { return C44file; }
    std::string getC55file() { return C55file; }
    std::string getC66file() { return C66file; }
    std::string getRfile() { return Rfile; }
    void setC11file(std::string name) { C11file = name; }
    void setC12file(std::string name) { C12file = name; }
    void setC13file(std::string name) { C13file = name; }
    void setC22file(std::string name) { C22file = name; }
    void setC23file(std::string name) { C23file = name; }
    void setC33file(std::string name) { C33file = name; }
    void setC44file(std::string name) { C44file = name; }
    void setC55file(std::string name) { C55file = name; }
    void setC66file(std::string name) { C66file = name; }
    void setRfile(std::string name) { Rfile = name; }
    T getMinVp();  ///< Returns min Vp
    T getMinVs();  ///< Returns min Vs
    T getMinR() {return this->getMin(R); } ///< Returns min R
    T getMaxVp();  ///< Returns max Vp
    T getMaxVs(); ///< Returns max Vs
    T getMaxR() {return this->getMax(R); } ///< Returns max R

    std::shared_ptr<ModelOrtho3D<T>> getLocal(std::shared_ptr<Data3D<T>>, T aperture_x, T aperture_y, bool map);
    std::shared_ptr<ModelOrtho3D<T>> getDomainmodel(std::shared_ptr<Data3D<T>>, T aperturex, T aperturey, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order); ///< Returns a model of a domain

    /** Stagger model functions. 
    It creates the padded Rx, Ry, Rz, L, L2M, M_xz, M_yz, M_xy from the non-padded models R, C11 and C44. 
    */
    void staggerModels();

    /** Create model
    It creates an empty model of C11, C12, C13, C22, C23, C33, C44, C55, C66  and R
    */
    void createModel();
    void createPaddedmodel();

private:
    T *C11;  ///< C11 stiffness
    T *C12;  ///< C12 stiffness
    T *C13;  ///< C13 stiffness
    T *C22;  ///< C22 stiffness
    T *C23;  ///< C23 stiffness
    T *C33;  ///< C33 stiffness
    T *C44;  ///< C44 stiffness
    T *C55;  ///< C55 stiffness
    T *C66;  ///< C66 stiffness
    T *R;   ///< Density 
    T *C11p;  ///< C11 stiffness padded
    T *C12p;  ///< C12 stiffness padded
    T *C13p;  ///< C13 stiffness padded
    T *C22p;  ///< C22 stiffness padded
    T *C23p;  ///< C23 stiffness padded
    T *C33p;  ///< C33 stiffness padded
    T *C44p;  ///< C44 stiffness padded
    T *C55p;  ///< C55 stiffness padded
    T *C66p;  ///< C66 stiffness padded
    T *Rx;  ///< Staggered inverse of density in x (padded)
    T *Ry;  ///< Staggered inverse of density in y (padded)
    T *Rz;  ///< Staggered inverse of density in z (padded)
    std::string C11file; ///< Filename to C11 file
    std::string C12file; ///< Filename to C12 file
    std::string C13file; ///< Filename to C13 file
    std::string C22file; ///< Filename to C22 file
    std::string C23file; ///< Filename to C23 file
    std::string C33file; ///< Filename to C33 file
    std::string C44file; ///< Filename to C44 file
    std::string C55file; ///< Filename to C55 file
    std::string C66file; ///< Filename to C66 file
    std::string Rfile; ///< Filename to density model
};

// =============== 2D POROELASTIC MODEL CLASS =============== //
/** The 2D Poroelastic model class
 *
 */
template<typename T>
class ModelPoroelastic2D: public Model<T> {
public:
    ModelPoroelastic2D();	///< Constructor
    ModelPoroelastic2D(const int _nx, const int _nz, const int lpml, const T _dx, const T _dz, const T _ox, const T _oz, const T _f0, const bool _fs);	///< Constructor
    ModelPoroelastic2D(std::string _Rhofile, std::string _Rhoffile, std::string _Porfile, std::string _Kdfile, std::string _Ksfile,  std::string _Kffile, std::string _Mufile, std::string _Mobfile, std::string _Psifile, const int _lpml, const T _f0, const bool _fs);	///< Constructor
    ~ModelPoroelastic2D();	///< Destructor
    
    // I/O functions
    void readModel();	///< Read a model from file

    void writeRho(); ///< Write only the Rho model to file
    void writeRhof(); ///< Write only the Rhof model to file
    void writePor(); ///< Write only the Por model to file
    void writeKd(); ///< Write only the Kd model to file
    void writeKs(); ///< Write only the Ks model to file
    void writeKf(); ///< Write only the Kf model to file
    void writeMu(); ///< Write only the Mu model to file
    void writeMob(); ///< Write only the Mob model to file
    void writePsi(); ///< Write only the Psi model to file
    void writeModel() { writeRho(); writeRhof(); writePor(); writeKd(); writeKs(); writeKf(); writeMu(); writeMob(); writePsi(); } ///< Write a model to file

    // Get functions
    T *getRho() { return Rho; }	///< Get Rho
    T *getRhof() { return Rhof; }	///< Get Rhof
    T *getPor() { return Por; }	///< Get Por
    T *getKd() { return Kd; }	///< Get Kd
    T *getKs() { return Ks; }	///< Get Ks
    T *getKf() { return Kf; }	///< Get Kf
    T *getMu() { return Mu; }	///< Get Mu
    T *getMob() { return Mob; }	///< Get Mob
    T *getPsi() { return Psi; }	///< Get Psi
    T *getLu() { return Lu; }		///< Get L undrained
    T *getLuM() { return LuM; }		///< Get L2M undrained
    T *getAlpha() { return Alpha; }		///< Get Alpha
    T *getBeta() { return Beta; }		///< Get Beta
    T *getM_xz() { return M_xz; }	///< Get M_xz staggered
    T *getRho_x() { return Rho_x; }		///< Get Rho_x staggered
    T *getRho_z() { return Rho_z; }		///< Get Rho_z staggered
    T *getRhof_x() { return Rhof_x; }		///< Get Rhof_x staggered
    T *getRhof_z() { return Rhof_z; }		///< Get Rhof_z staggered
    T *getMob_x() { return Mob_x; }		///< Get Mob_x 
    T *getMob_z() { return Mob_z; }		///< Get Mob_z 
    T *getPsi_x() { return Psi_x; }		///< Get Psi_x 
    T *getPsi_z() { return Psi_z; }		///< Get Psi_z 
    T getF0() { return f0; } ///< Get the dominant frequency
    std::string getRhofile() { return Rhofile; }
    std::string getRhoffile() { return Rhoffile; }
    std::string getPorfile() { return Porfile; }
    std::string getKdfile() { return Kdfile; }
    std::string getKsfile() { return Ksfile; }
    std::string getKffile() { return Kffile; }
    std::string getMufile() { return Mufile; }
    std::string getMobfile() { return Mobfile; }
    std::string getPsifile() { return Psifile; }
    void setRhofile(std::string name) { Rhofile = name; }
    void setRhoffile(std::string name) { Rhoffile = name; }
    void setPorfile(std::string name) { Porfile = name; }
    void setKdfile(std::string name) { Kdfile = name; }
    void setKsfile(std::string name) { Ksfile = name; }
    void setKffile(std::string name) { Kffile = name; }
    void setMufile(std::string name) { Mufile = name; }
    void setMobfile(std::string name) { Mobfile = name; }
    void setPsifile(std::string name) { Psifile = name; }
    void setF0(T val) { f0 = val; } ///< Set the dominant frequency
    T getMinVp();  ///< Returns min Vp
    T getMinVs();  ///< Returns min Vs
    T getMinR() {return this->getMax(Rho); }  ///< Returns min R
    T getMaxVp(); ///< Returns max Vp
    T getMaxVs(); ///< Returns max Vs
    T getMaxR() {return this->getMax(Rho); } ///< Returns max R
   //
    /** Stagger model functions. 
    It creates the padded models and other poroelastic variables
    */
    std::shared_ptr<ModelPoroelastic2D<T>> getLocal(std::shared_ptr<Data2D<T>>, T aperture, bool map);
    std::shared_ptr<ModelPoroelastic2D<T>> getDomainmodel(std::shared_ptr<Data2D<T>>, T aperture, bool map, const int d, const int nd0, const int nd1, const int order); ///< Returns a model of a domain

    /** Create model
    It creates an empty Poroelastic model
    */
    void createModel();
    void createPaddedmodel();

private:
    T *Rho; 
    T *Rhof; 
    T *Por; 
    T *Kd; 
    T *Ks; 
    T *Kf; 
    T *Mu; 
    T *Mob; 
    T *Psi; 

    T *Lu; // Lame undrained
    T *LuM; // Lame + 2Mu undrained
    T *Alpha; // Biot Willis constant
    T *Beta;   // Beta parameter
    T *M_xz; // Lame Mu  (staggered)
    T *Rho_x;  // Staggered inverse of density in x (padded)
    T *Rho_z;  // Staggered inverse of density in z (padded)
    T *Rhof_x;  // Staggered inverse of fluid density in x (padded)
    T *Rhof_z;  // Staggered inverse of fluid density in z (padded)
    T *Mob_x;  // Mobility Staggered 
    T *Mob_z;  // Mobility Staggered 
    T *Psi_x;  // (1 - O)F - rhof/rho Staggered 
    T *Psi_z;  // (1 - O)F - rhof/rho Staggered 
    T f0; // Center frequency
    std::string Rhofile; ///< Filename to Rho
    std::string Rhoffile; ///< Filename to Rhof
    std::string Porfile; ///< Filename to Por
    std::string Kdfile; ///< Filename to Kd
    std::string Ksfile; ///< Filename to Ks
    std::string Kffile; ///< Filename to Kf
    std::string Mufile; ///< Filename to Mu
    std::string Mobfile; ///< Filename to Mob
    std::string Psifile; ///< Filename to Psi
};


}
#endif //MODEL_H
