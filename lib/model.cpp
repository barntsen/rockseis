// Include statements
#include "model.h"
#include "file.h"

namespace rockseis {


// =============== ABSTRACT MODEL CLASS =============== //
template<typename T>
Model<T>::Model() {
    dim = 0;
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, 1);
    geometry->setN(2, 1);
    geometry->setN(3, 1);
    geometry->setD(1, 1.);
    geometry->setD(2, 1.);
    geometry->setD(3, 1.);
    geometry->setO(1, 0.);
    geometry->setO(2, 0.);
    geometry->setO(3, 0.);
    lpml = 1;
    fs = 0;
}
template<typename T>
Model<T>::Model(const int _dim) {
    dim = _dim;
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, 1);
    geometry->setN(2, 1);
    geometry->setN(3, 1);
    geometry->setD(1, 1.);
    geometry->setD(2, 1.);
    geometry->setD(3, 1.);
    geometry->setO(1, 0.);
    geometry->setO(2, 0.);
    geometry->setO(3, 0.);
    lpml = 1;
    fs = 0;
}

template<typename T>
Model<T>::Model(const int _dim, const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs) {
    dim = _dim;
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, _nx);
    geometry->setN(2, _ny);
    geometry->setN(3, _nz);
    geometry->setD(1, _dx);
    geometry->setD(2, _dy);
    geometry->setD(3, _dz);
    geometry->setO(1, _ox);
    geometry->setO(2, _oy);
    geometry->setO(3, _oz);
    lpml = _lpml;
    fs = _fs;
}

template<typename T>
Model<T>::~Model() {
    // Nothing
}

template <typename T>
void Model<T>::padmodel2d(T *padded, T *model, const int nx, const int ny, const int lpml){
    int ix,iy;
    int nx_pml, ny_pml;
    
    nx_pml = nx + 2*lpml;
    ny_pml = ny + 2*lpml;
    
    Index ind(nx, ny);
    Index ind_pml(nx_pml, ny_pml);
    
    T val;
    // Padding models
    for(ix=0; ix < nx; ix++){
        for(iy=0; iy < ny; iy++){
            val = model[ind(ix,iy)];
            padded[ind_pml(ix+lpml,iy+lpml)]=val;
        }
    }
    
    //Left and right
    for(ix=0; ix < lpml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            padded[ind_pml(ix,iy)]=padded[ind_pml(lpml,iy)];
            padded[ind_pml(nx_pml-lpml+ix,iy)]=padded[ind_pml(nx_pml-lpml-1,iy)];
        }
    }
    
    //Top and bottom
    for(ix=0; ix<nx_pml; ix++){
        for(iy=0; iy<lpml; iy++){
            padded[ind_pml(ix,iy)]=padded[ind_pml(ix,lpml)];
            padded[ind_pml(ix,ny_pml-lpml+iy)]=padded[ind_pml(ix,ny_pml-lpml-1)];
        }
    }
}

template <typename T>
void Model<T>::padmodel3d(T *padded, T *model, const int nx, const int ny, const int nz, const int lpml){
    int ix,iy,iz;
    int nx_pml, ny_pml, nz_pml;
    
    nx_pml = nx + 2*lpml;
    ny_pml = ny + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    Index ind(nx, ny, nz);
    Index ind_pml(nx_pml, ny_pml, nz_pml);
    
    T val;
    // Padding
    // Padding center
    for(ix=0; ix < nx; ix++){
        for(iy=0; iy < ny; iy++){
            for(iz=0; iz < nz; iz++){
                val = model[ind(ix,iy,iz)];
                padded[ind_pml(ix+lpml,iy+lpml,iz+lpml)]=val;
            }
        }
    }
    
    //Left and right
    for(ix=0; ix < lpml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                padded[ind_pml(ix,iy,iz)]=padded[ind_pml(lpml,iy,iz)];
                padded[ind_pml(nx_pml-lpml+ix,iy,iz)]=padded[ind_pml(nx_pml-lpml-1,iy,iz)];
            }
        }
    }
    
    //Top and bottom
    for(ix=0; ix<nx_pml; ix++){
        for(iy=0; iy<ny_pml; iy++){
            for(iz=0; iz<lpml; iz++){
                padded[ind_pml(ix,iy,iz)]=padded[ind_pml(ix,iy,lpml)];
                padded[ind_pml(ix,iy,nz_pml-lpml+iz)]=padded[ind_pml(ix,iy,nz_pml-lpml-1)];
            }
        }
    }
    
    //Front and back
    for(ix=0; ix<nx_pml; ix++){
        for(iy=0; iy<lpml; iy++){
            for(iz=0; iz<nz_pml; iz++){
                padded[ind_pml(ix,iy,iz)]=padded[ind_pml(ix,lpml,iy)];
                padded[ind_pml(ix,ny_pml-lpml+iy,iz)]=padded[ind_pml(ix,ny_pml-lpml-1,iz)];
            }
        }
    }
}

template <typename T>
void Model<T>::staggermodel_x(T *model, const int nx, const int ny, const int nz){
    int ix,iy,iz;
    
    Index ind(nx, ny, nz);
    
    T avg;
    //Stagger in x
    for(iz=0; iz<nz; iz++){
        for(iy=0; iy<ny; iy++){
            for(ix=0; ix<nx-1; ix++){
                avg = 2.0/(model[ind(ix+1,iy,iz)] + model[ind(ix,iy,iz)]);
                model[ind(ix,iy,iz)] = avg;
            }
            model[ind(nx-1,iy,iz)] = model[ind(nx-2,iy,iz)];
        }
    }
}

template <typename T>
void Model<T>::staggermodel_y(T *model, const int nx, const int ny, const int nz){
    int ix,iy,iz;
    
    Index ind(nx, ny, nz);
    
    T avg;
    //Stagger in y
    for(iz=0; iz<nz; iz++){
        for(ix=0; ix<nx; ix++){
            for(iy=0; iy<ny-1; iy++){
                avg = 2.0/(model[ind(ix,iy+1,iz)] + model[ind(ix,iy,iz)]);
                model[ind(ix,iy,iz)] = avg;
            }
            model[ind(ix,ny-1,iz)] = model[ind(ix,ny-2,iz)];
        }
    }
}

template <typename T>
void Model<T>::staggermodel_z(T *model, const int nx, const int ny, const int nz){
    int ix,iy,iz;
    
    Index ind(nx, ny, nz);
    
    T avg;
    //Stagger in z
    for(ix=0; ix<nx; ix++){
        for(iy=0; iy<ny; iy++){
            for(iz=0; iz<nz-1; iz++){
                avg = 2.0/(model[ind(ix,iy,iz+1)] + model[ind(ix,iy,iz)]);
                model[ind(ix,iy,iz)] = avg;
            }
            model[ind(ix,iy,nz-1)] = model[ind(ix,iy,nz-2)];
        }
    }
}


// =============== 2D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(): Model<T>(2) {
    /* Allocate variables */
    int nx, nz, nx_pml, nz_pml, lpml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    nx_pml = nx +2*lpml;
    nz_pml = nz +2*lpml;
    
    Vp = (T *) malloc(nx*nz*sizeof(T));
    R = (T *) malloc(nx*nz*sizeof(T));
    L = (T *) malloc(nx_pml*nz_pml*sizeof(T));
    Rx = (T *) malloc(nx_pml*nz_pml*sizeof(T));
    Rz = (T *) malloc(nx_pml*nz_pml*sizeof(T));
    
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    int nx, nz, nx_pml, nz_pml, lpml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    nx_pml = nx +2*lpml;
    nz_pml = nz +2*lpml;
    
    /* Allocate variables */
    Vp = (T *) calloc(_nx*_nz,sizeof(T));
    R = (T *) calloc(_nx*_nz,sizeof(T));
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz, nx_pml, nz_pml;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelAcoustic2D::Error reading from density file.");
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelAcoustic2D::Geometries in Vp and Density model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelAcoustic2D::Numerical precision mismatch in Vp and Density model files.");
    }

    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelAcoustic2D::Numerical precision in Vp and Density model files mismatch with constructor.");
    }
    
    // Read geometry from file
    nx = Fvp->getN(1);
    dx = (T) Fvp->getD(1);
    ox = (T) Fvp->getO(1);
    nz = Fvp->getN(3);
    dz = (T) Fvp->getD(3);
    oz = (T) Fvp->getO(3);
    
    // Close files
    Fvp->close();
    Frho->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);

    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();
    
    /* Allocate variables */
    Vp = (T *) calloc(nx*nz,sizeof(T));
    R = (T *) calloc(nx*nz,sizeof(T));
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
ModelAcoustic2D<T>::~ModelAcoustic2D() {
    // Freeing all variables
    free(Vp);
    free(R);
    free(L);
    free(Rx);
    free(Rz);
}

template<typename T>
void ModelAcoustic2D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::readModel : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::readModel : Error reading from Density file.");
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();

    Vp = this->getVp();
    Fvp->read(Vp, nx*nz);
    Fvp->close();

    R = this->getR();
    Frho->read(R, nx*nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic2D<T>::writeModel() {
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile.c_str());

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvp->setN(1,nx);
    Fvp->setN(3,nz);
    Fvp->setD(1,dx);
    Fvp->setD(3,dz);
    Fvp->setO(1,ox);
    Fvp->setO(3,oz);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*nz);
    Fvp->close();

    Frho->setN(1,nx);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(3,oz);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic2D<T>::staggerModels(){
    int ix,iz;
    int nx, nz, lpml, nx_pml, nz_pml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();
    
    Index ind(nx, nz);
    Index ind_pml(nx_pml, nz_pml);
    
    // Padding
    this->padmodel2d(Rx, R, nx, nz, lpml);
    this->padmodel2d(Rz, R, nx, nz, lpml);
    this->padmodel2d(L, Vp, nx, nz, lpml);
    
    // Computing moduli
    T _rho, _vp;
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            _vp=L[ind_pml(ix,iz)];
            _rho=Rx[ind_pml(ix,iz)];
            L[ind_pml(ix,iz)]=_rho*_vp*_vp;
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, 1, nz_pml);
    this->staggermodel_z(Rz, nx_pml, 1, nz_pml);
    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            Rx[ind_pml(ix,lpml)] *= 2.0;
            L[ind_pml(ix,lpml)] *= 0.0;
        }
    }
}

// =============== 3D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic3D<T>::ModelAcoustic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, _fs) {
    int nx, ny, nz, nx_pml, ny_pml, nz_pml, lpml;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    lpml = this->getLpml();
    nx_pml = nx +2*lpml;
    ny_pml = ny +2*lpml;
    nz_pml = nz +2*lpml;
    
    /* Allocate variables */
    Vp = (T *) calloc(_nx*_ny*_nz,sizeof(T));
    R = (T *) calloc(_nx*_ny*_nz,sizeof(T));
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
ModelAcoustic3D<T>::ModelAcoustic3D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz, nx_pml, ny_pml, nz_pml;
    T dx, dy, dz;
    T ox, oy, oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelAcoustic3D::Error reading from density file.");
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelAcoustic3D::Geometries in Vp and Density model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelAcoustic3D::Numerical precision mismatch in Vp and Density model files.");
    }
    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelAcoustic3D::Numerical precision in Vp and Density model files mismatch with constructor.");
    }
 
    
    // Read geometry from file
    nx = Fvp->getN(1);
    dx = (T) Fvp->getD(1);
    ox = (T) Fvp->getO(1);
    ny = Fvp->getN(2);
    dy = (T) Fvp->getD(2);
    oy = (T) Fvp->getO(2);
    nz = Fvp->getN(3);
    dz = (T) Fvp->getD(3);
    oz = (T) Fvp->getO(3);
    
    // Close files
    Fvp->close();
    Frho->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNy(ny);
    this->setDy(dy);
    this->setOy(oy);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);

    nx_pml = this->getNx_pml();
    ny_pml = this->getNy_pml();
    nz_pml = this->getNz_pml();
    
    /* Allocate variables */
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}


template<typename T>
void ModelAcoustic3D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::readModel : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::readModel : Error reading from Density file.");
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    Vp = this->getVp();
    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

    R = this->getR();
    Frho->read(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic3D<T>::writeModel() {
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile.c_str());

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int ny = this->getNy();
    T dy = this->getDy();
    T oy = this->getOy();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvp->setN(1,nx);
    Fvp->setN(2,ny);
    Fvp->setN(3,nz);
    Fvp->setD(1,dx);
    Fvp->setD(2,dy);
    Fvp->setD(3,dz);
    Fvp->setO(1,ox);
    Fvp->setO(2,oy);
    Fvp->setO(3,oz);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*ny*nz);
    Fvp->close();

    Frho->setN(1,nx);
    Frho->setN(2,ny);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(2,dy);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(2,oy);
    Frho->setO(3,oz);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
ModelAcoustic3D<T>::~ModelAcoustic3D() {
    // Freeing all variables
    free(Vp);
    free(R);
    free(L);
    free(Rx);
    free(Ry);
    free(Rz);
}

template<typename T>
void ModelAcoustic3D<T>::staggerModels(){
    int ix,iy,iz;
    int nx, ny, nz, lpml, nx_pml, ny_pml, nz_pml;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = nx + 2*lpml;
    ny_pml = ny + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    Index ind_pml(nx_pml, ny_pml, nz_pml);
    
    // Padding
    this->padmodel3d(Rx, R, nx, ny, nz, lpml);
    this->padmodel3d(Ry, R, nx, ny, nz, lpml);
    this->padmodel3d(Rz, R, nx, ny, nz, lpml);
    this->padmodel3d(L, Vp, nx, ny, nz, lpml);
    
    // Computing moduli
    T _rho, _vp;
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                _vp=L[ind_pml(ix,iy,iz)];
                _rho=Rx[ind_pml(ix,iy,iz)];
                L[ind_pml(ix,iy,iz)]=_rho*_vp*_vp;
            }
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(Ry, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(Rz, nx_pml, ny_pml, nz_pml);
    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy<ny_pml; iy++){
                Rx[ind_pml(ix,iy,lpml)] *= 2.0;
                Ry[ind_pml(ix,iy,lpml)] *= 2.0;
                L[ind_pml(ix,iy,lpml)] *= 0.0;
            }
        }
    }
}



// =============== 2D ELASTIC MODEL CLASS =============== //
template<typename T>
ModelElastic2D<T>::ModelElastic2D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelElastic2D<T>::ModelElastic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    int nx, nz, nx_pml, nz_pml, lpml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    nx_pml = nx +2*lpml;
    nz_pml = nz +2*lpml;
    
    /* Allocate variables */
    Vp = (T *) calloc(_nx*_nz,sizeof(T));
    Vs = (T *) calloc(_nx*_nz,sizeof(T));
    R = (T *) calloc(_nx*_nz,sizeof(T));
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}

template<typename T>
ModelElastic2D<T>::ModelElastic2D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz, nx_pml, nz_pml;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::Error reading from Vp file.");
	    exit(1);
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic2D::Error reading from Vs file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic2D::Error reading from density file. ");
        exit(1);
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Fvs) != 0)
    {
        rs_error("ModelElastic2D::Geometries in Vp and Vs model files do not match.");
    }

    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelElastic2D::Geometries in Vp and Density model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelElastic2D::Numerical precision mismatch in Vp and Density model files.");
    }
    if(Fvp->getData_format() != Fvs->getData_format())
    {
        rs_error("ModelElastic2D::Numerical precision mismatch in Vp and Vs model files.");
    }
    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelElastic2D::Numerical precision in Vp, Vs and Density model files mismatch with constructor.");
    }
 

    
    // Read geometry from file
    nx = Fvp->getN(1);
    dx = (T) Fvp->getD(1);
    ox = (T) Fvp->getO(1);
    nz = Fvp->getN(3);
    dz = (T) Fvp->getD(3);
    oz = (T) Fvp->getO(3);
    
    // Close files
    Fvp->close();
    Fvs->close();
    Frho->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);

    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();
    
    /* Allocate variables */
    Vp = (T *) calloc(nx*nz,sizeof(T));
    Vs = (T *) calloc(nx*nz,sizeof(T));
    R = (T *) calloc(nx*nz,sizeof(T));
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
}


template<typename T>
void ModelElastic2D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Vp file. ");
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Vs file. ");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Density file. ");
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();

    Vp = this->getVp();
    Fvp->read(Vp, nx*nz);
    Fvp->close();

    Vs = this->getVs();
    Fvs->read(Vs, nx*nz);
    Fvs->close();

    R = this->getR();
    Frho->read(R, nx*nz);
    Frho->close();
}

template<typename T>
void ModelElastic2D<T>::writeModel() {
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    Fvs->output(Vsfile.c_str());
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile.c_str());

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvp->setN(1,nx);
    Fvp->setN(3,nz);
    Fvp->setD(1,dx);
    Fvp->setD(3,dz);
    Fvp->setO(1,ox);
    Fvp->setO(3,oz);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*nz);
    Fvp->close();

    Fvs->setN(1,nx);
    Fvs->setN(3,nz);
    Fvs->setD(1,dx);
    Fvs->setD(3,dz);
    Fvs->setO(1,ox);
    Fvs->setO(3,oz);
    Fvs->setData_format(sizeof(T));
    Fvs->writeHeader();
    Vs = this->getVs();
    Fvs->write(Vs, nx*nz);
    Fvs->close();

    Frho->setN(1,nx);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(3,oz);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*nz);
    Frho->close();

}


template<typename T>
void ModelElastic2D<T>::staggerModels(){
    int ix,iz;
    int nx, nz, lpml, nx_pml, nz_pml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = nx + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    Index ind(nx, nz);
    Index ind_pml(nx_pml, nz_pml);
    
    // Padding
    this->padmodel2d(Rx, R, nx, nz, lpml);
    this->padmodel2d(Rz, R, nx, nz, lpml);
    this->padmodel2d(L, Vp, nx, nz, lpml);
    this->padmodel2d(M, Vs, nx, nz, lpml);
    
    // Computing moduli
    T _rho, _vp, _vs;
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            _vp=L[ind_pml(ix,iz)];
            _vs=M[ind_pml(ix,iz)];
            _rho=Rx[ind_pml(ix,iz)];
            L2M[ind_pml(ix,iz)]=_rho*_vp*_vp;
            M[ind_pml(ix,iz)]=_rho*_vs*_vs;
            L[ind_pml(ix,iz)]=_rho*_vp*_vp - 2.0*_rho*_vs*_vs; 
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, 1, nz_pml); 
    this->staggermodel_z(Rz, nx_pml, 1, nz_pml); 
    
    this->staggermodel_z(M, nx_pml, 1, nz_pml); 
    this->staggermodel_x(M, nx_pml, 1, nz_pml); 
    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            Rx[ind_pml(ix,lpml)] *= 2.0;
            L[ind_pml(ix,lpml)] *= 0.0;
            L2M[ind_pml(ix,lpml)] *= 0.0;
        }
    }
}

template<typename T>
ModelElastic2D<T>::~ModelElastic2D() {
    free(Vp);
    free(Vs);
    free(R);
    free(L);
    free(L2M);
    free(M);
    free(Rx);
    free(Rz);
}


// =============== 3D ELASTIC MODEL CLASS =============== //
template<typename T>
ModelElastic3D<T>::ModelElastic3D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelElastic3D<T>::ModelElastic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, _fs) {
    int nx, ny, nz, nx_pml, ny_pml, nz_pml, lpml;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    lpml = this->getLpml();
    nx_pml = nx +2*lpml;
    ny_pml = ny +2*lpml;
    nz_pml = nz +2*lpml;
    
    /* Allocate variables */
    Vp = (T *) calloc(_nx*_ny*_nz,sizeof(T));
    Vs = (T *) calloc(_nx*_ny*_nz,sizeof(T));
    R = (T *) calloc(_nx*_ny*_nz,sizeof(T));
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_yz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
ModelElastic3D<T>::ModelElastic3D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz, nx_pml, ny_pml, nz_pml;
    T dx, dy, dz;
    T ox, oy, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::Error reading from Vp file. ");
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::Error reading from Vs file. ");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::Error reading from density file. ");
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Fvs) != 0)
    {
        rs_error("ModelElastic3D::Geometries in Vp and Vs model files do not match.");
    }

    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelElastic3D::Geometries in Vp and Density model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelElastic3D::Numerical precision mismatch in Vp and Density model files.");
    }
    if(Fvp->getData_format() != Fvs->getData_format())
    {
        rs_error("ModelElastic3D::Numerical precision mismatch in Vp and Vs model files.");
    }
    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelElastic3D::Numerical precision in Vp, Vs and Density model files mismatch with constructor.");
    }

    // Read geometry from file
    nx = Fvp->getN(1);
    dx = (T) Fvp->getD(1);
    ox = (T) Fvp->getO(1);
    ny = Fvp->getN(2);
    dy = (T) Fvp->getD(2);
    oy = (T) Fvp->getO(2);
    nz = Fvp->getN(3);
    dz = (T) Fvp->getD(3);
    oz = (T) Fvp->getO(3);
    
    // Close files
    Fvp->close();
    Fvs->close();
    Frho->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNy(ny);
    this->setDy(dy);
    this->setOy(oy);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);

    nx_pml = this->getNx_pml();
    ny_pml = this->getNy_pml();
    nz_pml = this->getNz_pml();
    
    /* Allocate variables */
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    Vs = (T *) calloc(nx*ny*nz,sizeof(T));
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_yz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
}

template<typename T>
void ModelElastic3D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Vp file. ");
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Vs file. ");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Density file. ");
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    Vp = this->getVp();
    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

    Vs = this->getVs();
    Fvs->read(Vs, nx*ny*nz);
    Fvs->close();

    R = this->getR();
    Frho->read(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
void ModelElastic3D<T>::writeModel() {
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    Fvs->output(Vsfile.c_str());
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile.c_str());

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int ny = this->getNy();
    T dy = this->getDy();
    T oy = this->getOy();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvp->setN(1,nx);
    Fvp->setN(2,ny);
    Fvp->setN(3,nz);
    Fvp->setD(1,dx);
    Fvp->setD(2,dy);
    Fvp->setD(3,dz);
    Fvp->setO(1,ox);
    Fvp->setO(2,oy);
    Fvp->setO(3,oz);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*ny*nz);
    Fvp->close();

    Fvs->setN(1,nx);
    Fvs->setN(2,ny);
    Fvs->setN(3,nz);
    Fvs->setD(1,dx);
    Fvs->setD(2,dy);
    Fvs->setD(3,dz);
    Fvs->setO(1,ox);
    Fvs->setO(2,oy);
    Fvs->setO(3,oz);
    Fvs->setData_format(sizeof(T));
    Fvs->writeHeader();
    Vs = this->getVs();
    Fvs->write(Vs, nx*ny*nz);
    Fvs->close();

    Frho->setN(1,nx);
    Frho->setN(2,ny);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(2,dy);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(2,oy);
    Frho->setO(3,oz);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
void ModelElastic3D<T>::staggerModels(){
    int ix,iy,iz;
    int nx, ny, nz, lpml, nx_pml, ny_pml, nz_pml;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = nx + 2*lpml;
    ny_pml = ny + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    Index ind_pml(nx_pml, ny_pml, nz_pml);
    
    // Padding
    this->padmodel3d(Rx, R, nx, ny, nz, lpml);
    this->padmodel3d(Ry, R, nx, ny, nz, lpml);
    this->padmodel3d(Rz, R, nx, ny, nz, lpml);
    this->padmodel3d(L, Vp, nx, ny, nz, lpml);
    this->padmodel3d(M_xz, Vs, nx, ny, nz, lpml);
    
    // Computing moduli
    T _rho, _vp, _vs;
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                _vp=L[ind_pml(ix,iy,iz)];
                _vs=M_xz[ind_pml(ix,iy,iz)];
                _rho=Rx[ind_pml(ix,iy,iz)];
                L2M[ind_pml(ix,iy,iz)]=_rho*_vp*_vp;
                L[ind_pml(ix,iy,iz)]=_rho*_vp*_vp - 2.0*_rho*_vs*_vs;
                M_xz[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                M_yz[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                M_xy[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
            }
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(Ry, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(Rz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_x(M_xz, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(M_xz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_y(M_yz, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(M_yz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_x(M_xy, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(M_xy, nx_pml, ny_pml, nz_pml);
    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy<ny_pml; iy++){
                Rx[ind_pml(ix,iy,lpml)] *= 2.0;
                Ry[ind_pml(ix,iy,lpml)] *= 2.0;
                L[ind_pml(ix,iy,lpml)] *= 0.0;
            }
        }
    }
}

template<typename T>
ModelElastic3D<T>::~ModelElastic3D() {
    free(Vp);
    free(Vs);
    free(R);
    free(L);
    free(L2M);
    free(M_xz);
    free(M_yz);
    free(M_xy);
    free(Rx);
    free(Ry);
    free(Rz);
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class ModelAcoustic2D<float>;
template class ModelAcoustic3D<float>;
template class ModelElastic2D<float>;
template class ModelElastic3D<float>;

template class ModelAcoustic2D<double>;
template class ModelAcoustic3D<double>;
template class ModelElastic2D<double>;
template class ModelElastic3D<double>;


}
