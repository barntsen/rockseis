// Include statements
#include "model.h"

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
    fs = false;
    realized=false;
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
    fs = false;
    realized=false;
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
    realized=false;
}

template<typename T>
Model<T>::~Model() {
    // Nothing
}

template <typename T>
void Model<T>::padmodel1d(T *padded, T *model, const int nx, const int lpml){
    int ix;
    int nx_pml;
    
    nx_pml = nx + 2*lpml;
    
    T val;
    // Padding models
    for(ix=0; ix < nx; ix++){
            val = model[ix];
            padded[ix+lpml]=val;
    }
    
    for(ix=0; ix < lpml; ix++){
            padded[ix]=padded[lpml];
            padded[nx_pml-lpml+ix]=padded[nx_pml-lpml-1];
    }
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

template <typename T>
T Model<T>::getMax(T *model){
    int i;

    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    T max = model[0];
    for(i=1; i<nx*ny*nz; i++){
        if(model[i] > max) max = model[i];
    }
    return max;
}

template <typename T>
T Model<T>::getMin(T *model){
    int i;

    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    T min = model[0];
    for(i=1; i<nx*ny*nz; i++){
        if(model[i] < min) min = model[i];
    }
    return min;
}

// =============== 1D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic1D<T>::ModelAcoustic1D(): Model<T>(2) {
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
    
}

template<typename T>
ModelAcoustic1D<T>::ModelAcoustic1D(const int _nz, const int _lpml, const T _dz, const T _oz, const bool _fs): Model<T>(2, 1, 1, _nz,  _lpml, 1.0, 1.0, _dz, 0.0, 0.0, _oz, _fs) {
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelAcoustic1D<T>::ModelAcoustic1D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nz;
    T dz;
    T oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic1D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelAcoustic1D::Error reading from density file: ", Rfile);
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelAcoustic1D::Geometries in Vp and Density model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelAcoustic1D::Numerical precision mismatch in Vp and Density model files.");
    }

    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelAcoustic1D::Numerical precision in Vp and Density model files mismatch with constructor.");
    }
    
    // Read geometry from file
    nz = Fvp->getN(3);
    dz = (T) Fvp->getD(3);
    oz = (T) Fvp->getO(3);
    
    // Close files
    Fvp->close();
    Frho->close();

    // Store geometry in model class
    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);

    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelAcoustic1D<T>::~ModelAcoustic1D() {
    // Freeing all variables
    free(Vp);
    free(R);
    free(L);
    free(Rz);
}

template<typename T>
void ModelAcoustic1D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic1D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic1D::readModel : Error reading from Density file: ", Rfile );
    }

    // Read models
    int nz = this->getNz();

    // Reallocate variables to correct size
    free(Vp); free(R);
    Vp = (T *) calloc(nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic1D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic1D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Vp = this->getVp();
    Fvp->read(Vp, nz);
    Fvp->close();

    R = this->getR();
    Frho->read(R, nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic1D<T>::writeModel() {
    if(!this->getRealized()) {
        rs_error("ModelAcoustic1D::writeModel: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile);
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile);

    // Write models
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvp->setN(3,nz);
    Fvp->setD(3,dz);
    Fvp->setO(3,oz);
    Fvp->setType(REGULAR);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nz, 0);
    Fvp->close();

    Frho->setN(3,nz);
    Frho->setD(3,dz);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nz, 0);
    Frho->close();
}

template<typename T>
void ModelAcoustic1D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelAcoustic1D::staggerModels: Model is not allocated.");
    }
    int iz;
    int nz, lpml, nz_pml;
    nz = this->getNz();
    lpml = this->getLpml();
    nz_pml = this->getNz_pml();
    
    // Reallocate necessary variables 
    free(L); free(Rz);
    L = (T *) calloc(nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic1D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelAcoustic1D::staggerModels: Failed to allocate memory.");
    
    // Padding
    this->padmodel1d(Rz, R, nz, lpml);
    this->padmodel1d(L, Vp, nz, lpml);
    
    // Computing moduli
    T _rho, _vp;
    for(iz=0; iz < nz_pml; iz++){
        _vp=L[iz];
        _rho=Rz[iz];
        L[iz]=_rho*_vp*_vp;
    }

    // Staggering using arithmetic average
    this->staggermodel_z(Rz, 1, 1, nz_pml);
    
    // In case of free surface
    if(this->getFs()){
            L[lpml] *= 0.0;
    }
}


// =============== 2D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(): Model<T>(2) {
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
    
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 0.0, _oz, _fs) {
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelAcoustic2D::Error reading from density file: ", Rfile);
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

    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
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
	    rs_error("ModelAcoustic2D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();

    // Reallocate variables to correct size
    free(Vp); free(R);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic2D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Vp = this->getVp();
    Fvp->read(Vp, nx*nz);
    Fvp->close();

    R = this->getR();
    Frho->read(R, nx*nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic2D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelAcoustic2D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile);
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
    Fvp->setType(REGULAR);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*nz, 0);
    Fvp->close();
}

template<typename T>
void ModelAcoustic2D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelAcoustic2D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Frho->setN(1,nx);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*nz, 0);
    Frho->close();
}

template<typename T>
void ModelAcoustic2D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelAcoustic2D::staggerModels: Model is not allocated.");
    }
    int ix,iz;
    int nx, nz, lpml, nx_pml, nz_pml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();
    
    Index ind(nx, nz);
    Index ind_pml(nx_pml, nz_pml);

    // Reallocate necessary variables 
    free(L); free(Rx); free(Rz);
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic2D::staggerModels: Failed to allocate memory.");
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelAcoustic2D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelAcoustic2D::staggerModels: Failed to allocate memory.");
    
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

template<typename T>
void ModelAcoustic2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate Vp and R */
    free(Vp); free(R);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic2D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<rockseis::ModelAcoustic2D<T>> ModelAcoustic2D<T>::getLocal(std::shared_ptr<rockseis::Data2D<T>> data, T aperture, bool map) {
    std::shared_ptr<rockseis::ModelAcoustic2D<T>> local;
    /* Get source or receiver min and max positions */
    Point2D<T> *coords;
    size_t ntr = data->getNtrace();
    T min, max; 
    if(map == SMAP){
        coords = (data->getGeom())->getScoords();
        min = coords[0].x;
        max = coords[0].x;
        for (int i=1; i < ntr; i++){
            if(coords[i].x < min) min = coords[i].x;
            if(coords[i].x > max) max = coords[i].x;
        }
    }else{
        coords = (data->getGeom())->getGcoords();
        min = coords[0].x;
        max = coords[0].x;
        for (size_t i=1; i < ntr; i++){
            if(coords[i].x < min) min = coords[i].x;
            if(coords[i].x > max) max = coords[i].x;
        }
    }

    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
	/* Determine grid positions and sizes */
    size_t size = (size_t) (rintf((max-min + aperture)/dx) + 1);
    if( size % 2 == 0 ) size++; // Get odd size due to symmetry
    off_t start = (off_t) (rintf((min - ox)/dx) - (size - 1)/2); 

    /* Create local model */
    local = std::make_shared<rockseis::ModelAcoustic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic2d::getlocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic2d::getlocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getLocal : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    rockseis::Index l2d(size,nz);
    rockseis::Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic2D::getLocal: Error reading from vp file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic2D::getLocal: Error reading from rho file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
        }
    }

    /* Free traces */
    free(vptrace);
    free(rhotrace);

    return local;
}

// =============== 3D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic3D<T>::ModelAcoustic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, _fs) {
    
    /* Allocate minimally the variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelAcoustic3D<T>::ModelAcoustic3D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz;
    T dx, dy, dz;
    T ox, oy, oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelAcoustic3D::Error reading from density file: ", Rfile);
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

    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
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
	    rs_error("ModelAcoustic3D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    // Allocate models before reading
    free(Vp); free(R);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic3D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic3D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Vp = this->getVp();
    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

    R = this->getR();
    Frho->read(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
void ModelAcoustic3D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelAcoustic3D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());

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
    Fvp->setType(REGULAR);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*ny*nz, 0);
    Fvp->close();
}

template<typename T>
void ModelAcoustic3D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelAcoustic3D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
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

    Frho->setN(1,nx);
    Frho->setN(2,ny);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(2,dy);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(2,oy);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*ny*nz, 0);
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
    if(!this->getRealized()) {
        rs_error("ModelAcoustic3D::staggerModels: Model is not allocated.");
    }
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

    // Reallocate necessary variables 
    free(L); free(Rx); free(Ry); free(Rz);
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic3D::staggerModels: Failed to allocate memory.");
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelAcoustic3D::staggerModels: Failed to allocate memory.");
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    if(Ry == NULL) rs_error("ModelAcoustic3D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelAcoustic3D::staggerModels: Failed to allocate memory.");
    
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

template<typename T>
void ModelAcoustic3D<T>::createModel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Vp and R */
    free(Vp); free(R);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic3D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic3D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<rockseis::ModelAcoustic3D<T>> ModelAcoustic3D<T>::getLocal(std::shared_ptr<rockseis::Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<rockseis::ModelAcoustic3D<T>> local;
    /* Get source or receiver min and max positions */
    Point3D<T> *coords;
    size_t ntr = data->getNtrace();
    T min_x, max_x; 
    T min_y, max_y; 
    if(map == SMAP){
        coords = (data->getGeom())->getScoords();
        min_x = coords[0].x;
        max_x = coords[0].x;
        min_y = coords[0].y;
        max_y = coords[0].y;
        for (int i=1; i < ntr; i++){
            if(coords[i].x < min_x) min_x = coords[i].x;
            if(coords[i].x > max_x) max_x = coords[i].x;
            if(coords[i].y < min_y) min_y = coords[i].y;
            if(coords[i].y > max_y) max_y = coords[i].y;
        }
    }else{
        coords = (data->getGeom())->getGcoords();
        min_x = coords[0].x;
        max_x = coords[0].x;
        min_y = coords[0].y;
        max_y = coords[0].y;
        for (size_t i=1; i < ntr; i++){
            if(coords[i].x < min_x) min_x = coords[i].x;
            if(coords[i].x > max_x) max_x = coords[i].x;
            if(coords[i].y < min_y) min_y = coords[i].y;
            if(coords[i].y > max_y) max_y = coords[i].y;
        }
    }

    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();

	/* Determine grid positions and sizes */
    size_t size_x = (size_t) (rintf((max_x-min_x + aperture_x)/dx) + 1);
    if( size_x % 2 == 0 ) size_x++; // Get odd size due to symmetry
    off_t start_x = (off_t) (rintf((min_x - ox)/dx) - (size_x - 1)/2); 

    size_t size_y = (size_t) (rintf((max_y-min_y + aperture_y)/dy) + 1);
    if( size_y % 2 == 0 ) size_y++; // Get odd size due to symmetry
    off_t start_y = (off_t) (rintf((min_y - oy)/dy) - (size_y - 1)/2); 

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<rockseis::ModelAcoustic3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic3D::getlocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic3D::getlocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    rockseis::Index l3d(size_x, size_y, nz);
    rockseis::Index f3d(nx, ny, nz);
    rockseis::Index l2d(nx, ny);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f3d(0, 0, i1)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from vp file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from rho file");
        for(size_t i3=0; i3<size_y; i3++) {
            lpos_y = j + i3;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<size_x; i2++) {
                lpos_x = i + i2;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
            }
        }
    }

    /* Free traces */
    free(vptrace);
    free(rhotrace);

    return local;
}



// =============== 2D ELASTIC MODEL CLASS =============== //
template<typename T>
ModelElastic2D<T>::ModelElastic2D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelElastic2D<T>::ModelElastic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

}

template<typename T>
ModelElastic2D<T>::ModelElastic2D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::Error reading from Vp file: ", Vpfile);
	    exit(1);
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic2D::Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic2D::Error reading from density file: ", Rfile);
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

    
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

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
	    rs_error("ModelElastic2D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();
    
    /* Reallocate Vp, Vs  and R */
    free(Vp); free(Vs); free(R);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelElastic2D::readModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelElastic2D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelElastic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

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
void ModelElastic2D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelElastic2D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    Frho->output(Rfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Frho->setN(1,nx);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*nz, 0);
    Frho->close();
}


template<typename T>
void ModelElastic2D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelElastic2D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile);

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
    Fvp->setType(REGULAR);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*nz, 0);
    Fvp->close();
}


template<typename T>
void ModelElastic2D<T>::writeVs() {
    if(!this->getRealized()) {
        rs_error("ModelElastic2D::writeVs: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    Fvs->output(Vsfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fvs->setN(1,nx);
    Fvs->setN(3,nz);
    Fvs->setD(1,dx);
    Fvs->setD(3,dz);
    Fvs->setO(1,ox);
    Fvs->setO(3,oz);
    Fvs->setType(REGULAR);
    Fvs->setData_format(sizeof(T));
    Fvs->writeHeader();
    Vs = this->getVs();
    Fvs->write(Vs, nx*nz, 0);
    Fvs->close();
}

template<typename T>
void ModelElastic2D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelElastic2D::staggerModels: Model is not allocated.");
    }
    int ix,iz;
    int nx, nz, lpml, nx_pml, nz_pml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = nx + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    Index ind(nx, nz);
    Index ind_pml(nx_pml, nz_pml);
    
    // Reallocate necessary variables 
    free(L); free(L2M); free(M); free(Rx); free(Rz);
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelElastic2D::staggerModels: Failed to allocate memory.");
    L2M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(L2M == NULL) rs_error("ModelElastic2D::staggerModels: Failed to allocate memory.");
    M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(M == NULL) rs_error("ModelElastic2D::staggerModels: Failed to allocate memory.");
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelElastic2D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelElastic2D::staggerModels: Failed to allocate memory.");

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
void ModelElastic2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate Vp, Vs and R */
    free(Vp); free(Vs); free(R);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelElastic2D::createModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelElastic2D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelElastic2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<rockseis::ModelElastic2D<T>> ModelElastic2D<T>::getLocal(std::shared_ptr<rockseis::Data2D<T>> data, T aperture, bool map) {

    std::shared_ptr<rockseis::ModelElastic2D<T>> local;
    /* Get source or receiver min and max positions */
    Point2D<T> *coords;
    size_t ntr = data->getNtrace();
    T min, max; 
    if(map == SMAP){
        coords = (data->getGeom())->getScoords();
        min = coords[0].x;
        max = coords[0].x;
        for (int i=1; i < ntr; i++){
            if(coords[i].x < min) min = coords[i].x;
            if(coords[i].x > max) max = coords[i].x;
        }
    }else{
        coords = (data->getGeom())->getGcoords();
        min = coords[0].x;
        max = coords[0].x;
        for (size_t i=1; i < ntr; i++){
            if(coords[i].x < min) min = coords[i].x;
            if(coords[i].x > max) max = coords[i].x;
        }
    }

    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
	/* Determine grid positions and sizes */
    size_t size = (size_t) (rintf((max-min + aperture)/dx) + 1);
    if( size % 2 == 0 ) size++; // Get odd size due to symmetry
    off_t start = (off_t) (rintf((min - ox)/dx) - (size - 1)/2); 

    /* Create local model */
    local = std::make_shared<rockseis::ModelElastic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic2d::getlocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic2d::getlocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic2d::getlocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Vs file.");
    }

    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    rockseis::Index l2d(size,nz);
    rockseis::Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelElastic2D::getLocal: Error reading from vp file");
        Fvs->read(vstrace, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic2D::getLocal: Error reading from vs file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelElastic2D::getLocal: Error reading from rho file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            Vs[l2d(i2,i1)] = vstrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(rhotrace);

    return local;
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
    
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    M_yz = (T *) calloc(1,1);
    M_xy = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelElastic3D<T>::ModelElastic3D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz;
    T dx, dy, dz;
    T ox, oy, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;

    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::Error reading from Vp file:  ", _Vpfile);
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::Error reading from Vs file: ", _Vsfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::Error reading from density file: ", _Rfile);
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

    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    M_yz = (T *) calloc(1,1);
    M_xy = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
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
	    rs_error("ModelElastic3D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Vp, Vs  and R */
    free(Vp); free(Vs); free(R);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelElastic3D::readModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelElastic3D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelElastic3D::readModel: Failed to allocate memory.");
    this->setRealized(true);

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
void ModelElastic3D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelElastic3D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    Fvp->output(Vpfile.c_str());

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
    Fvp->setType(REGULAR);
    Fvp->setData_format(sizeof(T));
    Fvp->writeHeader();
    Vp = this->getVp();
    Fvp->write(Vp, nx*ny*nz, 0);
    Fvp->close();
}

template<typename T>
void ModelElastic3D<T>::writeVs() {
    if(!this->getRealized()) {
        rs_error("ModelElastic3D::writeVs: Model is not allocated.");
    }
    // Get file names
    std::string Vsfile = this->getVsfile();
    // Open files for writting
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    Fvs->output(Vsfile.c_str());

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

    Fvs->setN(1,nx);
    Fvs->setN(2,ny);
    Fvs->setN(3,nz);
    Fvs->setD(1,dx);
    Fvs->setD(2,dy);
    Fvs->setD(3,dz);
    Fvs->setO(1,ox);
    Fvs->setO(2,oy);
    Fvs->setO(3,oz);
    Fvs->setType(REGULAR);
    Fvs->setData_format(sizeof(T));
    Fvs->writeHeader();
    Vs = this->getVs();
    Fvs->write(Vs, nx*ny*nz, 0);
    Fvs->close();
}

template<typename T>
void ModelElastic3D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelElastic3D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
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

    Frho->setN(1,nx);
    Frho->setN(2,ny);
    Frho->setN(3,nz);
    Frho->setD(1,dx);
    Frho->setD(2,dy);
    Frho->setD(3,dz);
    Frho->setO(1,ox);
    Frho->setO(2,oy);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    R = this->getR();
    Frho->write(R, nx*ny*nz, 0);
    Frho->close();
}

template<typename T>
void ModelElastic3D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelElastic3D::staggerModels: Model is not allocated.");
    }
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

    // Reallocate necessary variables 
    free(L); free(L2M); free(M_xz); free(M_yz); 
    free(M_xy); free(Rx); free(Ry); free(Rz);
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_yz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    
    // Padding
    this->padmodel3d(Rx, R, nx, ny, nz, lpml);
    this->padmodel3d(Ry, R, nx, ny, nz, lpml);
    this->padmodel3d(Rz, R, nx, ny, nz, lpml);
    this->padmodel3d(L, Vp, nx, ny, nz, lpml);
    this->padmodel3d(M_xz, Vs, nx, ny, nz, lpml);
    this->padmodel3d(M_yz, Vs, nx, ny, nz, lpml);
    this->padmodel3d(M_xy, Vs, nx, ny, nz, lpml);
    
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
                L2M[ind_pml(ix,iy,lpml)] *= 0.0;
            }
        }
    }
}

template<typename T>
void ModelElastic3D<T>::createModel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Vp and R */
    free(Vp); free(Vs); free(R);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelElastic3D::createModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelElastic3D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelElastic3D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<rockseis::ModelElastic3D<T>> ModelElastic3D<T>::getLocal(std::shared_ptr<rockseis::Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<rockseis::ModelElastic3D<T>> local;
    /* Get source or receiver min and max positions */
    Point3D<T> *coords;
    size_t ntr = data->getNtrace();
    T min_x, max_x; 
    T min_y, max_y; 
    if(map == SMAP){
        coords = (data->getGeom())->getScoords();
        min_x = coords[0].x;
        max_x = coords[0].x;
        min_y = coords[0].y;
        max_y = coords[0].y;
        for (int i=1; i < ntr; i++){
            if(coords[i].x < min_x) min_x = coords[i].x;
            if(coords[i].x > max_x) max_x = coords[i].x;
            if(coords[i].y < min_y) min_y = coords[i].y;
            if(coords[i].y > max_y) max_y = coords[i].y;
        }
    }else{
        coords = (data->getGeom())->getGcoords();
        min_x = coords[0].x;
        max_x = coords[0].x;
        min_y = coords[0].y;
        max_y = coords[0].y;
        for (size_t i=1; i < ntr; i++){
            if(coords[i].x < min_x) min_x = coords[i].x;
            if(coords[i].x > max_x) max_x = coords[i].x;
            if(coords[i].y < min_y) min_y = coords[i].y;
            if(coords[i].y > max_y) max_y = coords[i].y;
        }
    }

    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();

	/* Determine grid positions and sizes */
    size_t size_x = (size_t) (rintf((max_x-min_x + aperture_x)/dx) + 1);
    if( size_x % 2 == 0 ) size_x++; // Get odd size due to symmetry
    off_t start_x = (off_t) (rintf((min_x - ox)/dx) - (size_x - 1)/2); 

    size_t size_y = (size_t) (rintf((max_y-min_y + aperture_y)/dy) + 1);
    if( size_y % 2 == 0 ) size_y++; // Get odd size due to symmetry
    off_t start_y = (off_t) (rintf((min_y - oy)/dy) - (size_y - 1)/2); 

    /* Create local model */
    local = std::make_shared<rockseis::ModelElastic3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), (ox + start_x*dx), (oy + start_y*dy), this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic3D::getlocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic3D::getlocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic3D::getlocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<rockseis::File> Fvp (new rockseis::File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<rockseis::File> Fvs (new rockseis::File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::getLocal : Error reading from Vs file.");
    }
    std::shared_ptr<rockseis::File> Frho (new rockseis::File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    rockseis::Index l3d(size_x, size_y, nz);
    rockseis::Index f3d(nx, ny, nz);
    rockseis::Index l2d(nx, ny);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f3d(0, 0, i1)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelElastic3D::getLocal: Error reading from vp file");
        Fvs->read(vstrace, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic3D::getLocal: Error reading from vs file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelElastic3D::getLocal: Error reading from rho file");
        for(size_t i3=0; i3<size_y; i3++) {
            lpos_y = j + i3;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<size_x; i2++) {
                lpos_x = i + i2;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
                Vs[l3d(i2,i3,i1)] = vstrace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
            }
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(rhotrace);

    return local;
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
template class Model<float>;
template class Model<double>;

template class ModelAcoustic2D<float>;
template class ModelAcoustic3D<float>;
template class ModelElastic2D<float>;
template class ModelElastic3D<float>;

template class ModelAcoustic2D<double>;
template class ModelAcoustic3D<double>;
template class ModelElastic2D<double>;
template class ModelElastic3D<double>;


}
