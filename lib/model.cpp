// Include statements
#include "model.h"
#include "balloc.h"

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

    domain = std::make_shared<Domain<T>>(); 
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

    domain = std::make_shared<Domain<T>>(); 
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

    domain = std::make_shared<Domain<T>>(); 
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

    //Front and back
    for(ix=0; ix<nx_pml; ix++){
        for(iy=0; iy<lpml; iy++){
            for(iz=0; iz<nz_pml; iz++){
                padded[ind_pml(ix,iy,iz)]=padded[ind_pml(ix,lpml,iz)];
                padded[ind_pml(ix,ny_pml-lpml+iy,iz)]=padded[ind_pml(ix,ny_pml-lpml-1,iz)];
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
                avg = (model[ind(ix+1,iy,iz)] + model[ind(ix,iy,iz)])/2.0;
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
                avg = (model[ind(ix,iy+1,iz)] + model[ind(ix,iy,iz)])/2.0;
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
                avg = (model[ind(ix,iy,iz+1)] + model[ind(ix,iy,iz)])/2.0;
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

template<typename T>
void Model<T>::getLocalsize2d(std::shared_ptr<Data2D<T>> data, T aperture, bool map, off_t *start, size_t *size) {
    /* Get source or receiver min and max positions */
    Point2D<T> *scoords;
    Point2D<T> *gcoords;
    size_t ntr = data->getNtrace();
    double min, max; 
    double off;
    double sx, gx;
    double daperture = aperture;
    T dx = this->getDx();
    T ox = this->getOx();

    /* Determine grid positions and sizes */
    if(aperture >= 0){
        if(map == SMAP){
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sx = scoords[0].x;
            gx = gcoords[0].x;
            min = sx;
            max = sx;
            off = fabs(gx - sx);
            for (long long i=1; i < ntr; i++){
                sx = scoords[i].x;
                gx = gcoords[i].x;
                if(sx < min) min = sx;
                if(sx > max) max = sx;
                if(fabs(gx - sx) > off) off = fabs(gx - sx);
            }
        }else{
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sx = scoords[0].x;
            gx = gcoords[0].x;
            min = gx;
            max = gx;
            off = fabs(gx - sx);
            for (long long i=1; i < ntr; i++){
                sx = scoords[i].x;
                gx = gcoords[i].x;
                if(gx < min) min = gx;
                if(gx > max) max = gx;
                if(fabs(gx - sx) > off) off = fabs(gx - sx);
            }
        }
        if(aperture > 0){
            *size = (size_t) (rintf((max-min + aperture)/dx) + 1);
            if( *size % 2 == 0 ) *size = *size + 1; // Get odd size due to symmetry
        }else{
            *size = (size_t) (2*rintf(off/dx) + 1);
        }
        *start = (off_t) (rintf((min - ox)/dx) - (*size - 1)/2); 
    }else{
        scoords = (data->getGeom())->getScoords();
        gcoords = (data->getGeom())->getGcoords();
        sx = scoords[0].x;
        gx = gcoords[0].x;
        min = sx;
        max = sx;
        for (long long i=0; i < ntr; i++){
            sx = scoords[i].x;
            gx = gcoords[i].x;
            if(sx < min) min = sx;
            if(sx > max) max = sx;
            if(gx < min) min = gx;
            if(gx > max) max = gx;
        }
        *size = (size_t) (rintf((max-min + 2*fabs(daperture))/dx) + 2);
        *start = (off_t) (rintf((min - ox)/dx) - rintf(fabs(daperture/dx))) - 1; 
    }
}

template<typename T>
void Model<T>::getLocalsize3d(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, off_t *start_x, size_t *size_x, off_t *start_y, size_t *size_y)
{
    /* Get source or receiver min and max positions */
    Point3D<T> *scoords;
    Point3D<T> *gcoords;
    size_t ntr = data->getNtrace();
    double sx, gx, sy, gy;
    double min_x, max_x; 
    double min_y, max_y; 
    double off_x, off_y;
    double daperture_x = aperture_x;
    double daperture_y = aperture_y;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();

    /* Determine grid positions and sizes */
    if(aperture_x >= 0){
        if(map == SMAP){
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sx = scoords[0].x;
            gx = gcoords[0].x;
            min_x = sx;
            max_x = sx;
            off_x = fabs(gx - sx);
            for (int i=1; i < ntr; i++){
                sx = scoords[i].x;
                gx = gcoords[i].x;
                if(sx < min_x) min_x = sx;
                if(sx > max_x) max_x = sx;
                if(fabs(gx - sx) > off_x) off_x = fabs(gx - sx);
            }
        }else{
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sx = scoords[0].x;
            gx = gcoords[0].x;
            min_x = gx;
            max_x = gx;
            off_x = fabs(gx - sx);
            for (size_t i=1; i < ntr; i++){
                sx = scoords[i].x;
                gx = gcoords[i].x;
                if(gx < min_x) min_x = gx;
                if(gx > max_x) max_x = gx;
                if(fabs(gx - sx) > off_x) off_x = fabs(gx - sx);
            }
        }
        if(aperture_x > 0){
            *size_x = (size_t) (rintf((max_x-min_x + aperture_x)/dx) + 1);
            if( *size_x % 2 == 0 ) *size_x = *size_x + 1; // Get odd size due to symmetry
        }else{
            *size_x = (size_t) (2*rintf(off_x/dx) + 1);
        }
        *start_x = (off_t) (rintf((min_x - ox)/dx) - (*size_x - 1)/2); 
    }else{
        scoords = (data->getGeom())->getScoords();
        gcoords = (data->getGeom())->getGcoords();
        sx = scoords[0].x;
        gx = gcoords[0].x;
        min_x = sx;
        max_x = sx;
        for (int i=0; i < ntr; i++){
            sx = scoords[i].x;
            gx = gcoords[i].x;
            if(scoords[i].x < min_x) min_x = scoords[i].x;
            if(scoords[i].x > max_x) max_x = scoords[i].x;
            if(gcoords[i].x < min_x) min_x = gcoords[i].x;
            if(gcoords[i].x > max_x) max_x = gcoords[i].x;
        }
        *size_x = (size_t) (rintf((max_x-min_x + 2*fabs(daperture_x))/dx) + 2);
        *start_x = (off_t) (rintf((min_x - ox)/dx) - rintf(fabs(daperture_x/dx))) - 1; 
    }
    if(aperture_y >= 0){
        if(map == SMAP){
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sy = scoords[0].y;
            gy = gcoords[0].y;
            min_y = sy;
            max_y = sy;
            off_y = fabs(gy - sy);
            for (int i=1; i < ntr; i++){
                sy = scoords[i].y;
                gy = gcoords[i].y;
                if(sy < min_y) min_y = sy;
                if(sy > max_y) max_y = sy;
                if(fabs(gy - sy) > off_y) off_y = fabs(gy - sy);
            }
        }else{
            scoords = (data->getGeom())->getScoords();
            gcoords = (data->getGeom())->getGcoords();
            sy = scoords[0].y;
            gy = gcoords[0].y;
            min_y = gy;
            max_y = gy;
            off_y = fabs(gy - sy);
            for (size_t i=1; i < ntr; i++){
                sy = scoords[i].y;
                gy = gcoords[i].y;
                if(gy < min_y) min_y = gy;
                if(gy > max_y) max_y = gy;
                if(fabs(gy - sy) > off_y) off_y = fabs(gy - sy);
            }
        }
        if(aperture_y > 0){
            *size_y = (size_t) (rintf((max_y-min_y + aperture_y)/dy) + 1);
            if( *size_y % 2 == 0 ) *size_y = *size_y + 1; // Get odd size due to symmetry
        }else{
            *size_y = (size_t) (2*rintf(off_y/dy) + 1);
        }
        *start_y = (off_t) (rintf((min_y - oy)/dy) - (*size_y - 1)/2); 
    }else{
        scoords = (data->getGeom())->getScoords();
        gcoords = (data->getGeom())->getGcoords();
            sy = scoords[0].y;
            gy = gcoords[0].y;
            min_y = sy;
            max_y = sy;
        for (int i=0; i < ntr; i++){
            sy = scoords[i].y;
            gy = gcoords[i].y;
            if(sy < min_y) min_y = sy;
            if(sy > max_y) max_y = sy;
            if(gy < min_y) min_y = gy;
            if(gy > max_y) max_y = gy;
        }
        *size_y = (size_t) (rintf((max_y-min_y + 2*fabs(daperture_y))/dy) + 2);
        *start_y = (off_t) (rintf((min_y - oy)/dy) - rintf(fabs(daperture_y/dy))) - 1; 
    }
}

// =============== 2D EIKONAL MODEL CLASS =============== //
template<typename T>
ModelEikonal2D<T>::ModelEikonal2D(): Model<T>(2) {
    /* Allocate variables */
    Velocity = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
}

template<typename T>
ModelEikonal2D<T>::ModelEikonal2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 0.0, _oz, false) {
    /* Allocate variables */
    Velocity = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
}

template<typename T>
ModelEikonal2D<T>::ModelEikonal2D(std::string _Velocityfile, const int _lpml): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Velocityfile = _Velocityfile;

    std::shared_ptr<File> Fmod (new File());
    status = Fmod->input(Velocityfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelEikonal2D::Error reading from Velocity file: ", Velocityfile);
    }

    // Read geometry from file
    nx = Fmod->getN(1);
    dx = (T) Fmod->getD(1);
    ox = (T) Fmod->getO(1);
    nz = Fmod->getN(3);
    dz = (T) Fmod->getD(3);
    oz = (T) Fmod->getO(3);
    
    // Close files
    Fmod->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(false);

    /* Allocate variables */
    Velocity = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
}

template<typename T>
ModelEikonal2D<T>::~ModelEikonal2D() {
    // Freeing all variables
    free(Velocity);
    free(L);
}

template<typename T>
void ModelEikonal2D<T>::readVelocity() {
    bool status;
    // Get file names
    std::string Velocityfile = this->getVelocityfile();
    // Open files for reading
    std::shared_ptr<File> Fmod (new File());
    status = Fmod->input(Velocityfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelEikonal2D::readVelocity : Error reading from Velocity file: ", Velocityfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();

    // Reallocate variables to correct size
    free(Velocity); 
    Velocity = (T *) calloc(nx*nz,sizeof(T));
    if(Velocity == NULL) rs_error("ModelEikonal2D::readVelocity: Failed to allocate memory.");
    this->setRealized(true);

    Velocity = this->getVelocity();
    Fmod->read(Velocity, nx*nz);
    Fmod->close();
}

template<typename T>
void ModelEikonal2D<T>::writeVelocity() {
    if(!this->getRealized()) {
        rs_error("ModelEikonal2D::writeVelocity: Model is not allocated.");
    }
    // Get file names
    std::string Velocityfile = this->getVelocityfile();
    // Open files for writting
    std::shared_ptr<File> Fmod (new File());
    Fmod->output(Velocityfile);
    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fmod->setN(1,nx);
    Fmod->setN(3,nz);
    Fmod->setD(1,dx);
    Fmod->setD(3,dz);
    Fmod->setO(1,ox);
    Fmod->setO(3,oz);
    Fmod->setType(REGULAR);
    Fmod->setData_format(sizeof(T));
    Fmod->writeHeader();
    Velocity = this->getVelocity();
    Fmod->write(Velocity, nx*nz, 0);
    Fmod->close();
}

template<typename T>
void ModelEikonal2D<T>::Expand(){
    if(!this->getRealized()) {
        rs_error("ModelEikonal2D::Expand: Model is not allocated.");
    }
    int nx, nz, lpml, nx_pml, nz_pml;
    nx = this->getNx();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = this->getNx_pml();
    nz_pml = this->getNz_pml();
    
    // Reallocate necessary variables 
    free(L); 
    L = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelEikonal2D::Expand: Failed to allocate memory.");
    
    // Padding
    this->padmodel2d(L, Velocity, nx, nz, lpml);
}

template<typename T>
void ModelEikonal2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate Model and R */
    free(Velocity); 
    Velocity = (T *) calloc(nx*nz,sizeof(T));
    if(Velocity == NULL) rs_error("ModelEikonal2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<ModelEikonal2D<T>> ModelEikonal2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {
    std::shared_ptr<ModelEikonal2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelEikonal2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Velocity = local->getVelocity();

    /* Allocate two traces to read models from file */
    T *veltrace = (T *) calloc(nx, sizeof(T));
    if(veltrace == NULL) rs_error("ModelEikonal2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fmod (new File());
    status = Fmod->input(Velocityfile);
    if(status == FILE_ERR){
	    rs_error("ModelEikonal2D::getLocal : Error reading from Velocity file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        Fmod->read(veltrace, nx, fpos);
        if(Fmod->getFail()) rs_error("ModelEikonal2D::getLocal: Error reading from model file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Velocity[l2d(i2,i1)] = veltrace[lpos];
        }
    }

    /* Free traces */
    free(veltrace);

    return local;
}

// =============== 3D EIKONAL MODEL CLASS =============== //
template<typename T>
ModelEikonal3D<T>::ModelEikonal3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, false) {
    
    /* Allocate minimally the variables */
    Velocity = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
}

template<typename T>
ModelEikonal3D<T>::ModelEikonal3D(std::string _Velocityfile, const int _lpml): Model<T>(3) {
    bool status;
    int nx, ny, nz;
    T dx, dy, dz;
    T ox, oy, oz;
    Velocityfile = _Velocityfile;

    std::shared_ptr<File> Fmodel (new File());
    status = Fmodel->input(Velocityfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelEikonal3D::Error reading from Velocity file: ", Velocityfile);
    }

    // Read geometry from file
    nx = Fmodel->getN(1);
    dx = (T) Fmodel->getD(1);
    ox = (T) Fmodel->getO(1);
    ny = Fmodel->getN(2);
    dy = (T) Fmodel->getD(2);
    oy = (T) Fmodel->getO(2);
    nz = Fmodel->getN(3);
    dz = (T) Fmodel->getD(3);
    oz = (T) Fmodel->getO(3);
    
    // Close files
    Fmodel->close();

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
    this->setFs(false);

    /* Allocate variables */
    Velocity = (T *) calloc(1,1);
    L = (T *) calloc(1,1);
}


template<typename T>
void ModelEikonal3D<T>::readVelocity() {
    bool status;
    // Get file names
    std::string Velocityfile = this->getVelocityfile();
    // Open files for reading
    std::shared_ptr<File> Fmodel (new File());
    status = Fmodel->input(Velocityfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelEikonal3D::readVelocity : Error reading from Velocity file: ", Velocityfile);
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    // Allocate models before reading
    free(Velocity); 
    Velocity = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Velocity == NULL) rs_error("ModelEikonal3D::readVelocity: Failed to allocate memory.");
    this->setRealized(true);

    Velocity = this->getVelocity();
    Fmodel->read(Velocity, nx*ny*nz);
    Fmodel->close();
}

template<typename T>
void ModelEikonal3D<T>::writeVelocity() {
    if(!this->getRealized()) {
        rs_error("ModelEikonal3D::writeVelocity: Model is not allocated.");
    }
    // Get file names
    std::string Velocityfile = this->getVelocityfile();
    // Open files for writting
    std::shared_ptr<File> Fmodel (new File());
    Fmodel->output(Velocityfile.c_str());

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

    Fmodel->setN(1,nx);
    Fmodel->setN(2,ny);
    Fmodel->setN(3,nz);
    Fmodel->setD(1,dx);
    Fmodel->setD(2,dy);
    Fmodel->setD(3,dz);
    Fmodel->setO(1,ox);
    Fmodel->setO(2,oy);
    Fmodel->setO(3,oz);
    Fmodel->setType(REGULAR);
    Fmodel->setData_format(sizeof(T));
    Fmodel->writeHeader();
    Velocity = this->getVelocity();
    Fmodel->write(Velocity, nx*ny*nz, 0);
    Fmodel->close();
}

template<typename T>
ModelEikonal3D<T>::~ModelEikonal3D() {
    // Freeing all variables
    free(Velocity);
    free(L);
}

template<typename T>
void ModelEikonal3D<T>::Expand(){
    if(!this->getRealized()) {
        rs_error("ModelEikonal3D::Expand: Model is not allocated.");
    }
    int nx, ny, nz, lpml, nx_pml, ny_pml, nz_pml;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    lpml = this->getLpml();
    
    nx_pml = nx + 2*lpml;
    ny_pml = ny + 2*lpml;
    nz_pml = nz + 2*lpml;
    
    // Reallocate necessary variables 
    free(L); 
    L = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelEikonal3D::Expand: Failed to allocate memory.");
    
    // Padding
    this->padmodel3d(L, Velocity, nx, ny, nz, lpml);
}


template<typename T>
void ModelEikonal3D<T>::createModel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Model and R */
    free(Velocity);
    Velocity = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Velocity == NULL) rs_error("ModelEikonal3D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
std::shared_ptr<ModelEikonal3D<T>> ModelEikonal3D<T>::getLocal(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<ModelEikonal3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<ModelEikonal3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Velocity = local->getVelocity();

    /* Allocate two traces to read models from file */
    T *veltrace = (T *) calloc(nx*ny, sizeof(T));
    if(veltrace == NULL) rs_error("ModelEikonal3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fmodel (new File());
    status = Fmodel->input(Velocityfile);
    if(status == FILE_ERR){
	    rs_error("ModelEikonal3D::getLocal : Error reading from Velocity file.");
    }
    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    Index l3d(size_x, size_y, nz);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f3d(0, 0, i1)*sizeof(T);
        Fmodel->read(veltrace, nx*ny, fpos);
        if(Fmodel->getFail()) rs_error("ModelEikonal3D::getLocal: Error reading from model file");
        for(size_t i3=0; i3<size_y; i3++) {
            lpos_y = j + i3;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<size_x; i2++) {
                lpos_x = i + i2;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Velocity[l3d(i2,i3,i1)] = veltrace[l2d(lpos_x, lpos_y)];
            }
        }
    }

    /* Free traces */
    free(veltrace);

    return local;
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

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic1D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
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
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic1D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
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

    Fvp->read(Vp, nz);
    Fvp->close();

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
    std::shared_ptr<File> Fvp (new File());
    Fvp->output(Vpfile);
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nz, 0);
    Fvp->close();

    Frho->setN(3,nz);
    Frho->setD(3,dz);
    Frho->setO(3,oz);
    Frho->setType(REGULAR);
    Frho->setData_format(sizeof(T));
    Frho->writeHeader();
    Mod = this->getR();
    Frho->write(Mod, nz, 0);
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

    for(iz=0; iz < nz_pml; iz++){
        if(R[iz] == 0.0) rs_error("staggerModels: Zero density found.");
        Rz[iz] = 1.0/R[iz];
    }

    // In case of free surface
    if(this->getFs()){
            L[lpml] *= 0.0;
    }
}


// =============== 2D ACOUSTIC MODEL CLASS =============== //
template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(): Model<T>(2) {
    /* Allocate variables */
    Vp = (T *) BallocNew(1,1);
    R = (T *) BallocNew(1,1);
    L = (T *) BallocNew(1,1);
    Rx = (T *) BallocNew(1,1);
    Rz = (T *) BallocNew(1,1);
    
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 0.0, _oz, _fs) {
    /* Allocate variables */
    Vp = (T *) BallocNew(1,1);
    R = (T *) BallocNew(1,1);
    L = (T *) BallocNew(1,1);
    Rx = (T *) BallocNew(1,1);
    Rz = (T *) BallocNew(1,1);
}

template<typename T>
ModelAcoustic2D<T>::ModelAcoustic2D(std::string _Vpfile, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Rfile = _Rfile;

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
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
    Vp = (T *) BallocNew(1,1);
    R = (T *) BallocNew(1,1);
    L = (T *) BallocNew(1,1);
    Rx = (T *) BallocNew(1,1);
    Rz = (T *) BallocNew(1,1);
}

template<typename T>
ModelAcoustic2D<T>::~ModelAcoustic2D() {
    // Freeing all variables
    BallocDelete(Vp);
    BallocDelete(R);
    BallocDelete(L);
    BallocDelete(Rx);
    BallocDelete(Rz);
}

template<typename T>
void ModelAcoustic2D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();

    // Reallocate variables to correct size
    BallocDelete(Vp); BallocDelete(R);
    Vp = (T *) BallocNew(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic2D::readModel: Failed to allocate memory.");
    R = (T *) BallocNew(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Fvp->read(Vp, nx*nz);
    Fvp->close();

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
    std::shared_ptr<File> Fvp (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nx*nz, 0);
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
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*nz, 0);
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
    BallocDelete(L); BallocDelete(Rx); BallocDelete(Rz);
    L = (T *) BallocNew(nx_pml*nz_pml,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic2D::staggerModels: Failed to allocate memory.");
    Rx = (T *) BallocNew(nx_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelAcoustic2D::staggerModels: Failed to allocate memory.");
    Rz = (T *) BallocNew(nx_pml*nz_pml,sizeof(T));
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

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            if(Rx[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            if(Rz[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            Rx[ind_pml(ix,iz)] = 1.0/Rx[ind_pml(ix,iz)];
            Rz[ind_pml(ix,iz)] = 1.0/Rz[ind_pml(ix,iz)];
        }
    }
    
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
    BallocDelete(Vp); BallocDelete(R);
    Vp = (T *) BallocNew(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelAcoustic2D::createModel: Failed to allocate memory.");
    R = (T *) BallocNew(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelAcoustic2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
void ModelAcoustic2D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate L Rx, and Rz */
    BallocDelete(L); BallocDelete(Rx); BallocDelete(Rz);
    L = (T *) BallocNew(nx*nz,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic2D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) BallocNew(nx*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelAcoustic2D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) BallocNew(nx*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelAcoustic2D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelAcoustic2D<T>> ModelAcoustic2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {
    std::shared_ptr<ModelAcoustic2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelAcoustic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getLocal : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
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

template<typename T>
std::shared_ptr<ModelAcoustic2D<T>> ModelAcoustic2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd, const int order) {
    std::shared_ptr<ModelAcoustic2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain1D(size+2*lpml,1,nz+2*this->getLpml(),d,nd,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelAcoustic2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getFs());
    (local->getDomain())->setupDomain1D(size+2*lpml,1,nz+2*lpml,d,nd,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();
    T *L = local->getL();
    T *Rx = local->getRx();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from vp file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from rho file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Frho->read(rhotrace_adv, nx, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from rho file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            L[l2d(i2,i1)] = rhotrace[lpos]*vptrace[lpos]*vptrace[lpos];
            if(rhotrace[lpos] <= 0.0) rs_error("ModelAcoustic2D::getDomainmodel: Zero density found.");
            if(lpos < nx-1){
               Rx[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace[lpos+1]);
            }else{
               Rx[l2d(i2,i1)] = 1.0/(rhotrace[lpos]);
            }
            Rz[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace_adv[lpos]);
        }
    }

       
    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rx[l2d(ix,lpml-iz0)] *= 2.0;
            L[l2d(ix,lpml-iz0)] *= 0.0;
        }
    }


    /* Free traces */
    free(vptrace);
    free(rhotrace);
    free(rhotrace_adv);

    
    return local;
}

template<typename T>
std::shared_ptr<ModelAcoustic2D<T>> ModelAcoustic2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd0, const int nd1, const int order) {
    std::shared_ptr<ModelAcoustic2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*this->getLpml(),d,nd0,1,nd1,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelAcoustic2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getFs());
    (local->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*lpml,d,nd0,1,nd1,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();
    T *L = local->getL();
    T *Rx = local->getRx();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelAcoustic2d::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic2D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from vp file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from rho file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Frho->read(rhotrace_adv, nx, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic2D::getDomainmodel: Error reading from rho file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            L[l2d(i2,i1)] = rhotrace[lpos]*vptrace[lpos]*vptrace[lpos];
            if(rhotrace[lpos] <= 0.0) rs_error("ModelAcoustic2D::getDomainmodel: Zero density found.");
            if(lpos < nx-1){
               Rx[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace[lpos+1]);
            }else{
               Rx[l2d(i2,i1)] = 1.0/(rhotrace[lpos]);
            }
            Rz[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace_adv[lpos]);
        }
    }

       
    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rx[l2d(ix,lpml-iz0)] *= 2.0;
            L[l2d(ix,lpml-iz0)] *= 0.0;
        }
    }


    /* Free traces */
    free(vptrace);
    free(rhotrace);
    free(rhotrace_adv);

    
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

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
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
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Frho (new File());
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

    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

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
    std::shared_ptr<File> Fvp (new File());
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
    T * Mod = this->getVp();
    Fvp->write(Mod, nx*ny*nz, 0);
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
    std::shared_ptr<File> Frho (new File());
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
    T * Mod = this->getR();
    Frho->write(Mod, nx*ny*nz, 0);
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

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                if(Rx[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Ry[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Rz[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                Rx[ind_pml(ix,iy,iz)] = 1.0/Rx[ind_pml(ix,iy,iz)];
                Ry[ind_pml(ix,iy,iz)] = 1.0/Ry[ind_pml(ix,iy,iz)];
                Rz[ind_pml(ix,iy,iz)] = 1.0/Rz[ind_pml(ix,iy,iz)];
            }
        }
    }

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
void ModelAcoustic3D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate L Rx, and Rz */
    free(L); free(Rx); free(Ry); free(Rz);
    L = (T *) calloc(nx*ny*nz,sizeof(T));
    if(L == NULL) rs_error("ModelAcoustic3D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelAcoustic3D::createPaddedmodel: Failed to allocate memory.");
    Ry = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ry == NULL) rs_error("ModelAcoustic3D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelAcoustic3D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelAcoustic3D<T>> ModelAcoustic3D<T>::getDomainmodel(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, const int d, const int nd, const int order) {

    std::shared_ptr<ModelAcoustic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();
    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;
    int nxd,nyd,nzd;
    int ix0,iy0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);
    (this->getDomain())->setupDomain1D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd,order);

    nxd = (this->getDomain())->getNx_pad();
    nyd = (this->getDomain())->getNy_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iy0 = (this->getDomain())->getIy0();
    iz0 = (this->getDomain())->getIz0();

    T oxl, oyl, ozl; 
    oxl = (ox + (start_x+ix0-lpml)*dx);
    oyl = (oy + (start_y+iy0-lpml)*dy);
    ozl = (oz + (iz0-lpml)*dz);

    /* Create local model */
    local = std::make_shared<ModelAcoustic3D<T>>(nxd, nyd, nzd, lpml, dx, dy, dz, oxl, oyl, ozl, this->getFs());
    (local->getDomain())->setupDomain1D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();
    T *L = local->getL();
    T *Rx = local->getRx();
    T *Ry = local->getRy();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, lpos_z, fpos;
    Index l3d(nxd, nyd, nzd);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nzd; i1++) {
        lpos_z = iz0 + i1 - lpml;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from vp file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from rho file");

        lpos_z = iz0 + i1 - lpml + 1;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);

        Frho->read(rhotrace_adv, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from rho file");

        for(size_t i3=0; i3<nyd; i3++) {
            lpos_y = j + i3 + iy0 - lpml;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<nxd; i2++) {
                lpos_x = i + i2 + ix0 - lpml;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
                L[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)];

                if(rhotrace[l2d(lpos_x, lpos_y)] <= 0.0) rs_error("ModelAcoustic3D::getDomainmodel: Zero density found.");

                if(lpos_x < nx-1){
                    Rx[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x+1, lpos_y)]);
                }else{
                    Rx[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                if(lpos_y < ny-1){
                    Ry[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x, lpos_y+1)]);
                }else{
                    Ry[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                Rz[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace_adv[l2d(lpos_x, lpos_y)]);
            }
        }
    }

    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            for(size_t iy=0; iy<nyd; iy++){
                Rx[l3d(ix,iy,lpml-iz0)] *= 2.0;
                Ry[l3d(ix,iy,lpml-iz0)] *= 2.0;
                L[l3d(ix,iy,lpml-iz0)] *= 0.0;
            }
        }
    }


    /* Free traces */
    free(vptrace);
    free(rhotrace);
    free(rhotrace_adv);

    return local;
}

template<typename T>
std::shared_ptr<ModelAcoustic3D<T>> ModelAcoustic3D<T>::getDomainmodel(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order) {

    std::shared_ptr<ModelAcoustic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();
    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;
    int nxd,nyd,nzd;
    int ix0,iy0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);
    (this->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    nxd = (this->getDomain())->getNx_pad();
    nyd = (this->getDomain())->getNy_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iy0 = (this->getDomain())->getIy0();
    iz0 = (this->getDomain())->getIz0();

    T oxl, oyl, ozl; 
    oxl = (ox + (start_x+ix0-lpml)*dx);
    oyl = (oy + (start_y+iy0-lpml)*dy);
    ozl = (oz + (iz0-lpml)*dz);

    /* Create local model */
    local = std::make_shared<ModelAcoustic3D<T>>(nxd, nyd, nzd, lpml, dx, dy, dz, oxl, oyl, ozl, this->getFs());
    (local->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();
    T *L = local->getL();
    T *Rx = local->getRx();
    T *Ry = local->getRy();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, lpos_z, fpos;
    Index l3d(nxd, nyd, nzd);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nzd; i1++) {
        lpos_z = iz0 + i1 - lpml;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from vp file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from rho file");

        lpos_z = iz0 + i1 - lpml + 1;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);

        Frho->read(rhotrace_adv, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelAcoustic3D::getLocal: Error reading from rho file");

        for(size_t i3=0; i3<nyd; i3++) {
            lpos_y = j + i3 + iy0 - lpml;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<nxd; i2++) {
                lpos_x = i + i2 + ix0 - lpml;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
                L[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)];

                if(rhotrace[l2d(lpos_x, lpos_y)] <= 0.0) rs_error("ModelAcoustic3D::getDomainmodel: Zero density found.");

                if(lpos_x < nx-1){
                    Rx[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x+1, lpos_y)]);
                }else{
                    Rx[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                if(lpos_y < ny-1){
                    Ry[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x, lpos_y+1)]);
                }else{
                    Ry[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                Rz[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace_adv[l2d(lpos_x, lpos_y)]);
            }
        }
    }

    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            for(size_t iy=0; iy<nyd; iy++){
                Rx[l3d(ix,iy,lpml-iz0)] *= 2.0;
                Ry[l3d(ix,iy,lpml-iz0)] *= 2.0;
                L[l3d(ix,iy,lpml-iz0)] *= 0.0;
            }
        }
    }


    /* Free traces */
    free(vptrace);
    free(rhotrace);
    free(rhotrace_adv);

    return local;
}


template<typename T>
std::shared_ptr<ModelAcoustic3D<T>> ModelAcoustic3D<T>::getLocal(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<ModelAcoustic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<ModelAcoustic3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelAcoustic3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelAcoustic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    Index l3d(size_x, size_y, nz);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
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

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::Error reading from Vp file: ", Vpfile);
	    exit(1);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic2D::Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
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
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
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

    Fvp->read(Vp, nx*nz);
    Fvp->close();

    Fvs->read(Vs, nx*nz);
    Fvs->close();

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
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*nz, 0);
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
    std::shared_ptr<File> Fvp (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nx*nz, 0);
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
    std::shared_ptr<File> Fvs (new File());
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
    T *Mod = this->getVs();
    Fvs->write(Mod, nx*nz, 0);
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

    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            L[ind_pml(ix,iz)] = 0.0;
            L2M[ind_pml(ix,iz)] = M[ind_pml(ix,iz)];
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, 1, nz_pml); 
    this->staggermodel_z(Rz, nx_pml, 1, nz_pml); 
    
    this->staggermodel_z(M, nx_pml, 1, nz_pml); 
    this->staggermodel_x(M, nx_pml, 1, nz_pml); 

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            if(Rx[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            if(Rz[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            Rx[ind_pml(ix,iz)] = 1.0/Rx[ind_pml(ix,iz)];
            Rz[ind_pml(ix,iz)] = 1.0/Rz[ind_pml(ix,iz)];
        }
    }
    
    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            Rx[ind_pml(ix,iz)] *= 2.0;
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
void ModelElastic2D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate L2M, L , M, Rx, and Rz */
    free(L2M); free(L); free(M); free(Rx); free(Rz);

    L2M = (T *) calloc(nx*nz,sizeof(T));
    if(L2M == NULL) rs_error("ModelElastic2D::createPaddedmodel: Failed to allocate memory.");
    L = (T *) calloc(nx*nz,sizeof(T));
    if(L == NULL) rs_error("ModelElastic2D::createPaddedmodel: Failed to allocate memory.");
    M = (T *) calloc(nx*nz,sizeof(T));
    if(M == NULL) rs_error("ModelElastic2D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelElastic2D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelElastic2D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelElastic2D<T>> ModelElastic2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {

    std::shared_ptr<ModelElastic2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelElastic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic2d::getLocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Vs file.");
    }

    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
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
std::shared_ptr<ModelElastic2D<T>> ModelElastic2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd0, const int nd1, const int order) {
    std::shared_ptr<ModelElastic2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*this->getLpml(),d,nd0,1,nd1,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelElastic2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getFs());
    (local->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*lpml,d,nd0,1,nd1,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();
    T *M = local->getM();
    T *L = local->getL();
    T *L2M = local->getL2M();
    T *Rx = local->getRx();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic2d::getDomainmodel: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic2d::getLocal: Failed to allocate memory.");
    T *vstrace_adv = (T *) calloc(nx, sizeof(T));
    if(vstrace_adv == NULL) rs_error("ModelElastic2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelElastic2d::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getLocal : Error reading from Vs file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic2D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);
    T M1, M2, M3, M4;

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelElastic2D::getDomainmodel: Error reading from vp file");
        Fvs->read(vstrace, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic2D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelElastic2D::getDomainmodel: Error reading from rho file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvs->read(vstrace_adv, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic2D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace_adv, nx, fpos);
        if(Frho->getFail()) rs_error("ModelElastic2D::getDomainmodel: Error reading from rho file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            Vs[l2d(i2,i1)] = vstrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            L2M[l2d(i2,i1)] = rhotrace[lpos]*vptrace[lpos]*vptrace[lpos];
            L[l2d(i2,i1)] = rhotrace[lpos]*vptrace[lpos]*vptrace[lpos] - 2.0*rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];
            if(rhotrace[lpos] <= 0.0) rs_error("ModelElastic2D::getDomainmodel: Zero density found.");
            if(lpos < nx-1){
               Rx[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace[lpos+1]);
            }else{
               Rx[l2d(i2,i1)] = 1.0/(rhotrace[lpos]);
            }
            Rz[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace_adv[lpos]);

            M1 = rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];
            M3 = rhotrace_adv[lpos]*vstrace_adv[lpos]*vstrace_adv[lpos];
            if(lpos < nx-1){
               M2 = rhotrace[lpos+1]*vstrace[lpos+1]*vstrace[lpos+1];
               M4 = rhotrace_adv[lpos+1]*vstrace_adv[lpos+1]*vstrace_adv[lpos+1];
            }else{
               M2 = rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];
               M4 = rhotrace_adv[lpos]*vstrace_adv[lpos]*vstrace_adv[lpos];
            }
            M[l2d(i2,i1)] = 0.25*(M1+M2+M3+M4);
        }
    }

       
    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rx[l2d(ix,lpml-iz0)] *= 2.0;
            L2M[l2d(ix,lpml-iz0)] = (L2M[l2d(ix,lpml-iz0)] - L[l2d(ix,lpml-iz0)])/2.0;
            L[l2d(ix,lpml-iz0)] = 0.0;
        }
    }


    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(vstrace_adv);
    free(rhotrace);
    free(rhotrace_adv);
    
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
ModelElastic3D<T>::ModelElastic3D(): Model<T>(3) {
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

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::Error reading from Vp file:  ", _Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::Error reading from Vs file: ", _Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
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
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
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

    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

    Fvs->read(Vs, nx*ny*nz);
    Fvs->close();

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
    std::shared_ptr<File> Fvp (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nx*ny*nz, 0);
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
    std::shared_ptr<File> Fvs (new File());
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
    T *Mod = this->getVs();
    Fvs->write(Mod, nx*ny*nz, 0);
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
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*ny*nz, 0);
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

    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy < ny_pml; iy++){
                L[ind_pml(ix,iy,iz)] = 0.0;
                L2M[ind_pml(ix,iy,iz)] = M_xz[ind_pml(ix,iy,iz)];
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

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                if(Rx[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Ry[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Rz[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                Rx[ind_pml(ix,iy,iz)] = 1.0/Rx[ind_pml(ix,iy,iz)];
                Ry[ind_pml(ix,iy,iz)] = 1.0/Ry[ind_pml(ix,iy,iz)];
                Rz[ind_pml(ix,iy,iz)] = 1.0/Rz[ind_pml(ix,iy,iz)];
            }
        }
    }

    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy<ny_pml; iy++){
                Rx[ind_pml(ix,iy,lpml)] *= 2.0;
                Ry[ind_pml(ix,iy,lpml)] *= 2.0;
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
void ModelElastic3D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate L Rx, and Rz */
    free(L); free(L2M); free(M_xz); free(M_yz); 
    free(M_xy); free(Rx); free(Ry); free(Rz);
    L2M = (T *) calloc(nx*ny*nz,sizeof(T));
    if(L2M == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    L = (T *) calloc(nx*ny*nz,sizeof(T));
    if(L == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    M_xz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_xz == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    M_yz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_yz == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    M_xy = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_xy == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    Ry = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ry == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelElastic3D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelElastic3D<T>> ModelElastic3D<T>::getLocal(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<ModelElastic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<ModelElastic3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic3D::getLocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
        rs_error("ModelElastic3D::getLocal : Error reading from Vs file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    Index l3d(size_x, size_y, nz);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
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
std::shared_ptr<ModelElastic3D<T>> ModelElastic3D<T>::getDomainmodel(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order) {

    std::shared_ptr<ModelElastic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();
    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;
    int nxd,nyd,nzd;
    int ix0,iy0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);
    (this->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    nxd = (this->getDomain())->getNx_pad();
    nyd = (this->getDomain())->getNy_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iy0 = (this->getDomain())->getIy0();
    iz0 = (this->getDomain())->getIz0();

    T oxl, oyl, ozl; 
    oxl = (ox + (start_x+ix0-lpml)*dx);
    oyl = (oy + (start_y+iy0-lpml)*dy);
    ozl = (oz + (iz0-lpml)*dz);

    /* Create local model */
    local = std::make_shared<ModelElastic3D<T>>(nxd, nyd, nzd, lpml, dx, dy, dz, oxl, oyl, ozl, this->getFs());
    (local->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();
    T *L2M = local->getL2M();
    T *L = local->getL();
    T *M_xz = local->getM_xz();
    T *M_yz = local->getM_yz();
    T *M_xy = local->getM_xy();
    T *Rx = local->getRx();
    T *Ry = local->getRy();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelElastic3D::getDomainmodel: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace == NULL) rs_error("ModelElastic3D::getDomainmodel: Failed to allocate memory.");
    T *vstrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace_adv == NULL) rs_error("ModelElastic3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelElastic3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelElastic3D::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getDomainmodel : Error reading from Density file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelElastic3D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, lpos_z, fpos;
    Index l3d(nxd, nyd, nzd);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nzd; i1++) {
        lpos_z = iz0 + i1 - lpml;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelElastic3D::getDomainmodel: Error reading from vp file");
        Fvs->read(vstrace, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic3D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelElastic3D::getDomainmodel: Error reading from rho file");

        lpos_z = iz0 + i1 - lpml + 1;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);

        Fvs->read(vstrace_adv, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelElastic3D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace_adv, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelElastic3D::getDomainmodel: Error reading from rho file");

        T M1, M2, M3, M4, M5, M6, M7;
        for(size_t i3=0; i3<nyd; i3++) {
            lpos_y = j + i3 + iy0 - lpml;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<nxd; i2++) {
                lpos_x = i + i2 + ix0 - lpml;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
                Vs[l3d(i2,i3,i1)] = vstrace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
                L2M[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)];
                L[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]- 2.0*rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];

                if(rhotrace[l2d(lpos_x, lpos_y)] <= 0.0) rs_error("ModelElastic3D::getDomainmodel: Zero density found.");

                if(lpos_x < nx-1){
                    Rx[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x+1, lpos_y)]);
                }else{
                    Rx[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                if(lpos_y < ny-1){
                    Ry[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x, lpos_y+1)]);
                }else{
                    Ry[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                Rz[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace_adv[l2d(lpos_x, lpos_y)]);
                M1 = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
                M5 = rhotrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)];

                if(lpos_x < nx-1){
                   M2 = rhotrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)];
                   M6 = rhotrace_adv[l2d(lpos_x+1, lpos_y)]*vstrace_adv[l2d(lpos_x+1, lpos_y)]*vstrace_adv[l2d(lpos_x+1, lpos_y)];
                   if(lpos_y < ny-1){
                      M4 = rhotrace[l2d(lpos_x+1, lpos_y+1)]*vstrace[l2d(lpos_x+1, lpos_y+1)]*vstrace[l2d(lpos_x+1, lpos_y+1)];
                   }else{
                      M4 = rhotrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)];
                   }
                }else{
                   M2 = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
                   M6 = rhotrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)];
                   if(lpos_y < ny-1){
                      M4 = rhotrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)];
                   }else{
                      M4 = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
                   }
                }

                if(lpos_y < ny-1){
                   M3 = rhotrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)];
                   M7 = rhotrace_adv[l2d(lpos_x, lpos_y+1)]*vstrace_adv[l2d(lpos_x, lpos_y+1)]*vstrace_adv[l2d(lpos_x, lpos_y+1)];
                }else{
                   M3 = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
                   M7 = rhotrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)];
                }
                M_xz[l3d(i2,i3,i1)] = 0.25*(M1+M2+M5+M6);
                M_yz[l3d(i2,i3,i1)] = 0.25*(M1+M3+M5+M7);
                M_xy[l3d(i2,i3,i1)] = 0.25*(M1+M2+M3+M4);

            }
        }
    }

    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            for(size_t iy=0; iy<nyd; iy++){
                Rx[l3d(ix,iy,lpml-iz0)] *= 2.0;
                Ry[l3d(ix,iy,lpml-iz0)] *= 2.0;
                L2M[l3d(ix,iy,lpml-iz0)] = (L2M[l3d(ix,iy,lpml-iz0)] - L[l3d(ix,iy,lpml-iz0)])/2;
                L[l3d(ix,iy,lpml-iz0)] = 0.0;
            }
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(vstrace_adv);
    free(rhotrace);
    free(rhotrace_adv);

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

// =============== 2D VISCOELASTIC MODEL CLASS =============== //
template<typename T>
ModelViscoelastic2D<T>::ModelViscoelastic2D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelViscoelastic2D<T>::ModelViscoelastic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const T _f0, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

    Qp = (T *) calloc(1,1);
    Qs = (T *) calloc(1,1);
    Tp = (T *) calloc(1,1);
    Ts = (T *) calloc(1,1);
    Ts_xz = (T *) calloc(1,1);
    f0  = _f0;

}

template<typename T>
ModelViscoelastic2D<T>::ModelViscoelastic2D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, std::string _Qpfile, std::string _Qsfile, const int _lpml, const T _f0, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;
    Qpfile = _Qpfile;
    Qsfile = _Qsfile;
    f0 = _f0;

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::Error reading from Vp file: ", Vpfile);
	    exit(1);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from density file: ", Rfile);
        exit(1);
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from Qp file: ", Qpfile);
        exit(1);
    }
    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from Qs file: ", Qsfile);
        exit(1);
    }



    // Compare geometry in the two files
    if(Fvp->compareGeometry(Fvs) != 0)
    {
        rs_error("ModelViscoelastic2D::Geometries in Vp and Vs model files do not match.");
    }

    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelViscoelastic2D::Geometries in Vp and Density model files do not match.");
    }
    if(Fvp->compareGeometry(Fqp) != 0)
    {
        rs_error("ModelViscoelastic2D::Geometries in Vp and Qp model files do not match.");
    }
    if(Fvp->compareGeometry(Fqs) != 0)
    {
        rs_error("ModelViscoelastic2D::Geometries in Vp and Qs model files do not match.");
    }


    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelViscoelastic2D::Numerical precision mismatch in Vp and Density model files.");
    }
    if(Fvp->getData_format() != Fvs->getData_format())
    {
        rs_error("ModelViscoelastic2D::Numerical precision mismatch in Vp and Vs model files.");
    }
    if(Fvp->getData_format() != Fqp->getData_format())
    {
        rs_error("ModelViscoelastic2D::Numerical precision mismatch in Vp and Qp model files.");
    }
    if(Fvp->getData_format() != Fqs->getData_format())
    {
        rs_error("ModelViscoelastic2D::Numerical precision mismatch in Vp and Qs model files.");
    }
    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelViscoelastic2D::Numerical precision in Vp, Vs and Density model files mismatch with constructor.");
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
    Fqp->close();
    Fqs->close();

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
    L2M = (T *) calloc(1,1);
    M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
    Qp = (T *) calloc(1,1);
    Qs = (T *) calloc(1,1);
    Tp = (T *) calloc(1,1);
    Ts = (T *) calloc(1,1);
    Ts_xz = (T *) calloc(1,1);
}

template<typename T>
void ModelViscoelastic2D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    std::string Qpfile = this->getQpfile();
    std::string Qsfile = this->getQsfile();
    // Open files for reading
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::readModel : Error reading from Density file: ", Rfile);
    }

    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::readModel : Error reading from Qp file: ", Qpfile);
    }

    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::readModel : Error reading from Qs file: ", Qsfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();
    
    /* Reallocate Vp, Vs  and R */
    free(Vp); free(Vs); free(R); free(Qp); free(Qs);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    Qp = (T *) calloc(nx*nz,sizeof(T));
    if(Qp == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    Qs = (T *) calloc(nx*nz,sizeof(T));
    if(Qs == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Fvp->read(Vp, nx*nz);
    Fvp->close();

    Fvs->read(Vs, nx*nz);
    Fvs->close();

    Frho->read(R, nx*nz);
    Frho->close();

    Fqp->read(Qp, nx*nz);
    Fqp->close();

    Fqs->read(Qs, nx*nz);
    Fqs->close();
}

template<typename T>
void ModelViscoelastic2D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*nz, 0);
    Frho->close();
}


template<typename T>
void ModelViscoelastic2D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    // Open files for writting
    std::shared_ptr<File> Fvp (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nx*nz, 0);
    Fvp->close();
}


template<typename T>
void ModelViscoelastic2D<T>::writeVs() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::writeVs: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<File> Fvs (new File());
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
    T *Mod = this->getVs();
    Fvs->write(Mod, nx*nz, 0);
    Fvs->close();
}

template<typename T>
void ModelViscoelastic2D<T>::writeQp() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::writeQp: Model is not allocated.");
    }
    // Get file names
    std::string Qpfile = this->getQpfile();
    // Open files for writting
    std::shared_ptr<File> Fqp (new File());
    Fqp->output(Qpfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fqp->setN(1,nx);
    Fqp->setN(3,nz);
    Fqp->setD(1,dx);
    Fqp->setD(3,dz);
    Fqp->setO(1,ox);
    Fqp->setO(3,oz);
    Fqp->setType(REGULAR);
    Fqp->setData_format(sizeof(T));
    Fqp->writeHeader();
    T *Mod = this->getQp();
    Fqp->write(Mod, nx*nz, 0);
    Fqp->close();
}

template<typename T>
void ModelViscoelastic2D<T>::writeQs() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::writeQs: Model is not allocated.");
    }
    // Get file names
    std::string Qsfile = this->getQsfile();
    // Open files for writting
    std::shared_ptr<File> Fqs (new File());
    Fqs->output(Qsfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    Fqs->setN(1,nx);
    Fqs->setN(3,nz);
    Fqs->setD(1,dx);
    Fqs->setD(3,dz);
    Fqs->setO(1,ox);
    Fqs->setO(3,oz);
    Fqs->setType(REGULAR);
    Fqs->setData_format(sizeof(T));
    Fqs->writeHeader();
    T *Mod = this->getQs();
    Fqs->write(Mod, nx*nz, 0);
    Fqs->close();
}



template<typename T>
void ModelViscoelastic2D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic2D::staggerModels: Model is not allocated.");
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
    free(M_xz); free(L2M); free(M); free(Rx); free(Rz); free(Tp); free(Ts); free(Ts_xz);
    M_xz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(M_xz == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    L2M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(L2M == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    M = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(M == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    Tp = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Tp == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    Ts = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Ts == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");
    Ts_xz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Ts_xz == NULL) rs_error("ModelViscoelastic2D::staggerModels: Failed to allocate memory.");

    // Padding
    this->padmodel2d(Rx, R, nx, nz, lpml);
    this->padmodel2d(Rz, R, nx, nz, lpml);
    this->padmodel2d(M_xz, Vp, nx, nz, lpml);
    this->padmodel2d(M, Vs, nx, nz, lpml);
    this->padmodel2d(Tp, Qp, nx, nz, lpml);
    this->padmodel2d(Ts, Qs, nx, nz, lpml);
    this->padmodel2d(Ts_xz, Qs, nx, nz, lpml);
    
    // Computing moduli
    T _rho, _vp, _vs;
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            _vp=M_xz[ind_pml(ix,iz)];
            _vs=M[ind_pml(ix,iz)];
            _rho=Rx[ind_pml(ix,iz)];
            L2M[ind_pml(ix,iz)]=_rho*_vp*_vp;
            M[ind_pml(ix,iz)]=_rho*_vs*_vs;
            M_xz[ind_pml(ix,iz)]=_rho*_vs*_vs;
            Tp[ind_pml(ix,iz)] = (Tp[ind_pml(ix,iz)]+1)/(Tp[ind_pml(ix,iz)]-1);
            Ts[ind_pml(ix,iz)] = (Ts[ind_pml(ix,iz)]+1)/(Ts[ind_pml(ix,iz)]-1);
            Ts_xz[ind_pml(ix,iz)] = (Ts_xz[ind_pml(ix,iz)]+1)/(Ts_xz[ind_pml(ix,iz)]-1);
        }
    }

    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            L2M[ind_pml(ix,iz)] = M[ind_pml(ix,iz)];
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, 1, nz_pml); 
    this->staggermodel_z(Rz, nx_pml, 1, nz_pml); 
    
    this->staggermodel_z(M_xz, nx_pml, 1, nz_pml); 
    this->staggermodel_x(M_xz, nx_pml, 1, nz_pml); 

    this->staggermodel_z(Ts_xz, nx_pml, 1, nz_pml); 
    this->staggermodel_x(Ts_xz, nx_pml, 1, nz_pml); 

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            if(Rx[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            if(Rz[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            Rx[ind_pml(ix,iz)] = 1.0/Rx[ind_pml(ix,iz)];
            Rz[ind_pml(ix,iz)] = 1.0/Rz[ind_pml(ix,iz)];
        }
    }
    
    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            Rx[ind_pml(ix,iz)] *= 2.0;
        }
    }
}

template<typename T>
void ModelViscoelastic2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();


    /* Reallocate Vp, Vs and R */
    free(Vp); free(Vs); free(R); free(Qp); free(Qs);
    Vp = (T *) calloc(nx*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    Qp = (T *) calloc(nx*nz,sizeof(T));
    if(Qp == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    Qs = (T *) calloc(nx*nz,sizeof(T));
    if(Qs == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
void ModelViscoelastic2D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int nz = this->getNz();

    // Reallocate necessary variables 
    free(M_xz); free(L2M); free(M); free(Rx); free(Rz); free(Tp); free(Ts); free(Ts_xz);
    M_xz = (T *) calloc(nx*nz,sizeof(T));
    if(M_xz == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    L2M = (T *) calloc(nx*nz,sizeof(T));
    if(L2M == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    M = (T *) calloc(nx*nz,sizeof(T));
    if(M == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    Tp = (T *) calloc(nx*nz,sizeof(T));
    if(Tp == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    Ts = (T *) calloc(nx*nz,sizeof(T));
    if(Ts == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");
    Ts_xz = (T *) calloc(nx*nz,sizeof(T));
    if(Ts_xz == NULL) rs_error("ModelViscoelastic2D::createPaddedmodel: Failed to allocate memory.");

}



template<typename T>
std::shared_ptr<ModelViscoelastic2D<T>> ModelViscoelastic2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {

    std::shared_ptr<ModelViscoelastic2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelViscoelastic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getF0(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();
    T *Qp = local->getQp();
    T *Qs = local->getQs();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx, sizeof(T));
    if(vstrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *qptrace = (T *) calloc(nx, sizeof(T));
    if(qptrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *qstrace = (T *) calloc(nx, sizeof(T));
    if(qstrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getLocal : Error reading from Vs file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getLocal : Error reading from Density file.");
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::getLocal : Error reading from Qp file: ", Qpfile);
    }

    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::getLocal : Error reading from Qs file: ", Qsfile);
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from vp file");
        Fvs->read(vstrace, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from vs file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from rho file");
        Fqp->read(qptrace, nx, fpos);
        if(Fqp->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qp file");
        Fqs->read(qstrace, nx, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qs file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            Vs[l2d(i2,i1)] = vstrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            Qp[l2d(i2,i1)] = qptrace[lpos];
            Qs[l2d(i2,i1)] = qstrace[lpos];
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(rhotrace);
    free(qptrace);
    free(qstrace);

    return local;
}

template<typename T>
std::shared_ptr<ModelViscoelastic2D<T>> ModelViscoelastic2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd0, const int nd1, const int order) {
    std::shared_ptr<ModelViscoelastic2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*this->getLpml(),d,nd0,1,nd1,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelViscoelastic2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getF0(), this->getFs());
    (local->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*lpml,d,nd0,1,nd1,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *Qp = local->getQp();
    T *Qs = local->getQs();
    T *R = local->getR();
    T *M = local->getM();
    T *M_xz = local->getM_xz();
    T *L2M = local->getL2M();
    T *Rx = local->getRx();
    T *Rz = local->getRz();

    T *Tp = local->getTp();
    T *Ts = local->getTs();
    T *Ts_xz = local->getTs_xz();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx, sizeof(T));
    if(vptrace == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx, sizeof(T));
    if(vstrace == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *vstrace_adv = (T *) calloc(nx, sizeof(T));
    if(vstrace_adv == NULL) rs_error("ModelViscoelastic2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");
    T *qptrace = (T *) calloc(nx, sizeof(T));
    if(qptrace == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");
    T *qstrace = (T *) calloc(nx, sizeof(T));
    if(qstrace == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");
    T *qstrace_adv = (T *) calloc(nx, sizeof(T));
    if(qstrace_adv == NULL) rs_error("ModelViscoelastic2d::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getDomainmodel : Error reading from Vs file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic2D::getDomainmodel : Error reading from Density file.");
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::getDomainmodel : Error reading from Qp file: ", Qpfile);
    }

    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::getDomainmodel : Error reading from Qs file: ", Qsfile);
    }



    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);
    T M1, M2, M3, M4;

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvp->read(vptrace, nx, fpos);
        if(Fvp->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from vp file");
        Fvs->read(vstrace, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from rho file");
        Fqp->read(qptrace, nx, fpos);
        if(Fqp->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qp file");
        Fqs->read(qstrace, nx, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qs file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        Fvs->read(vstrace_adv, nx, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace_adv, nx, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from rho file");
        Fqs->read(qstrace_adv, nx, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from qs file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Vp[l2d(i2,i1)] = vptrace[lpos];
            Vs[l2d(i2,i1)] = vstrace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            Qp[l2d(i2,i1)] = qptrace[lpos];
            Qs[l2d(i2,i1)] = qstrace[lpos];
            L2M[l2d(i2,i1)] = rhotrace[lpos]*vptrace[lpos]*vptrace[lpos];
            M[l2d(i2,i1)] = rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];

            Tp[l2d(i2,i1)] = (qptrace[lpos]+1)/(qptrace[lpos]-1);
            Ts[l2d(i2,i1)] = (qstrace[lpos]+1)/(qstrace[lpos]-1);
            if(rhotrace[lpos] <= 0.0) rs_error("ModelViscoelastic2D::getDomainmodel: Zero density found.");
            if(lpos < nx-1){
               Rx[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace[lpos+1]);
            }else{
               Rx[l2d(i2,i1)] = 1.0/(rhotrace[lpos]);
            }
            Rz[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace_adv[lpos]);

            M1 = rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];
            M3 = rhotrace_adv[lpos]*vstrace_adv[lpos]*vstrace_adv[lpos];
            if(lpos < nx-1){
               M2 = rhotrace[lpos+1]*vstrace[lpos+1]*vstrace[lpos+1];
               M4 = rhotrace_adv[lpos+1]*vstrace_adv[lpos+1]*vstrace_adv[lpos+1];
            }else{
               M2 = rhotrace[lpos]*vstrace[lpos]*vstrace[lpos];
               M4 = rhotrace_adv[lpos]*vstrace_adv[lpos]*vstrace_adv[lpos];
            }
            M_xz[l2d(i2,i1)] = 0.25*(M1+M2+M3+M4);

            M1 = (qstrace[lpos]+1)/(qstrace[lpos]-1);
            M3 = (qstrace_adv[lpos]+1)/(qstrace_adv[lpos]-1);
            if(lpos < nx-1){
                M2 = (qstrace[lpos+1]+1)/(qstrace[lpos+1]-1);
                M4 = (qstrace_adv[lpos+1]+1)/(qstrace_adv[lpos+1]-1);
            }else{
                M2 = M1;
                M4 = M3;
            }
            Ts_xz[l2d(i2,i1)] = 0.25*(M1+M2+M3+M4);
        }
    }

       
    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rx[l2d(ix,lpml-iz0)] *= 2.0;
            L2M[l2d(ix,lpml-iz0)] = M[l2d(ix,lpml-iz0)];
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(vstrace_adv);
    free(rhotrace);
    free(rhotrace_adv);
    free(qptrace);
    free(qstrace);
    free(qstrace_adv);
    
    return local;
}

template<typename T>
ModelViscoelastic2D<T>::~ModelViscoelastic2D() {
    free(Vp);
    free(Vs);
    free(R);
    free(M);
    free(L2M);
    free(M_xz);
    free(Rx);
    free(Rz);
    free(Qp);
    free(Qs);
    free(Tp);
    free(Ts);
    free(Ts_xz);
}

// =============== 3D VISCOELASTIC MODEL CLASS =============== //
template<typename T>
ModelViscoelastic3D<T>::ModelViscoelastic3D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelViscoelastic3D<T>::ModelViscoelastic3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const T _f0, const bool _fs): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, _fs) {
    
    /* Allocate variables */
    Vp = (T *) calloc(1,1);
    Vs = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    M = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    M_yz = (T *) calloc(1,1);
    M_xy = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

    Qp = (T *) calloc(1,1);
    Qs = (T *) calloc(1,1);
    Tp = (T *) calloc(1,1);
    Ts = (T *) calloc(1,1);
    Ts_xz = (T *) calloc(1,1);
    Ts_yz = (T *) calloc(1,1);
    Ts_xy = (T *) calloc(1,1);
    f0 = _f0;
}

template<typename T>
ModelViscoelastic3D<T>::ModelViscoelastic3D(std::string _Vpfile, std::string _Vsfile, std::string _Rfile, std::string _Qpfile, std::string _Qsfile, const T _f0, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz;
    T dx, dy, dz;
    T ox, oy, oz;
    Vpfile = _Vpfile;
    Vsfile = _Vsfile;
    Rfile = _Rfile;
    Qpfile = _Qpfile;
    Qsfile = _Qsfile;
    f0 = _f0;

    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::Error reading from Vp file:  ", _Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::Error reading from Vs file: ", _Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::Error reading from density file: ", _Rfile);
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::Error reading from Qp file: ", Qpfile);
        exit(1);
    }
    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::Error reading from Qs file: ", Qsfile);
        exit(1);
    }

    // Compare geometry in the two files
    if(Fvp->compareGeometry(Fvs) != 0)
    {
        rs_error("ModelViscoelastic3D::Geometries in Vp and Vs model files do not match.");
    }

    if(Fvp->compareGeometry(Frho) != 0)
    {
        rs_error("ModelViscoelastic3D::Geometries in Vp and Density model files do not match.");
    }
    if(Fvp->compareGeometry(Fqp) != 0)
    {
        rs_error("ModelViscoelastic3D::Geometries in Vp and Qp model files do not match.");
    }
    if(Fvp->compareGeometry(Fqs) != 0)
    {
        rs_error("ModelViscoelastic3D::Geometries in Vp and Qs model files do not match.");
    }

    if(Fvp->getData_format() != Frho->getData_format())
    {
        rs_error("ModelViscoelastic3D::Numerical precision mismatch in Vp and Density model files.");
    }
    if(Fvp->getData_format() != Fvs->getData_format())
    {
        rs_error("ModelViscoelastic3D::Numerical precision mismatch in Vp and Vs model files.");
    }
    if(Fvp->getData_format() != Fqp->getData_format())
    {
        rs_error("ModelViscoelastic3D::Numerical precision mismatch in Vp and Qp model files.");
    }
    if(Fvp->getData_format() != Fqs->getData_format())
    {
        rs_error("ModelViscoelastic3D::Numerical precision mismatch in Vp and Qs model files.");
    }
    if(Fvp->getData_format() != sizeof(T))
    {
        rs_error("ModelViscoelastic3D::Numerical precision in Vp, Vs and Density model files mismatch with constructor.");
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
    Fqp->close();
    Fqs->close();

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
    M = (T *) calloc(1,1);
    L2M = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    M_yz = (T *) calloc(1,1);
    M_xy = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

    Qp = (T *) calloc(1,1);
    Qs = (T *) calloc(1,1);
    Tp = (T *) calloc(1,1);
    Ts = (T *) calloc(1,1);
    Ts_xz = (T *) calloc(1,1);
    Ts_yz = (T *) calloc(1,1);
    Ts_xy = (T *) calloc(1,1);
    f0 = _f0;
}

template<typename T>
void ModelViscoelastic3D<T>::readModel() {
    bool status;
    // Get file names
    std::string Vpfile = this->getVpfile();
    std::string Vsfile = this->getVsfile();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::readModel : Error reading from Vp file: ", Vpfile);
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::readModel : Error reading from Vs file: ", Vsfile);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::readModel : Error reading from Density file: ", Rfile);
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from Qp file: ", Qpfile);
        exit(1);
    }
    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic2D::Error reading from Qs file: ", Qsfile);
        exit(1);
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Vp, Vs  and R */
    free(Vp); free(Vs); free(R); free(Qp); free(Qs);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelViscoelastic3D::readModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelViscoelastic3D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelViscoelastic3D::readModel: Failed to allocate memory.");
    Qp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Qp == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    Qs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Qs == NULL) rs_error("ModelViscoelastic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Fvp->read(Vp, nx*ny*nz);
    Fvp->close();

    Fvs->read(Vs, nx*ny*nz);
    Fvs->close();

    Frho->read(R, nx*ny*nz);
    Frho->close();

    Fqp->read(Qp, nx*ny*nz);
    Fqp->close();

    Fqs->read(Qs, nx*ny*nz);
    Fqs->close();
}

template<typename T>
void ModelViscoelastic3D<T>::writeVp() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic3D::writeVp: Model is not allocated.");
    }
    // Get file names
    std::string Vpfile = this->getVpfile();
    // Open files for writting
    std::shared_ptr<File> Fvp (new File());
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
    T *Mod = this->getVp();
    Fvp->write(Mod, nx*ny*nz, 0);
    Fvp->close();
}

template<typename T>
void ModelViscoelastic3D<T>::writeVs() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic3D::writeVs: Model is not allocated.");
    }
    // Get file names
    std::string Vsfile = this->getVsfile();
    // Open files for writting
    std::shared_ptr<File> Fvs (new File());
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
    T *Mod = this->getVs();
    Fvs->write(Mod, nx*ny*nz, 0);
    Fvs->close();
}

template<typename T>
void ModelViscoelastic3D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic3D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*ny*nz, 0);
    Frho->close();
}

template<typename T>
void ModelViscoelastic3D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelViscoelastic3D::staggerModels: Model is not allocated.");
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
    free(M); free(L2M); free(M_xz); free(M_yz); free(M_xy); 
    free(Rx); free(Ry); free(Rz);
    free(Tp); free(Ts); free(Ts_xz); free(Ts_yz); free(Ts_xy);
    M = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    L2M = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_yz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    M_xy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));

    Tp = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ts = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ts_xz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ts_yz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ts_xy = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    
    // Padding
    this->padmodel3d(Rx, R, nx, ny, nz, lpml);
    this->padmodel3d(Ry, R, nx, ny, nz, lpml);
    this->padmodel3d(Rz, R, nx, ny, nz, lpml);
    this->padmodel3d(L2M, Vp, nx, ny, nz, lpml);
    this->padmodel3d(M, Vs, nx, ny, nz, lpml);
    this->padmodel3d(M_xz, Vs, nx, ny, nz, lpml);
    this->padmodel3d(Tp, Qp, nx, ny, nz, lpml);
    this->padmodel3d(Ts, Qs, nx, ny, nz, lpml);
    
    // Computing moduli
    T _rho, _vp, _vs;
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                _vp=L2M[ind_pml(ix,iy,iz)];
                _vs=M[ind_pml(ix,iy,iz)];
                _rho=Rx[ind_pml(ix,iy,iz)];
                L2M[ind_pml(ix,iy,iz)]=_rho*_vp*_vp;
                M[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                M_xz[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                M_yz[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                M_xy[ind_pml(ix,iy,iz)]=_rho*_vs*_vs;
                Tp[ind_pml(ix,iy,iz)] = (Tp[ind_pml(ix,iy,iz)]+1)/(Tp[ind_pml(ix,iy,iz)]-1);
                Ts_xz[ind_pml(ix,iy,iz)] = (Ts[ind_pml(ix,iy,iz)]+1)/(Ts[ind_pml(ix,iy,iz)]-1);
                Ts_yz[ind_pml(ix,iy,iz)] = (Ts[ind_pml(ix,iy,iz)]+1)/(Ts[ind_pml(ix,iy,iz)]-1);
                Ts_xy[ind_pml(ix,iy,iz)] = (Ts[ind_pml(ix,iy,iz)]+1)/(Ts[ind_pml(ix,iy,iz)]-1);
                Ts[ind_pml(ix,iy,iz)] = (Ts[ind_pml(ix,iy,iz)]+1)/(Ts[ind_pml(ix,iy,iz)]-1);
            }
        }
    }

    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy < ny_pml; iy++){
                L2M[ind_pml(ix,iy,iz)] = M_xz[ind_pml(ix,iy,iz)];
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

    this->staggermodel_x(Ts_xz, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(Ts_xz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_y(Ts_yz, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(Ts_yz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_x(Ts_xy, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(Ts_xy, nx_pml, ny_pml, nz_pml);

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                if(Rx[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Ry[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Rz[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                Rx[ind_pml(ix,iy,iz)] = 1.0/Rx[ind_pml(ix,iy,iz)];
                Ry[ind_pml(ix,iy,iz)] = 1.0/Ry[ind_pml(ix,iy,iz)];
                Rz[ind_pml(ix,iy,iz)] = 1.0/Rz[ind_pml(ix,iy,iz)];
            }
        }
    }

    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy<ny_pml; iy++){
                Rx[ind_pml(ix,iy,lpml)] *= 2.0;
                Ry[ind_pml(ix,iy,lpml)] *= 2.0;
            }
        }
    }
}

template<typename T>
void ModelViscoelastic3D<T>::createModel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Vp and R */
    free(Vp); free(Vs); free(R); free(Qp); free(Qs);
    Vp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vp == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Vs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Vs == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Qp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Qp == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    Qs = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Qs == NULL) rs_error("ModelViscoelastic2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
void ModelViscoelastic3D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    // Reallocate necessary variables 
    free(M); free(L2M); free(M_xz); free(M_yz); free(M_xy); 
    free(Rx); free(Ry); free(Rz);
    free(Tp); free(Ts); free(Ts_xz); free(Ts_yz); free(Ts_xy);
    M = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    L2M = (T *) calloc(nx*ny*nz,sizeof(T));
    if(L2M == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    M_xz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_xz == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    M_yz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_yz == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    M_xy = (T *) calloc(nx*ny*nz,sizeof(T));
    if(M_xy == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Ry = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ry == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");

    Tp = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Tp == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Ts = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ts == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Ts_xz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ts_xz == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Ts_yz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ts_yz == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
    Ts_xy = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ts_xy == NULL) rs_error("ModelViscoelastic3D::createModel: Failed to allocate memory.");
}



template<typename T>
std::shared_ptr<ModelViscoelastic3D<T>> ModelViscoelastic3D<T>::getLocal(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

    std::shared_ptr<ModelViscoelastic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<ModelViscoelastic3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz(), this->getF0(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();
    T *Qp = local->getQp();
    T *Qs = local->getQs();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelViscoelastic3D::getLocal: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace == NULL) rs_error("ModelViscoelastic3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelViscoelastic3D::getLocal: Failed to allocate memory.");
    T *qptrace = (T *) calloc(nx*ny, sizeof(T));
    if(qptrace == NULL) rs_error("ModelViscoelastic3D::getLocal: Failed to allocate memory.");
    T *qstrace = (T *) calloc(nx*ny, sizeof(T));
    if(qstrace == NULL) rs_error("ModelViscoelastic3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::getLocal : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::getLocal : Error reading from Vs file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::getLocal : Error reading from Density file.");
    }
    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::getLocal : Error reading from Qp file: ", Qpfile);
    }
    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::getLocal : Error reading from Qs file: ", Qsfile);
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    Index l3d(size_x, size_y, nz);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f3d(0, 0, i1)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelViscoelastic3D::getLocal: Error reading from vp file");
        Fvs->read(vstrace, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic3D::getLocal: Error reading from vs file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic3D::getLocal: Error reading from rho file");
        Fqp->read(qptrace, nx*ny, fpos);
        if(Fqp->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qp file");
        Fqs->read(qstrace, nx*ny, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic2D::getLocal: Error reading from Qs file");
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
                Qp[l3d(i2,i3,i1)] = qptrace[l2d(lpos_x, lpos_y)];
                Qs[l3d(i2,i3,i1)] = qstrace[l2d(lpos_x, lpos_y)];
            }
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(rhotrace);
    free(qptrace);
    free(qstrace);

    return local;
}

template<typename T>
std::shared_ptr<ModelViscoelastic3D<T>> ModelViscoelastic3D<T>::getDomainmodel(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order) {

    std::shared_ptr<ModelViscoelastic3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();
    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;
    int nxd,nyd,nzd;
    int ix0,iy0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);
    (this->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    nxd = (this->getDomain())->getNx_pad();
    nyd = (this->getDomain())->getNy_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iy0 = (this->getDomain())->getIy0();
    iz0 = (this->getDomain())->getIz0();

    T oxl, oyl, ozl; 
    oxl = (ox + (start_x+ix0-lpml)*dx);
    oyl = (oy + (start_y+iy0-lpml)*dy);
    ozl = (oz + (iz0-lpml)*dz);

    /* Create local model */
    local = std::make_shared<ModelViscoelastic3D<T>>(nxd, nyd, nzd, lpml, dx, dy, dz, oxl, oyl, ozl, this->getF0(), this->getFs());
    (local->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

	/* Copying from big model into local model */
    T *Vp = local->getVp();
    T *Vs = local->getVs();
    T *R = local->getR();
    T *Qp = local->getQp();
    T *Qs = local->getQs();
    T *L2M = local->getL2M();
    T *M = local->getM();
    T *M_xz = local->getM_xz();
    T *M_yz = local->getM_yz();
    T *M_xy = local->getM_xy();
    T *Rx = local->getRx();
    T *Ry = local->getRy();
    T *Rz = local->getRz();

    T *Tp = local->getTp();
    T *Ts = local->getTs();

    T *Ts_xz = local->getTs_xz();
    T *Ts_yz = local->getTs_yz();
    T *Ts_xy = local->getTs_xy();

    /* Allocate two traces to read models from file */
    T *vptrace = (T *) calloc(nx*ny, sizeof(T));
    if(vptrace == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *vstrace = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *vstrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(vstrace_adv == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");

    T *qptrace = (T *) calloc(nx*ny, sizeof(T));
    if(qptrace == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *qstrace = (T *) calloc(nx*ny, sizeof(T));
    if(qstrace == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");
    T *qstrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(qstrace_adv == NULL) rs_error("ModelViscoelastic3D::getDomainmodel: Failed to allocate memory.");


    // Open files for reading
    bool status;
    std::shared_ptr<File> Fvp (new File());
    status = Fvp->input(Vpfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::getDomainmodel : Error reading from Vp file.");
    }
    std::shared_ptr<File> Fvs (new File());
    status = Fvs->input(Vsfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::getDomainmodel : Error reading from Density file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelViscoelastic3D::getDomainmodel : Error reading from Density file.");
    }

    std::shared_ptr<File> Fqp (new File());
    status = Fqp->input(Qpfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::getDomainmodel : Error reading from Qp file: ", Qpfile);
    }
    std::shared_ptr<File> Fqs (new File());
    status = Fqs->input(Qsfile);
    if(status == FILE_ERR){
        rs_error("ModelViscoelastic3D::getDomainmodel : Error reading from Qs file: ", Qsfile);
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, lpos_z, fpos;
    Index l3d(nxd, nyd, nzd);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nzd; i1++) {
        lpos_z = iz0 + i1 - lpml;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);
        Fvp->read(vptrace, nx*ny, fpos);
        if(Fvp->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from vp file");
        Fvs->read(vstrace, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from rho file");

        Fqp->read(qptrace, nx*ny, fpos);
        if(Fqp->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from Qp file");
        Fqs->read(qstrace, nx*ny, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic2D::getDomainmodel: Error reading from Qs file");

        lpos_z = iz0 + i1 - lpml + 1;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);

        Fvs->read(vstrace_adv, nx*ny, fpos);
        if(Fvs->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from vs file");
        Frho->read(rhotrace_adv, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from rho file");
        Fqs->read(qstrace_adv, nx*ny, fpos);
        if(Fqs->getFail()) rs_error("ModelViscoelastic3D::getDomainmodel: Error reading from qs file");

        T M1, M2, M3, M4, M5, M6, M7;
        for(size_t i3=0; i3<nyd; i3++) {
           lpos_y = j + i3 + iy0 - lpml;
           if(lpos_y < 0) lpos_y = 0;
           if(lpos_y > (ny-1)) lpos_y = ny - 1;
           for(size_t i2=0; i2<nxd; i2++) {
              lpos_x = i + i2 + ix0 - lpml;
              if(lpos_x < 0) lpos_x = 0;
              if(lpos_x > (nx-1)) lpos_x = nx - 1;
              Vp[l3d(i2,i3,i1)] = vptrace[l2d(lpos_x, lpos_y)];
              Vs[l3d(i2,i3,i1)] = vstrace[l2d(lpos_x, lpos_y)];
              R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
              Qp[l3d(i2,i3,i1)] = qptrace[l2d(lpos_x, lpos_y)];
              Qs[l3d(i2,i3,i1)] = qstrace[l2d(lpos_x, lpos_y)];
              L2M[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)]*vptrace[l2d(lpos_x, lpos_y)];
              M[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
              Tp[l3d(i2,i3,i1)] = (qptrace[l2d(lpos_x, lpos_y)] + 1)/(qptrace[l2d(lpos_x, lpos_y)] - 1);
              Ts[l3d(i2,i3,i1)] = (qstrace[l2d(lpos_x, lpos_y)] + 1)/(qstrace[l2d(lpos_x, lpos_y)] - 1);

              if(rhotrace[l2d(lpos_x, lpos_y)] <= 0.0) rs_error("ModelViscoelastic3D::getDomainmodel: Zero density found.");

              if(lpos_x < nx-1){
                 Rx[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x+1, lpos_y)]);
              }else{
                 Rx[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
              }
              if(lpos_y < ny-1){
                 Ry[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x, lpos_y+1)]);
              }else{
                 Ry[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
              }
              Rz[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace_adv[l2d(lpos_x, lpos_y)]);
              M1 = rhotrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)]*vstrace[l2d(lpos_x, lpos_y)];
              M5 = rhotrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)]*vstrace_adv[l2d(lpos_x, lpos_y)];

              if(lpos_x < nx-1){
                 M2 = rhotrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)]*vstrace[l2d(lpos_x+1, lpos_y)];
                 M6 = rhotrace_adv[l2d(lpos_x+1, lpos_y)]*vstrace_adv[l2d(lpos_x+1, lpos_y)]*vstrace_adv[l2d(lpos_x+1, lpos_y)];
                 if(lpos_y < ny-1){
                    M4 = rhotrace[l2d(lpos_x+1, lpos_y+1)]*vstrace[l2d(lpos_x+1, lpos_y+1)]*vstrace[l2d(lpos_x+1, lpos_y+1)];
                 }else{
                    M4 = M2;
                 }
              }else{
                 M2 = M1;
                 M6 = M5;
                 if(lpos_y < ny-1){
                    M4 = rhotrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)];
                 }else{
                    M4 = M1;
                 }
              }

              if(lpos_y < ny-1){
                 M3 = rhotrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)]*vstrace[l2d(lpos_x, lpos_y+1)];
                 M7 = rhotrace_adv[l2d(lpos_x, lpos_y+1)]*vstrace_adv[l2d(lpos_x, lpos_y+1)]*vstrace_adv[l2d(lpos_x, lpos_y+1)];
              }else{
                 M3 = M1;
                 M7 = M5;
              }
              M_xz[l3d(i2,i3,i1)] = 0.25*(M1+M2+M5+M6);
              M_yz[l3d(i2,i3,i1)] = 0.25*(M1+M3+M5+M7);
              M_xy[l3d(i2,i3,i1)] = 0.25*(M1+M2+M3+M4);


              M1 = (qstrace[l2d(lpos_x, lpos_y)] + 1)/(qstrace[l2d(lpos_x, lpos_y)] - 1);
              M5 = (qstrace_adv[l2d(lpos_x, lpos_y)] + 1)/(qstrace_adv[l2d(lpos_x, lpos_y)] - 1);

              if(lpos_x < nx-1){
                 M2 = (qstrace[l2d(lpos_x+1, lpos_y)] + 1)/(qstrace[l2d(lpos_x+1, lpos_y)] - 1);
                 M6 = (qstrace_adv[l2d(lpos_x+1, lpos_y)] + 1)/(qstrace_adv[l2d(lpos_x+1, lpos_y)] - 1);
                 if(lpos_y < ny-1){
                    M4 = (qstrace[l2d(lpos_x+1, lpos_y+1)] + 1)/(qstrace[l2d(lpos_x+1, lpos_y+1)] - 1);
                 }else{
                    M4 = M2;
                 }
              }else{
                 M2 = M1;
                 M6 = M5;
                 if(lpos_y < ny-1){
                    M4 = (qstrace[l2d(lpos_x, lpos_y+1)] + 1)/(qstrace[l2d(lpos_x, lpos_y+1)] - 1);
                 }else{
                    M4 = M1;
                 }
              }

              if(lpos_y < ny-1){
                 M3 = (qstrace[l2d(lpos_x, lpos_y+1)] + 1)/(qstrace[l2d(lpos_x, lpos_y+1)] - 1);
                 M7 = (qstrace_adv[l2d(lpos_x, lpos_y+1)] + 1)/(qstrace_adv[l2d(lpos_x, lpos_y+1)] - 1);
              }else{
                 M3 = M1;
                 M7 = M5;
              }
              Ts_xz[l3d(i2,i3,i1)] = 0.25*(M1+M2+M5+M6);
              Ts_yz[l3d(i2,i3,i1)] = 0.25*(M1+M3+M5+M7);
              Ts_xy[l3d(i2,i3,i1)] = 0.25*(M1+M2+M3+M4);
           }
        }
    }

    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            for(size_t iy=0; iy<nyd; iy++){
                Rx[l3d(ix,iy,lpml-iz0)] *= 2.0;
                Ry[l3d(ix,iy,lpml-iz0)] *= 2.0;
                L2M[l3d(ix,iy,lpml-iz0)] = M[l3d(ix,iy,lpml-iz0)];
            }
        }
    }

    /* Free traces */
    free(vptrace);
    free(vstrace);
    free(vstrace_adv);
    free(rhotrace);
    free(rhotrace_adv);
    free(qptrace);
    free(qstrace);
    free(qstrace_adv);

    return local;
}





template<typename T>
ModelViscoelastic3D<T>::~ModelViscoelastic3D() {
    free(Vp);
    free(Vs);
    free(R);
    free(Qp);
    free(Qs);
    free(M);
    free(L2M);
    free(Tp);
    free(Ts);
    free(Ts_xz);
    free(Ts_yz);
    free(Ts_xy);
    free(M_xz);
    free(M_yz);
    free(M_xy);
    free(Rx);
    free(Ry);
    free(Rz);
}

// =============== 2D VTI MODEL CLASS =============== //
template<typename T>
ModelVti2D<T>::ModelVti2D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelVti2D<T>::ModelVti2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    
    /* Allocate variables */
    C11 = (T *) calloc(1,1);
    C13 = (T *) calloc(1,1);
    C33 = (T *) calloc(1,1);
    C55 = (T *) calloc(1,1);
    C11p = (T *) calloc(1,1);
    C13p = (T *) calloc(1,1);
    C33p = (T *) calloc(1,1);
    C55p = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

}

template<typename T>
ModelVti2D<T>::ModelVti2D(std::string _C11file, std::string _C13file, std::string _C33file, std::string _C55file, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    C11file = _C11file;
    C13file = _C13file;
    C33file = _C33file;
    C55file = _C55file;
    Rfile = _Rfile;

    std::shared_ptr<File> FC11 (new File());
    status = FC11->input(C11file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::Error reading from C11 file: ", C11file);
	    exit(1);
    }
    std::shared_ptr<File> FC13 (new File());
    status = FC13->input(C13file.c_str());
    if(status == FILE_ERR){
        rs_error("ModelVti2D::Error reading from C13 file: ", C13file);
    }
    std::shared_ptr<File> FC33 (new File());
    status = FC33->input(C33file.c_str());
    if(status == FILE_ERR){
        rs_error("ModelVti2D::Error reading from C33 file: ", C33file);
    }
    std::shared_ptr<File> FC55 (new File());
    status = FC55->input(C55file.c_str());
    if(status == FILE_ERR){
        rs_error("ModelVti2D::Error reading from C55 file: ", C55file);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelVti2D::Error reading from density file: ", Rfile);
        exit(1);
    }

    // Compare geometry in the two files
    if(FC11->compareGeometry(FC55) != 0)
    {
        rs_error("ModelVti2D::Geometries in C11 and C55 model files do not match.");
    }

    if(FC11->compareGeometry(FC13) != 0)
    {
        rs_error("ModelVti2D::Geometries in C11 and C13 model files do not match.");
    }

    if(FC11->compareGeometry(FC33) != 0)
    {
        rs_error("ModelVti2D::Geometries in C11 and C33 model files do not match.");
    }

    if(FC11->compareGeometry(Frho) != 0)
    {
        rs_error("ModelVti2D::Geometries in C11 and Density model files do not match.");
    }

    if(FC11->getData_format() != Frho->getData_format())
    {
        rs_error("ModelVti2D::Numerical precision mismatch in C11 and Density model files.");
    }
    if(FC11->getData_format() != FC13->getData_format())
    {
        rs_error("ModelVti2D::Numerical precision mismatch in C11 and C13 model files.");
    }
    if(FC11->getData_format() != FC33->getData_format())
    {
        rs_error("ModelVti2D::Numerical precision mismatch in C11 and C33 model files.");
    }
    if(FC11->getData_format() != FC55->getData_format())
    {
        rs_error("ModelVti2D::Numerical precision mismatch in C11 and C55 model files.");
    }
    if(FC11->getData_format() != sizeof(T))
    {
        rs_error("ModelVti2D::Numerical precision in C11, C55 and Density model files mismatch with constructor.");
    }
 

    
    // Read geometry from file
    nx = FC11->getN(1);
    dx = (T) FC11->getD(1);
    ox = (T) FC11->getO(1);
    nz = FC11->getN(3);
    dz = (T) FC11->getD(3);
    oz = (T) FC11->getO(3);
    
    // Close files
    FC11->close();
    FC13->close();
    FC33->close();
    FC55->close();
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
    C11 = (T *) calloc(1,1);
    C13 = (T *) calloc(1,1);
    C33 = (T *) calloc(1,1);
    C55 = (T *) calloc(1,1);
    C11p = (T *) calloc(1,1);
    C13p = (T *) calloc(1,1);
    C33p = (T *) calloc(1,1);
    C55p = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);

}

template<typename T>
void ModelVti2D<T>::readModel() {
    bool status;
    // Get file names
    std::string C11file = this->getC11file();
    std::string C13file = this->getC13file();
    std::string C33file = this->getC33file();
    std::string C55file = this->getC55file();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<File> FC11 (new File());
    status = FC11->input(C11file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::readModel : Error reading from C11 file: ", C11file);
    }
    std::shared_ptr<File> FC13 (new File());
    status = FC13->input(C13file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::readModel : Error reading from C13 file: ", C13file);
    }
    std::shared_ptr<File> FC33 (new File());
    status = FC33->input(C33file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::readModel : Error reading from C33 file: ", C33file);
    }
    std::shared_ptr<File> FC55 (new File());
    status = FC55->input(C55file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::readModel : Error reading from C55 file: ", C55file);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int nz = this->getNz();
    
    /* Reallocate C11, C13, C33, C55  and R */
    free(C11); free(C13); free(C33); free(C55); free(R);
    C11 = (T *) calloc(nx*nz,sizeof(T));
    if(C11 == NULL) rs_error("ModelVti2D::readModel: Failed to allocate memory.");
    C13 = (T *) calloc(nx*nz,sizeof(T));
    if(C13 == NULL) rs_error("ModelVti2D::readModel: Failed to allocate memory.");
    C33 = (T *) calloc(nx*nz,sizeof(T));
    if(C33 == NULL) rs_error("ModelVti2D::readModel: Failed to allocate memory.");
    C55 = (T *) calloc(nx*nz,sizeof(T));
    if(C55 == NULL) rs_error("ModelVti2D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelVti2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    FC11->read(C11, nx*nz);
    FC11->close();

    FC13->read(C13, nx*nz);
    FC13->close();

    FC33->read(C33, nx*nz);
    FC33->close();

    FC55->read(C55, nx*nz);
    FC55->close();

    Frho->read(R, nx*nz);
    Frho->close();
}

template<typename T>
void ModelVti2D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelVti2D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*nz, 0);
    Frho->close();
}


template<typename T>
void ModelVti2D<T>::writeC11() {
    if(!this->getRealized()) {
        rs_error("ModelVti2D::writeC11: Model is not allocated.");
    }
    // Get file names
    std::string C11file = this->getC11file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C11file);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC11();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelVti2D<T>::writeC13() {
    if(!this->getRealized()) {
        rs_error("ModelVti2D::writeC13: Model is not allocated.");
    }
    // Get file names
    std::string C13file = this->getC13file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C13file);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC13();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelVti2D<T>::writeC33() {
    if(!this->getRealized()) {
        rs_error("ModelVti2D::writeC33: Model is not allocated.");
    }
    // Get file names
    std::string C33file = this->getC33file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C33file);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC33();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelVti2D<T>::writeC55() {
    if(!this->getRealized()) {
        rs_error("ModelVti2D::writeC55: Model is not allocated.");
    }
    // Get file names
    std::string C55file = this->getC55file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C55file);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC55();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelVti2D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelVti2D::staggerModels: Model is not allocated.");
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
    free(C11p); free(C13p); free(C33p); free(C55p); free(Rx); free(Rz);
    C11p = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(C11p == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");
    C13p = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(C13p == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");
    C33p = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(C33p == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");
    C55p = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(C55p == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");
    Rx = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rx == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");
    Rz = (T *) calloc(nx_pml*nz_pml,sizeof(T));
    if(Rz == NULL) rs_error("ModelVti2D::staggerModels: Failed to allocate memory.");

    // Padding
    this->padmodel2d(Rx, R, nx, nz, lpml);
    this->padmodel2d(Rz, R, nx, nz, lpml);
    this->padmodel2d(C11p, C11, nx, nz, lpml);
    this->padmodel2d(C13p, C13, nx, nz, lpml);
    this->padmodel2d(C33p, C33, nx, nz, lpml);
    this->padmodel2d(C55p, C55, nx, nz, lpml);
    
    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            C13p[ind_pml(ix,iz)] = 0.0;
            C11p[ind_pml(ix,iz)] = C55p[ind_pml(ix,iz)];
        }
    }

    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, 1, nz_pml); 
    this->staggermodel_z(Rz, nx_pml, 1, nz_pml); 
    
    this->staggermodel_z(C55p, nx_pml, 1, nz_pml); 
    this->staggermodel_x(C55p, nx_pml, 1, nz_pml); 

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iz=0; iz < nz_pml; iz++){
            if(Rx[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            if(Rz[ind_pml(ix,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
            Rx[ind_pml(ix,iz)] = 1.0/Rx[ind_pml(ix,iz)];
            Rz[ind_pml(ix,iz)] = 1.0/Rz[ind_pml(ix,iz)];
        }
    }
    
    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            Rx[ind_pml(ix,iz)] *= 2.0;
        }
    }
}

template<typename T>
void ModelVti2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate C11, C13, C33, C55 and R */
    free(C11); free(C13); free(C33); free(C55); free(R);
    C11 = (T *) calloc(nx*nz,sizeof(T));
    if(C11 == NULL) rs_error("ModelVti2D::createModel: Failed to allocate memory.");
    C13 = (T *) calloc(nx*nz,sizeof(T));
    if(C13 == NULL) rs_error("ModelVti2D::createModel: Failed to allocate memory.");
    C33 = (T *) calloc(nx*nz,sizeof(T));
    if(C33 == NULL) rs_error("ModelVti2D::createModel: Failed to allocate memory.");
    C55 = (T *) calloc(nx*nz,sizeof(T));
    if(C55 == NULL) rs_error("ModelVti2D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*nz,sizeof(T));
    if(R == NULL) rs_error("ModelVti2D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
void ModelVti2D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate C11p, C13p ,C33p, C55p, Rx, and Rz */
    free(C11p); free(C13p); free(C33p); free(C55p); free(Rx); free(Rz);

    C11p = (T *) calloc(nx*nz,sizeof(T));
    if(C11p == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
    C13p = (T *) calloc(nx*nz,sizeof(T));
    if(C13p == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
    C33p = (T *) calloc(nx*nz,sizeof(T));
    if(C33p == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
    C55p = (T *) calloc(nx*nz,sizeof(T));
    if(C55p == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
    Rx = (T *) calloc(nx*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelVti2D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelVti2D<T>> ModelVti2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {

    std::shared_ptr<ModelVti2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelVti2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *C11 = local->getC11();
    T *C13 = local->getC13();
    T *C33 = local->getC33();
    T *C55 = local->getC55();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *C11trace = (T *) calloc(nx, sizeof(T));
    if(C11trace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *C13trace = (T *) calloc(nx, sizeof(T));
    if(C13trace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *C33trace = (T *) calloc(nx, sizeof(T));
    if(C33trace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *C55trace = (T *) calloc(nx, sizeof(T));
    if(C55trace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> FC11 (new File());
    status = FC11->input(C11file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from C11 file.");
    }
    std::shared_ptr<File> FC13 (new File());
    status = FC13->input(C13file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from C13 file.");
    }
    std::shared_ptr<File> FC33 (new File());
    status = FC33->input(C33file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from C33 file.");
    }
    std::shared_ptr<File> FC55 (new File());
    status = FC55->input(C55file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from C55 file.");
    }

    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        FC11->read(C11trace, nx, fpos);
        if(FC11->getFail()) rs_error("ModelVti2D::getLocal: Error reading from C11 file");
        FC13->read(C13trace, nx, fpos);
        if(FC13->getFail()) rs_error("ModelVti2D::getLocal: Error reading from C13 file");
        FC33->read(C33trace, nx, fpos);
        if(FC33->getFail()) rs_error("ModelVti2D::getLocal: Error reading from C33 file");
        FC55->read(C55trace, nx, fpos);
        if(FC55->getFail()) rs_error("ModelVti2D::getLocal: Error reading from C55 file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelVti2D::getLocal: Error reading from rho file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            C11[l2d(i2,i1)] = C11trace[lpos];
            C13[l2d(i2,i1)] = C13trace[lpos];
            C33[l2d(i2,i1)] = C33trace[lpos];
            C55[l2d(i2,i1)] = C55trace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
        }
    }

    /* Free traces */
    free(C11trace);
    free(C13trace);
    free(C33trace);
    free(C55trace);
    free(rhotrace);

    return local;
}

template<typename T>
std::shared_ptr<ModelVti2D<T>> ModelVti2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd0, const int nd1, const int order) {
    std::shared_ptr<ModelVti2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*this->getLpml(),d,nd0,1,nd1,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelVti2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getFs());
    (local->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*lpml,d,nd0,1,nd1,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *C11 = local->getC11();
    T *C13 = local->getC13();
    T *C33 = local->getC33();
    T *C55 = local->getC55();
    T *C11p = local->getC11p();
    T *C13p = local->getC13p();
    T *C33p = local->getC33p();
    T *C55p = local->getC55p();
    T *R = local->getR();
    T *Rx = local->getRx();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *C11trace = (T *) calloc(nx, sizeof(T));
    if(C11trace == NULL) rs_error("ModelVti2d::getDomainmodel: Failed to allocate memory.");
    T *C13trace = (T *) calloc(nx, sizeof(T));
    if(C13trace == NULL) rs_error("ModelVti2d::getDomainmodel: Failed to allocate memory.");
    T *C33trace = (T *) calloc(nx, sizeof(T));
    if(C33trace == NULL) rs_error("ModelVti2d::getDomainmodel: Failed to allocate memory.");
    T *C55trace = (T *) calloc(nx, sizeof(T));
    if(C55trace == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *C55trace_adv = (T *) calloc(nx, sizeof(T));
    if(C55trace_adv == NULL) rs_error("ModelVti2d::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelVti2d::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelVti2d::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> FC11 (new File());
    status = FC11->input(C11file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getDomainmodel : Error reading from C11 file.");
    }
    std::shared_ptr<File> FC13 (new File());
    status = FC13->input(C13file);
    if(status == FILE_ERR){
       rs_error("ModelVti2D::getDomainmodel : Error reading from C13 file.");
    }
    std::shared_ptr<File> FC33 (new File());
    status = FC33->input(C33file);
    if(status == FILE_ERR){
       rs_error("ModelVti2D::getDomainmodel : Error reading from C33 file.");
    }
    std::shared_ptr<File> FC55 (new File());
    status = FC55->input(C55file);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getLocal : Error reading from C55 file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelVti2D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);
    T M1, M2, M3, M4;

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        FC11->read(C11trace, nx, fpos);
        if(FC11->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from C11 file");
        FC13->read(C13trace, nx, fpos);
        if(FC13->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from C13 file");
        FC33->read(C33trace, nx, fpos);
        if(FC33->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from C33 file");
        FC55->read(C55trace, nx, fpos);
        if(FC55->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from C55 file");
        Frho->read(rhotrace, nx, fpos);
        if(Frho->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from rho file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        FC55->read(C55trace_adv, nx, fpos);
        if(FC55->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from C55 file");
        Frho->read(rhotrace_adv, nx, fpos);
        if(Frho->getFail()) rs_error("ModelVti2D::getDomainmodel: Error reading from rho file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            C11[l2d(i2,i1)] = C11trace[lpos];
            C13[l2d(i2,i1)] = C13trace[lpos];
            C33[l2d(i2,i1)] = C33trace[lpos];
            C55[l2d(i2,i1)] = C55trace[lpos];
            R[l2d(i2,i1)] = rhotrace[lpos];
            C11p[l2d(i2,i1)] = C11trace[lpos];
            C13p[l2d(i2,i1)] = C13trace[lpos];
            C33p[l2d(i2,i1)] = C33trace[lpos];
            if(rhotrace[lpos] <= 0.0) rs_error("ModelVti2D::getDomainmodel: Zero density found.");
            if(lpos < nx-1){
               Rx[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace[lpos+1]);
            }else{
               Rx[l2d(i2,i1)] = 1.0/(rhotrace[lpos]);
            }
            Rz[l2d(i2,i1)] = 2.0/(rhotrace[lpos]+rhotrace_adv[lpos]);

            M1 = C55trace[lpos];
            M3 = C55trace_adv[lpos];
            if(lpos < nx-1){
               M2 = C55trace[lpos+1];
               M4 = C55trace_adv[lpos+1];
            }else{
               M2 = C55trace[lpos];
               M4 = C55trace_adv[lpos];
            }
            C55p[l2d(i2,i1)] = 0.25*(M1+M2+M3+M4);
        }
    }

       
    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rx[l2d(ix,lpml-iz0)] *= 2.0;
            C11p[l2d(ix,lpml-iz0)] = (C11p[l2d(ix,lpml-iz0)] - C13p[l2d(ix,lpml-iz0)])/2.0;
            C13p[l2d(ix,lpml-iz0)] = 0.0;
        }
    }


    /* Free traces */
    free(C11trace);
    free(C13trace);
    free(C33trace);
    free(C55trace);
    free(C55trace_adv);
    free(rhotrace);
    free(rhotrace_adv);
    
    return local;
}


template<typename T>
ModelVti2D<T>::~ModelVti2D() {
    free(C11);
    free(C13);
    free(C33);
    free(C55);
    free(C11p);
    free(C13p);
    free(C33p);
    free(C55p);
    free(R);
    free(Rx);
    free(Rz);
}

// =============== 3D ORTHOROMBIC MODEL CLASS =============== //
template<typename T>
ModelOrtho3D<T>::ModelOrtho3D(): Model<T>(3) {
    // Nothing here
}

template<typename T>
ModelOrtho3D<T>::ModelOrtho3D(const int _nx, const int _ny, const int _nz, const int _lpml, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz, const bool _fs): Model<T>(3, _nx, _ny, _nz, _lpml, _dx, _dy, _dz, _ox, _oy, _oz, _fs) {
    
    /* Allocate variables */
    C11 = (T *) calloc(1,1);
    C12 = (T *) calloc(1,1);
    C13 = (T *) calloc(1,1);
    C22 = (T *) calloc(1,1);
    C23 = (T *) calloc(1,1);
    C33 = (T *) calloc(1,1);
    C44 = (T *) calloc(1,1);
    C55 = (T *) calloc(1,1);
    C66 = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    C11p = (T *) calloc(1,1);
    C12p = (T *) calloc(1,1);
    C13p = (T *) calloc(1,1);
    C22p = (T *) calloc(1,1);
    C23p = (T *) calloc(1,1);
    C33p = (T *) calloc(1,1);
    C44p = (T *) calloc(1,1);
    C55p = (T *) calloc(1,1);
    C66p = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
ModelOrtho3D<T>::ModelOrtho3D(std::string _C11file, std::string _C12file, std::string _C13file, std::string _C22file, std::string _C23file, std::string _C33file, std::string _C44file, std::string _C55file, std::string _C66file, std::string _Rfile, const int _lpml, const bool _fs): Model<T>(3) {
    bool status;
    int nx, ny, nz;
    T dx, dy, dz;
    T ox, oy, oz;
    C11file = _C11file;
    C12file = _C12file;
    C13file = _C13file;
    C22file = _C22file;
    C23file = _C23file;
    C33file = _C33file;
    C44file = _C44file;
    C55file = _C55file;
    C66file = _C66file;
    Rfile = _Rfile;

    std::shared_ptr<File> Fc11 (new File());
    status = Fc11->input(C11file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C11 file:  ", _C11file);
    }
    std::shared_ptr<File> Fc12 (new File());
    status = Fc12->input(C12file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C12 file:  ", _C12file);
    }
    std::shared_ptr<File> Fc13 (new File());
    status = Fc13->input(C13file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C13 file:  ", _C13file);
    }
    std::shared_ptr<File> Fc22 (new File());
    status = Fc22->input(C22file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C22 file:  ", _C22file);
    }
    std::shared_ptr<File> Fc23 (new File());
    status = Fc23->input(C23file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C23 file:  ", _C23file);
    }
    std::shared_ptr<File> Fc33 (new File());
    status = Fc33->input(C33file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C33 file:  ", _C33file);
    }
    std::shared_ptr<File> Fc44 (new File());
    status = Fc44->input(C44file.c_str());
    if(status == FILE_ERR){
        rs_error("ModelOrtho3D::Error reading from C44 file: ", _C44file);
    }
    std::shared_ptr<File> Fc55 (new File());
    status = Fc55->input(C55file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C55 file:  ", _C55file);
    }
    std::shared_ptr<File> Fc66 (new File());
    status = Fc66->input(C66file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::Error reading from C66 file:  ", _C66file);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
        rs_error("ModelOrtho3D::Error reading from density file: ", _Rfile);
    }

    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc12) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C12 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc13) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C13 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc22) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C22 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc23) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C23 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc33) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C33 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc44) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C44 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc55) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C55 model files do not match.");
    }
    // Compare geometry in the two files
    if(Fc11->compareGeometry(Fc66) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and C66 model files do not match.");
    }
    if(Fc11->compareGeometry(Frho) != 0)
    {
        rs_error("ModelOrtho3D::Geometries in C11 and Density model files do not match.");
    }

    if(Fc11->getData_format() != Frho->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and Density model files.");
    }
    if(Fc11->getData_format() != Fc12->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C12 model files.");
    }
    if(Fc11->getData_format() != Fc13->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C13 model files.");
    }
    if(Fc11->getData_format() != Fc23->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C23 model files.");
    }
    if(Fc11->getData_format() != Fc33->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C33 model files.");
    }
    if(Fc11->getData_format() != Fc44->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C44 model files.");
    }
    if(Fc11->getData_format() != Fc55->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C55 model files.");
    }
    if(Fc11->getData_format() != Fc66->getData_format())
    {
        rs_error("ModelOrtho3D::Numerical precision mismatch in C11 and C66 model files.");
    }
    if(Fc11->getData_format() != sizeof(T))
    {
        rs_error("ModelOrtho3D::Numerical precision in C11, C44 and Density model files mismatch with constructor.");
    }

    // Read geometry from file
    nx = Fc11->getN(1);
    dx = (T) Fc11->getD(1);
    ox = (T) Fc11->getO(1);
    ny = Fc11->getN(2);
    dy = (T) Fc11->getD(2);
    oy = (T) Fc11->getO(2);
    nz = Fc11->getN(3);
    dz = (T) Fc11->getD(3);
    oz = (T) Fc11->getO(3);
    
    // Close files
    Fc11->close();
    Fc12->close();
    Fc13->close();
    Fc22->close();
    Fc23->close();
    Fc33->close();
    Fc44->close();
    Fc55->close();
    Fc66->close();
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
    C11 = (T *) calloc(1,1);
    C12 = (T *) calloc(1,1);
    C13 = (T *) calloc(1,1);
    C22 = (T *) calloc(1,1);
    C23 = (T *) calloc(1,1);
    C33 = (T *) calloc(1,1);
    C44 = (T *) calloc(1,1);
    C55 = (T *) calloc(1,1);
    C66 = (T *) calloc(1,1);
    R = (T *) calloc(1,1);
    C11p = (T *) calloc(1,1);
    C12p = (T *) calloc(1,1);
    C13p = (T *) calloc(1,1);
    C22p = (T *) calloc(1,1);
    C23p = (T *) calloc(1,1);
    C33p = (T *) calloc(1,1);
    C44p = (T *) calloc(1,1);
    C55p = (T *) calloc(1,1);
    C66p = (T *) calloc(1,1);
    Rx = (T *) calloc(1,1);
    Ry = (T *) calloc(1,1);
    Rz = (T *) calloc(1,1);
}

template<typename T>
void ModelOrtho3D<T>::readModel() {
    bool status;
    // Get file names
    std::string C11file = this->getC11file();
    std::string C12file = this->getC12file();
    std::string C13file = this->getC13file();
    std::string C22file = this->getC22file();
    std::string C23file = this->getC23file();
    std::string C33file = this->getC33file();
    std::string C44file = this->getC44file();
    std::string C55file = this->getC55file();
    std::string C66file = this->getC66file();
    std::string Rfile = this->getRfile();
    // Open files for reading
    std::shared_ptr<File> Fc11 (new File());
    status = Fc11->input(C11file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C11 file: ", C11file);
    }
    std::shared_ptr<File> Fc12 (new File());
    status = Fc12->input(C12file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C12 file: ", C12file);
    }
    std::shared_ptr<File> Fc13 (new File());
    status = Fc13->input(C13file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C13 file: ", C13file);
    }
    std::shared_ptr<File> Fc22 (new File());
    status = Fc22->input(C22file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C22 file: ", C22file);
    }
    std::shared_ptr<File> Fc23 (new File());
    status = Fc23->input(C23file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C23 file: ", C23file);
    }
    std::shared_ptr<File> Fc33 (new File());
    status = Fc33->input(C33file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C33 file: ", C33file);
    }
    std::shared_ptr<File> Fc44 (new File());
    status = Fc44->input(C44file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C44 file: ", C44file);
    }
    std::shared_ptr<File> Fc55 (new File());
    status = Fc55->input(C55file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C55 file: ", C55file);
    }
    std::shared_ptr<File> Fc66 (new File());
    status = Fc66->input(C66file.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from C66 file: ", C66file);
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::readModel : Error reading from Density file: ", Rfile);
    }

    // Read models
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate C11, C12, C13, C22, C23, C44, C55, C66 and R */
    free(C11); free(C12); free(C13); free(C22); free(C23); free(C33); free(C44); free(C55); free(C66); free(R);
    C11 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C11 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C12 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C12 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C13 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C13 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C22 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C22 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C23 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C23 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C33 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C33 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C44 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C44 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C55 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C55 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    C66 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C66 == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelOrtho3D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    Fc11->read(C11, nx*ny*nz);
    Fc11->close();

    Fc12->read(C12, nx*ny*nz);
    Fc12->close();

    Fc13->read(C13, nx*ny*nz);
    Fc13->close();

    Fc22->read(C22, nx*ny*nz);
    Fc22->close();

    Fc23->read(C23, nx*ny*nz);
    Fc23->close();

    Fc33->read(C33, nx*ny*nz);
    Fc33->close();

    Fc44->read(C44, nx*ny*nz);
    Fc44->close();

    Fc55->read(C55, nx*ny*nz);
    Fc55->close();

    Fc66->read(C66, nx*ny*nz);
    Fc66->close();

    Frho->read(R, nx*ny*nz);
    Frho->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC11() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC11: Model is not allocated.");
    }
    // Get file names
    std::string C11file = this->getC11file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C11file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC11();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC12() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC12: Model is not allocated.");
    }
    // Get file names
    std::string C12file = this->getC12file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C12file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC12();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC13() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC13: Model is not allocated.");
    }
    // Get file names
    std::string C13file = this->getC13file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C13file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC13();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC22() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC22: Model is not allocated.");
    }
    // Get file names
    std::string C22file = this->getC22file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C22file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC22();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC23() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC23: Model is not allocated.");
    }
    // Get file names
    std::string C23file = this->getC23file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C23file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC23();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC33() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC33: Model is not allocated.");
    }
    // Get file names
    std::string C33file = this->getC33file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C33file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC33();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC44() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC44: Model is not allocated.");
    }
    // Get file names
    std::string C44file = this->getC44file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C44file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC44();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC55() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC55: Model is not allocated.");
    }
    // Get file names
    std::string C55file = this->getC55file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C55file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC55();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeC66() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeC66: Model is not allocated.");
    }
    // Get file names
    std::string C66file = this->getC66file();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(C66file.c_str());

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

    F->setN(1,nx);
    F->setN(2,ny);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(2,dy);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(2,oy);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getC66();
    F->write(Mod, nx*ny*nz, 0);
    F->close();
}

template<typename T>
void ModelOrtho3D<T>::writeR() {
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::writeR: Model is not allocated.");
    }
    // Get file names
    std::string Rfile = this->getRfile();
    // Open files for writting
    std::shared_ptr<File> Frho (new File());
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
    T *Mod = this->getR();
    Frho->write(Mod, nx*ny*nz, 0);
    Frho->close();
}

template<typename T>
void ModelOrtho3D<T>::staggerModels(){
    if(!this->getRealized()) {
        rs_error("ModelOrtho3D::staggerModels: Model is not allocated.");
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
    free(C11p); free(C12p); free(C13p); free(C22p); free(C23p); free(C33p); free(C44p); free(C55p); free(C66p); 
    free(Rx); free(Ry); free(Rz);
    C11 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C12 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C13 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C22 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C23 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C33 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C55 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C44 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    C66 = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rx = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Ry = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    Rz = (T *) calloc(nx_pml*ny_pml*nz_pml,sizeof(T));
    
    // Padding
    this->padmodel3d(Rx, R, nx, ny, nz, lpml);
    this->padmodel3d(Ry, R, nx, ny, nz, lpml);
    this->padmodel3d(Rz, R, nx, ny, nz, lpml);
    this->padmodel3d(C11p, C11, nx, ny, nz, lpml);
    this->padmodel3d(C12p, C12, nx, ny, nz, lpml);
    this->padmodel3d(C13p, C13, nx, ny, nz, lpml);
    this->padmodel3d(C22p, C22, nx, ny, nz, lpml);
    this->padmodel3d(C23p, C23, nx, ny, nz, lpml);
    this->padmodel3d(C33p, C33, nx, ny, nz, lpml);
    this->padmodel3d(C44p, C44, nx, ny, nz, lpml);
    this->padmodel3d(C55p, C55, nx, ny, nz, lpml);
    this->padmodel3d(C66p, C66, nx, ny, nz, lpml);
    
    // In case of free surface
    if(this->getFs()){
        iz=lpml;
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy < ny_pml; iy++){
                C11p[ind_pml(ix,iy,iz)] = C12p[ind_pml(ix,iy,iz)] + C66p[ind_pml(ix,iy,iz)];
                C22p[ind_pml(ix,iy,iz)] = C12p[ind_pml(ix,iy,iz)] + C66p[ind_pml(ix,iy,iz)];
                C66p[ind_pml(ix,iy,iz)] *= 0.5;
            }
        }
    }


    // Staggering using arithmetic average
    this->staggermodel_x(Rx, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(Ry, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(Rz, nx_pml, ny_pml, nz_pml);

    this->staggermodel_x(C55, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(C55, nx_pml, ny_pml, nz_pml);

    this->staggermodel_y(C44, nx_pml, ny_pml, nz_pml);
    this->staggermodel_z(C44, nx_pml, ny_pml, nz_pml);

    this->staggermodel_x(C66, nx_pml, ny_pml, nz_pml);
    this->staggermodel_y(C66, nx_pml, ny_pml, nz_pml);

    // Inverting the density
    for(ix=0; ix < nx_pml; ix++){
        for(iy=0; iy < ny_pml; iy++){
            for(iz=0; iz < nz_pml; iz++){
                if(Rx[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Ry[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                if(Rz[ind_pml(ix,iy,iz)] == 0.0) rs_error("staggerModels: Zero density found.");
                Rx[ind_pml(ix,iy,iz)] = 1.0/Rx[ind_pml(ix,iy,iz)];
                Ry[ind_pml(ix,iy,iz)] = 1.0/Ry[ind_pml(ix,iy,iz)];
                Rz[ind_pml(ix,iy,iz)] = 1.0/Rz[ind_pml(ix,iy,iz)];
            }
        }
    }

    
    // In case of free surface
    if(this->getFs()){
        for(ix=0; ix<nx_pml; ix++){
            for(iy=0; iy<ny_pml; iy++){
                Rx[ind_pml(ix,iy,lpml)] *= 2.0;
                Ry[ind_pml(ix,iy,lpml)] *= 2.0;
            }
        }
    }
}

template<typename T>
void ModelOrtho3D<T>::createModel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate Stiffnesses and R */
    free(C11); free(C12); free(C13); free(C22); free(C23); free(C33); free(C44); free(C55); free(C66); free(R);
    C11 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C11 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C12 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C12 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C13 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C13 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C22 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C22 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C23 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C23 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C33 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C33 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C44 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C44 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C55 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C55 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C66 = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C66 == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    R = (T *) calloc(nx*ny*nz,sizeof(T));
    if(R == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    this->setRealized(true);
}

template<typename T>
void ModelOrtho3D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();

    /* Reallocate L Rx, and Rz */
    free(C11p); free(C12p); free(C13p); free(C22p); free(C23p); free(C33p); free(C44p); free(C55p); free(C66p); 
    free(Rx); free(Ry); free(Rz);
    C11p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C11p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C12p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C12p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C13p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C13p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C22p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C22p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C23p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C23p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C33p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C33p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C44p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C44p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C55p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C55p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");
    C66p = (T *) calloc(nx*ny*nz,sizeof(T));
    if(C66p == NULL) rs_error("ModelOrtho3D::createModel: Failed to allocate memory.");

    Rx = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rx == NULL) rs_error("ModelOrtho3D::createPaddedmodel: Failed to allocate memory.");
    Ry = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Ry == NULL) rs_error("ModelOrtho3D::createPaddedmodel: Failed to allocate memory.");
    Rz = (T *) calloc(nx*ny*nz,sizeof(T));
    if(Rz == NULL) rs_error("ModelOrtho3D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelOrtho3D<T>> ModelOrtho3D<T>::getLocal(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map) {
    std::shared_ptr<ModelOrtho3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T ox = this->getOx();
    T oy = this->getOy();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);

    double oxl, oyl; 
    oxl = (ox + start_x*dx);
    oyl = (oy + start_y*dy);

    /* Create local model */
    local = std::make_shared<ModelOrtho3D<T>>(size_x, size_y, nz, this->getLpml(), dx, dy, this->getDz(), oxl, oyl, this->getOz(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *C11 = local->getC11();
    T *C12 = local->getC12();
    T *C13 = local->getC13();
    T *C22 = local->getC22();
    T *C23 = local->getC23();
    T *C33 = local->getC33();
    T *C44 = local->getC44();
    T *C55 = local->getC55();
    T *C66 = local->getC66();
    T *R = local->getR();

    /* Allocate two traces to read models from file */
    T *c11trace = (T *) calloc(nx*ny, sizeof(T));
    if(c11trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c12trace = (T *) calloc(nx*ny, sizeof(T));
    if(c12trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c13trace = (T *) calloc(nx*ny, sizeof(T));
    if(c13trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c22trace = (T *) calloc(nx*ny, sizeof(T));
    if(c22trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c23trace = (T *) calloc(nx*ny, sizeof(T));
    if(c23trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c33trace = (T *) calloc(nx*ny, sizeof(T));
    if(c33trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c44trace = (T *) calloc(nx*ny, sizeof(T));
    if(c44trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c55trace = (T *) calloc(nx*ny, sizeof(T));
    if(c55trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *c66trace = (T *) calloc(nx*ny, sizeof(T));
    if(c66trace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelOrtho3D::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fc11 (new File());
    status = Fc11->input(C11file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C11 file.");
    }
    std::shared_ptr<File> Fc12 (new File());
    status = Fc12->input(C12file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C12 file.");
    }
    std::shared_ptr<File> Fc13 (new File());
    status = Fc13->input(C13file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C13 file.");
    }
    std::shared_ptr<File> Fc22 (new File());
    status = Fc22->input(C22file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C22 file.");
    }
    std::shared_ptr<File> Fc23 (new File());
    status = Fc23->input(C23file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C23 file.");
    }
    std::shared_ptr<File> Fc33 (new File());
    status = Fc33->input(C33file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C33 file.");
    }
    std::shared_ptr<File> Fc44 (new File());
    status = Fc44->input(C44file);
    if(status == FILE_ERR){
        rs_error("ModelOrtho3D::getLocal : Error reading from C44 file.");
    }
    std::shared_ptr<File> Fc55 (new File());
    status = Fc55->input(C55file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C55 file.");
    }
    std::shared_ptr<File> Fc66 (new File());
    status = Fc66->input(C66file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from C66 file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getLocal : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, fpos;
    Index l3d(size_x, size_y, nz);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f3d(0, 0, i1)*sizeof(T);
        Fc11->read(c11trace, nx*ny, fpos);
        if(Fc11->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c11 file");
        Fc12->read(c12trace, nx*ny, fpos);
        if(Fc12->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c12 file");
        Fc13->read(c13trace, nx*ny, fpos);
        if(Fc13->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c13 file");
        Fc22->read(c22trace, nx*ny, fpos);
        if(Fc22->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c22 file");
        Fc23->read(c23trace, nx*ny, fpos);
        if(Fc23->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c23 file");
        Fc33->read(c33trace, nx*ny, fpos);
        if(Fc33->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c33 file");
        Fc44->read(c44trace, nx*ny, fpos);
        if(Fc44->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c44 file");
        Fc55->read(c55trace, nx*ny, fpos);
        if(Fc55->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c55 file");
        Fc66->read(c66trace, nx*ny, fpos);
        if(Fc66->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from c66 file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelOrtho3D::getLocal: Error reading from rho file");
        for(size_t i3=0; i3<size_y; i3++) {
            lpos_y = j + i3;
            if(lpos_y < 0) lpos_y = 0;
            if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<size_x; i2++) {
                lpos_x = i + i2;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                C11[l3d(i2,i3,i1)] = c11trace[l2d(lpos_x, lpos_y)];
                C12[l3d(i2,i3,i1)] = c12trace[l2d(lpos_x, lpos_y)];
                C13[l3d(i2,i3,i1)] = c13trace[l2d(lpos_x, lpos_y)];
                C22[l3d(i2,i3,i1)] = c22trace[l2d(lpos_x, lpos_y)];
                C23[l3d(i2,i3,i1)] = c23trace[l2d(lpos_x, lpos_y)];
                C33[l3d(i2,i3,i1)] = c33trace[l2d(lpos_x, lpos_y)];
                C44[l3d(i2,i3,i1)] = c44trace[l2d(lpos_x, lpos_y)];
                C55[l3d(i2,i3,i1)] = c55trace[l2d(lpos_x, lpos_y)];
                C66[l3d(i2,i3,i1)] = c66trace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
            }
        }
    }

    /* Free traces */
    free(c11trace);
    free(c12trace);
    free(c13trace);
    free(c22trace);
    free(c23trace);
    free(c33trace);
    free(c44trace);
    free(c55trace);
    free(c66trace);
    free(rhotrace);

    return local;
}

template<typename T>
std::shared_ptr<ModelOrtho3D<T>> ModelOrtho3D<T>::getDomainmodel(std::shared_ptr<Data3D<T>> data, T aperture_x, T aperture_y, bool map, const int d, const int nd0, const int nd1, const int nd2, const int order) {

    std::shared_ptr<ModelOrtho3D<T>> local;
    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();
    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();
    size_t nx = this->getNx();
    size_t ny = this->getNy();
    size_t nz = this->getNz();
    size_t size_x;
    off_t start_x;
    size_t size_y;
    off_t start_y;
    int nxd,nyd,nzd;
    int ix0,iy0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize3d(data, aperture_x, aperture_y, map, &start_x, &size_x, &start_y, &size_y);
    (this->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    nxd = (this->getDomain())->getNx_pad();
    nyd = (this->getDomain())->getNy_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iy0 = (this->getDomain())->getIy0();
    iz0 = (this->getDomain())->getIz0();

    T oxl, oyl, ozl; 
    oxl = (ox + (start_x+ix0-lpml)*dx);
    oyl = (oy + (start_y+iy0-lpml)*dy);
    ozl = (oz + (iz0-lpml)*dz);

    /* Create local model */
    local = std::make_shared<ModelOrtho3D<T>>(nxd, nyd, nzd, lpml, dx, dy, dz, oxl, oyl, ozl, this->getFs());
    (local->getDomain())->setupDomain3D(size_x+2*lpml,size_y+2*lpml,nz+2*lpml,d,nd0,nd1,nd2,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

	/* Copying from big model into local model */
    T *C11 = local->getC11();
    T *C12 = local->getC12();
    T *C13 = local->getC13();
    T *C22 = local->getC22();
    T *C23 = local->getC23();
    T *C33 = local->getC33();
    T *C44 = local->getC44();
    T *C55 = local->getC55();
    T *C66 = local->getC66();

    T *R = local->getR();
    T *C11p = local->getC11p();
    T *C12p = local->getC12p();
    T *C13p = local->getC13p();
    T *C22p = local->getC22p();
    T *C23p = local->getC23p();
    T *C33p = local->getC33p();
    T *C44p = local->getC44p();
    T *C55p = local->getC55p();
    T *C66p = local->getC66p();
    T *Rx = local->getRx();
    T *Ry = local->getRy();
    T *Rz = local->getRz();

    /* Allocate two traces to read models from file */
    T *c11trace = (T *) calloc(nx*ny, sizeof(T));
    if(c11trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c12trace = (T *) calloc(nx*ny, sizeof(T));
    if(c12trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c13trace = (T *) calloc(nx*ny, sizeof(T));
    if(c13trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c22trace = (T *) calloc(nx*ny, sizeof(T));
    if(c22trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c23trace = (T *) calloc(nx*ny, sizeof(T));
    if(c23trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c33trace = (T *) calloc(nx*ny, sizeof(T));
    if(c33trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c44trace = (T *) calloc(nx*ny, sizeof(T));
    if(c44trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c44trace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(c44trace_adv == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c55trace = (T *) calloc(nx*ny, sizeof(T));
    if(c55trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c55trace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(c55trace_adv == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *c66trace = (T *) calloc(nx*ny, sizeof(T));
    if(c66trace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");
    T *rhotrace_adv = (T *) calloc(nx*ny, sizeof(T));
    if(rhotrace_adv == NULL) rs_error("ModelOrtho3D::getDomainmodel: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> Fc11 (new File());
    status = Fc11->input(C11file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C11 file.");
    }
    std::shared_ptr<File> Fc12 (new File());
    status = Fc12->input(C12file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C12 file.");
    }
    std::shared_ptr<File> Fc13 (new File());
    status = Fc13->input(C13file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C13 file.");
    }
    std::shared_ptr<File> Fc22 (new File());
    status = Fc22->input(C22file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C22 file.");
    }
    std::shared_ptr<File> Fc23 (new File());
    status = Fc23->input(C23file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C23 file.");
    }
    std::shared_ptr<File> Fc33 (new File());
    status = Fc33->input(C33file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C33 file.");
    }
    std::shared_ptr<File> Fc44 (new File());
    status = Fc44->input(C44file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C44 file.");
    }
    std::shared_ptr<File> Fc55 (new File());
    status = Fc55->input(C55file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C55 file.");
    }
    std::shared_ptr<File> Fc66 (new File());
    status = Fc66->input(C66file);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from C66 file.");
    }
    std::shared_ptr<File> Frho (new File());
    status = Frho->input(Rfile);
    if(status == FILE_ERR){
	    rs_error("ModelOrtho3D::getDomainmodel : Error reading from Density file.");
    }

    off_t i = start_x;
    off_t j = start_y;
    off_t lpos_x, lpos_y, lpos_z, fpos;
    Index l3d(nxd, nyd, nzd);
    Index f3d(nx, ny, nz);
    Index l2d(nx, ny);
    for(size_t i1=0; i1<nzd; i1++) {
        lpos_z = iz0 + i1 - lpml;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);
        Fc11->read(c11trace, nx*ny, fpos);
        if(Fc11->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c11 file");
        Fc12->read(c12trace, nx*ny, fpos);
        if(Fc12->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c12 file");
        Fc13->read(c13trace, nx*ny, fpos);
        if(Fc13->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c13 file");
        Fc22->read(c22trace, nx*ny, fpos);
        if(Fc22->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c22 file");
        Fc23->read(c23trace, nx*ny, fpos);
        if(Fc23->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c23 file");
        Fc33->read(c33trace, nx*ny, fpos);
        if(Fc33->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c33 file");
        Fc44->read(c44trace, nx*ny, fpos);
        if(Fc44->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c44 file");
        Fc55->read(c55trace, nx*ny, fpos);
        if(Fc55->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c55 file");
        Fc66->read(c66trace, nx*ny, fpos);
        if(Fc66->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c66 file");
        Frho->read(rhotrace, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from rho file");

        lpos_z = iz0 + i1 - lpml + 1;
        if(lpos_z < 0) lpos_z = 0;
        if(lpos_z > (nz-1)) lpos_z = nz - 1;
        fpos = f3d(0, 0, lpos_z)*sizeof(T);

        Fc44->read(c44trace_adv, nx*ny, fpos);
        Fc55->read(c55trace_adv, nx*ny, fpos);
        if(Fc55->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from c55 file");
        Frho->read(rhotrace_adv, nx*ny, fpos);
        if(Frho->getFail()) rs_error("ModelOrtho3D::getDomainmodel: Error reading from rho file");

        T M1, M2, M3, M4, M5, M6, M7;
        for(size_t i3=0; i3<nyd; i3++) {
            lpos_y = j + i3 + iy0 - lpml;
                if(lpos_y < 0) lpos_y = 0;
                if(lpos_y > (ny-1)) lpos_y = ny - 1;
            for(size_t i2=0; i2<nxd; i2++) {
                lpos_x = i + i2 + ix0 - lpml;
                if(lpos_x < 0) lpos_x = 0;
                if(lpos_x > (nx-1)) lpos_x = nx - 1;
                C11[l3d(i2,i3,i1)] = c11trace[l2d(lpos_x, lpos_y)];
                C12[l3d(i2,i3,i1)] = c12trace[l2d(lpos_x, lpos_y)];
                C13[l3d(i2,i3,i1)] = c13trace[l2d(lpos_x, lpos_y)];
                C22[l3d(i2,i3,i1)] = c22trace[l2d(lpos_x, lpos_y)];
                C23[l3d(i2,i3,i1)] = c23trace[l2d(lpos_x, lpos_y)];
                C33[l3d(i2,i3,i1)] = c33trace[l2d(lpos_x, lpos_y)];
                C44[l3d(i2,i3,i1)] = c44trace[l2d(lpos_x, lpos_y)];
                C55[l3d(i2,i3,i1)] = c55trace[l2d(lpos_x, lpos_y)];
                C66[l3d(i2,i3,i1)] = c66trace[l2d(lpos_x, lpos_y)];
                R[l3d(i2,i3,i1)] = rhotrace[l2d(lpos_x, lpos_y)];
                C11p[l3d(i2,i3,i1)] = c11trace[l2d(lpos_x, lpos_y)];
                C12p[l3d(i2,i3,i1)] = c12trace[l2d(lpos_x, lpos_y)];
                C13p[l3d(i2,i3,i1)] = c13trace[l2d(lpos_x, lpos_y)];
                C22p[l3d(i2,i3,i1)] = c22trace[l2d(lpos_x, lpos_y)];
                C23p[l3d(i2,i3,i1)] = c23trace[l2d(lpos_x, lpos_y)];
                C33p[l3d(i2,i3,i1)] = c33trace[l2d(lpos_x, lpos_y)];

                if(rhotrace[l2d(lpos_x, lpos_y)] <= 0.0) rs_error("ModelOrtho3D::getDomainmodel: Zero density found.");

                if(lpos_x < nx-1){
                    Rx[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x+1, lpos_y)]);
                }else{
                    Rx[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                if(lpos_y < ny-1){
                    Ry[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace[l2d(lpos_x, lpos_y+1)]);
                }else{
                    Ry[l3d(i2,i3,i1)] = 1.0/(rhotrace[l2d(lpos_x, lpos_y)]);
                }
                Rz[l3d(i2,i3,i1)] = 2.0/(rhotrace[l2d(lpos_x, lpos_y)] + rhotrace_adv[l2d(lpos_x, lpos_y)]);

                M1 = c55trace[l2d(lpos_x, lpos_y)];
                M5 = c55trace_adv[l2d(lpos_x, lpos_y)];
                if(lpos_x < nx-1){
                   M2 = c55trace[l2d(lpos_x+1, lpos_y)];
                   M6 = c55trace_adv[l2d(lpos_x+1, lpos_y)];
                }else{
                   M2 = c55trace[l2d(lpos_x, lpos_y)];
                   M6 = c55trace_adv[l2d(lpos_x, lpos_y)];
                }
                C55p[l3d(i2,i3,i1)] = 0.25*(M1+M2+M5+M6);

                M1 = c44trace[l2d(lpos_x, lpos_y)];
                M5 = c44trace_adv[l2d(lpos_x, lpos_y)];
                if(lpos_y < ny-1){
                    M3 = c44trace[l2d(lpos_x, lpos_y+1)];
                    M7 = c44trace_adv[l2d(lpos_x, lpos_y+1)];
                }else{
                    M3 = c44trace[l2d(lpos_x, lpos_y)];
                    M7 = c44trace_adv[l2d(lpos_x, lpos_y)];
                }
                C44p[l3d(i2,i3,i1)] = 0.25*(M1+M3+M5+M7);

                M1 = c66trace[l2d(lpos_x, lpos_y)];
                if(lpos_x < nx-1){
                    M2 = c66trace[l2d(lpos_x+1, lpos_y)];
                    if(lpos_y < ny-1){
                        M4 = c66trace[l2d(lpos_x+1, lpos_y+1)];
                    }else{
                        M4 = c66trace[l2d(lpos_x+1, lpos_y)];
                    }
                }else{
                    M2 = c66trace[l2d(lpos_x, lpos_y)];
                    if(lpos_y < ny-1){
                        M4 = c66trace[l2d(lpos_x, lpos_y+1)];
                    }else{
                        M4 = c66trace[l2d(lpos_x, lpos_y)];
                    }
                }
                if(lpos_y < ny-1){
                    M3 = c66trace[l2d(lpos_x, lpos_y+1)];
                }else{
                    M3 = c66trace[l2d(lpos_x, lpos_y)];
                }
                C66p[l3d(i2,i3,i1)] = 0.25*(M1+M2+M3+M4);
            }
        }
    }

    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            for(size_t iy=0; iy<nyd; iy++){
                Rx[l3d(ix,iy,lpml-iz0)] *= 2.0;
                Ry[l3d(ix,iy,lpml-iz0)] *= 2.0;
                C11p[l3d(ix,iy,lpml-iz0)] = C12p[l3d(ix,iy,lpml-iz0)]+C66p[l3d(ix,iy,lpml-iz0)];
                C22p[l3d(ix,iy,lpml-iz0)] = C12p[l3d(ix,iy,lpml-iz0)]+C66p[l3d(ix,iy,lpml-iz0)];
                C66p[l3d(ix,iy,lpml-iz0)] *= 0.5;
            }
        }
    }

    /* Free traces */
    free(c11trace);
    free(c12trace);
    free(c13trace);
    free(c22trace);
    free(c23trace);
    free(c33trace);
    free(c44trace);
    free(c44trace_adv);
    free(c55trace);
    free(c55trace_adv);
    free(c66trace);
    free(rhotrace);
    free(rhotrace_adv);

    return local;
}

template<typename T>
ModelOrtho3D<T>::~ModelOrtho3D() {
    free(C11);
    free(C12);
    free(C13);
    free(C22);
    free(C23);
    free(C33);
    free(C44);
    free(C55);
    free(C66);
    free(C11p);
    free(C12p);
    free(C13p);
    free(C22p);
    free(C23p);
    free(C33p);
    free(C44p);
    free(C55p);
    free(C66p);
    free(R);
    free(Rx);
    free(Ry);
    free(Rz);
}

// =============== 2D POROELASTIC MODEL CLASS =============== //
template<typename T>
ModelPoroelastic2D<T>::ModelPoroelastic2D(): Model<T>(2) {
    // Nothing here
}

template<typename T>
ModelPoroelastic2D<T>::ModelPoroelastic2D(const int _nx, const int _nz, const int _lpml, const T _dx, const T _dz, const T _ox, const T _oz, const T _f0, const bool _fs): Model<T>(2, _nx, 1, _nz,  _lpml, _dx, 1.0, _dz, _ox, 1.0, _oz, _fs) {
    
    /* Allocate variables */
    Rho = (T *) calloc(1,1);
    Rhof = (T *) calloc(1,1);
    Por = (T *) calloc(1,1);
    Kd = (T *) calloc(1,1);
    Ks = (T *) calloc(1,1);
    Kf = (T *) calloc(1,1);
    Mu = (T *) calloc(1,1);
    Mob = (T *) calloc(1,1);
    Psi = (T *) calloc(1,1);
    Lu = (T *) calloc(1,1);
    LuM = (T *) calloc(1,1);
    Alpha = (T *) calloc(1,1);
    Beta = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    Rho_x = (T *) calloc(1,1);
    Rho_z = (T *) calloc(1,1);
    Rhof_x = (T *) calloc(1,1);
    Rhof_z = (T *) calloc(1,1);
    Mob_x = (T *) calloc(1,1);
    Mob_z = (T *) calloc(1,1);
    Psi_x = (T *) calloc(1,1);
    Psi_z = (T *) calloc(1,1);
    f0 = _f0;
}

template<typename T>
ModelPoroelastic2D<T>::ModelPoroelastic2D(std::string _Rhofile, std::string _Rhoffile, std::string _Porfile, std::string _Kdfile, std::string _Ksfile,  std::string _Kffile, std::string _Mufile, std::string _Mobfile, std::string _Psifile, const int _lpml, const T _f0, const bool _fs): Model<T>(2) {
    bool status;
    int nx, nz;
    T dx, dz;
    T ox, oz;
    Rhofile = _Rhofile;
    Rhoffile = _Rhoffile;
    Porfile = _Porfile;
    Kdfile = _Kdfile;
    Ksfile = _Ksfile;
    Kffile = _Kffile;
    Mufile = _Mufile;
    Mobfile = _Mobfile;
    Psifile = _Psifile;

    std::shared_ptr<File> FRho (new File());
    status = FRho->input(Rhofile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Rho file: ", Rhofile);
	    exit(1);
    }
    std::shared_ptr<File> FRhof (new File());
    status = FRhof->input(Rhoffile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Rhof file: ", Rhoffile);
	    exit(1);
    }
    std::shared_ptr<File> FPor (new File());
    status = FPor->input(Porfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Por file: ", Porfile);
	    exit(1);
    }
    std::shared_ptr<File> FKd (new File());
    status = FKd->input(Kdfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Kd file: ", Kdfile);
	    exit(1);
    }
    std::shared_ptr<File> FKs (new File());
    status = FKs->input(Ksfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Ks file: ", Ksfile);
	    exit(1);
    }
    std::shared_ptr<File> FKf (new File());
    status = FKf->input(Kffile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Kf file: ", Kffile);
	    exit(1);
    }
    std::shared_ptr<File> FMu (new File());
    status = FMu->input(Mufile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Mu file: ", Mufile);
	    exit(1);
    }
    std::shared_ptr<File> FMob (new File());
    status = FMob->input(Mobfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Mob file: ", Mobfile);
	    exit(1);
    }
    std::shared_ptr<File> FPsi (new File());
    status = FPsi->input(Psifile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::Error reading from Psi file: ", Psifile);
	    exit(1);
    }

    // Compare geometry in the two files
    if(FRho->compareGeometry(FRhof) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Rhof model files do not match.");
    }
    if(FRho->compareGeometry(FPor) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Por model files do not match.");
    }
    if(FRho->compareGeometry(FKd) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Kd model files do not match.");
    }
    if(FRho->compareGeometry(FKs) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Ks model files do not match.");
    }
    if(FRho->compareGeometry(FKf) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Kf model files do not match.");
    }
    if(FRho->compareGeometry(FMu) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Mu model files do not match.");
    }
    if(FRho->compareGeometry(FMob) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Mob model files do not match.");
    }
    if(FRho->compareGeometry(FPsi) != 0)
    {
        rs_error("ModelPoroelastic2D::Geometries in Rho and Psi model files do not match.");
    }
    if(FRho->getData_format() != sizeof(T))
    {
        rs_error("ModelPoroelastic2D::Numerical precision in Rho mismatch with constructor.");
    }
 

    
    // Read geometry from file
    nx = FRho->getN(1);
    dx = (T) FRho->getD(1);
    ox = (T) FRho->getO(1);
    nz = FRho->getN(3);
    dz = (T) FRho->getD(3);
    oz = (T) FRho->getO(3);
    
    // Close files
    FRho->close();
    FRhof->close();
    FPor->close();
    FKd->close();
    FKs->close();
    FKf->close();
    FMu->close();
    FMob->close();
    FPsi->close();

    // Store geometry in model class
    this->setNx(nx);
    this->setDx(dx);
    this->setOx(ox);

    this->setNz(nz);
    this->setDz(dz);
    this->setOz(oz);

    this->setLpml(_lpml);
    this->setFs(_fs);
    this->setF0(_f0);

    
    /* Allocate variables */
    Rho = (T *) calloc(1,1);
    Rhof = (T *) calloc(1,1);
    Por = (T *) calloc(1,1);
    Kd = (T *) calloc(1,1);
    Ks = (T *) calloc(1,1);
    Kf = (T *) calloc(1,1);
    Mu = (T *) calloc(1,1);
    Mob = (T *) calloc(1,1);
    Psi = (T *) calloc(1,1);
    Lu = (T *) calloc(1,1);
    LuM = (T *) calloc(1,1);
    Alpha = (T *) calloc(1,1);
    Beta = (T *) calloc(1,1);
    M_xz = (T *) calloc(1,1);
    Rho_x = (T *) calloc(1,1);
    Rho_z = (T *) calloc(1,1);
    Rhof_x = (T *) calloc(1,1);
    Rhof_z = (T *) calloc(1,1);
    Mob_x = (T *) calloc(1,1);
    Mob_z = (T *) calloc(1,1);
    Psi_x = (T *) calloc(1,1);
    Psi_z = (T *) calloc(1,1);
}

template<typename T>
void ModelPoroelastic2D<T>::readModel() {
    bool status;
    // Get file names
    std::string Rhofile = this->getRhofile();
    std::string Rhoffile = this->getRhoffile();
    std::string Porfile = this->getPorfile();
    std::string Kdfile = this->getKdfile();
    std::string Ksfile = this->getKsfile();
    std::string Kffile = this->getKffile();
    std::string Mufile = this->getMufile();
    std::string Mobfile = this->getMobfile();
    std::string Psifile = this->getPsifile();
    // Open files for reading
    std::shared_ptr<File> FRho (new File());
    status = FRho->input(Rhofile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Rho file: ", Rhofile);
    }
    std::shared_ptr<File> FRhof (new File());
    status = FRhof->input(Rhoffile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Rhof file: ", Rhoffile);
    }

    std::shared_ptr<File> FPor (new File());
    status = FPor->input(Porfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Por file: ", Porfile);
    }

    std::shared_ptr<File> FKd (new File());
    status = FKd->input(Kdfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Kd file: ", Kdfile);
    }

    std::shared_ptr<File> FKs (new File());
    status = FKs->input(Ksfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Ks file: ", Ksfile);
    }

    std::shared_ptr<File> FKf (new File());
    status = FKf->input(Kffile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Kf file: ", Kffile);
    }

    std::shared_ptr<File> FMu (new File());
    status = FMu->input(Mufile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Mu file: ", Mufile);
    }

    std::shared_ptr<File> FMob (new File());
    status = FMob->input(Mobfile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Mob file: ", Mobfile);
    }

    std::shared_ptr<File> FPsi (new File());
    status = FPsi->input(Psifile.c_str());
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::readModel : Error reading from Psi file: ", Psifile);
    }


    // Read models
    int nx = this->getNx();
    int nz = this->getNz();
    
    /* Reallocate Models */
    free(Rho); 
    Rho = (T *) calloc(nx*nz,sizeof(T));
    if(Rho == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Rhof); 
    Rhof = (T *) calloc(nx*nz,sizeof(T));
    if(Rhof == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Por); 
    Por = (T *) calloc(nx*nz,sizeof(T));
    if(Por == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Kd); 
    Kd = (T *) calloc(nx*nz,sizeof(T));
    if(Kd == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Ks); 
    Ks = (T *) calloc(nx*nz,sizeof(T));
    if(Ks == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Kf); 
    Kf = (T *) calloc(nx*nz,sizeof(T));
    if(Kf == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Mu); 
    Mu = (T *) calloc(nx*nz,sizeof(T));
    if(Mu == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Mob); 
    Mob = (T *) calloc(nx*nz,sizeof(T));
    if(Mob == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    free(Psi); 
    Psi = (T *) calloc(nx*nz,sizeof(T));
    if(Psi == NULL) rs_error("ModelPoroelastic2D::readModel: Failed to allocate memory.");
    this->setRealized(true);

    FRho->read(Rho, nx*nz);
    FRho->close();
    FRhof->read(Rhof, nx*nz);
    FRhof->close();
    FPor->read(Por, nx*nz);
    FPor->close();
    FKd->read(Kd, nx*nz);
    FKd->close();
    FKs->read(Ks, nx*nz);
    FKs->close();
    FKf->read(Kf, nx*nz);
    FKf->close();
    FMu->read(Mu, nx*nz);
    FMu->close();
    FMob->read(Mob, nx*nz);
    FMob->close();
    FPsi->read(Psi, nx*nz);
    FPsi->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeRho() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeRho: Model is not allocated.");
    }
    // Get file names
    std::string Rhofile = this->getRhofile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Rhofile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getRho();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeRhof() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeRhof: Model is not allocated.");
    }
    // Get file names
    std::string Rhoffile = this->getRhoffile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Rhoffile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getRhof();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writePor() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writePor: Model is not allocated.");
    }
    // Get file names
    std::string Porfile = this->getPorfile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Porfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getPor();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeKd() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeKd: Model is not allocated.");
    }
    // Get file names
    std::string Kdfile = this->getKdfile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Kdfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getKd();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeKs() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeKs: Model is not allocated.");
    }
    // Get file names
    std::string Ksfile = this->getKsfile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Ksfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getKs();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeKf() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeKf: Model is not allocated.");
    }
    // Get file names
    std::string Kffile = this->getKffile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Kffile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getKf();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeMu() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeMu: Model is not allocated.");
    }
    // Get file names
    std::string Mufile = this->getMufile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Mufile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getMu();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writeMob() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writeMob: Model is not allocated.");
    }
    // Get file names
    std::string Mobfile = this->getMobfile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Mobfile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getMob();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::writePsi() {
    if(!this->getRealized()) {
        rs_error("ModelPoroelastic2D::writePsi: Model is not allocated.");
    }
    // Get file names
    std::string Psifile = this->getPsifile();
    // Open files for writting
    std::shared_ptr<File> F (new File());
    F->output(Psifile);

    // Write models
    int nx = this->getNx();
    T dx = this->getDx();
    T ox = this->getOx();
    int nz = this->getNz();
    T dz = this->getDz();
    T oz = this->getOz();

    F->setN(1,nx);
    F->setN(3,nz);
    F->setD(1,dx);
    F->setD(3,dz);
    F->setO(1,ox);
    F->setO(3,oz);
    F->setType(REGULAR);
    F->setData_format(sizeof(T));
    F->writeHeader();
    T *Mod = this->getPsi();
    F->write(Mod, nx*nz, 0);
    F->close();
}

template<typename T>
void ModelPoroelastic2D<T>::createModel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate models */
    free(Rho); 
    Rho = (T *) calloc(nx*nz,sizeof(T));
    if(Rho == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Rhof); 
    Rhof = (T *) calloc(nx*nz,sizeof(T));
    if(Rhof == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Por); 
    Por = (T *) calloc(nx*nz,sizeof(T));
    if(Por == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Kd); 
    Kd = (T *) calloc(nx*nz,sizeof(T));
    if(Kd == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Ks); 
    Ks = (T *) calloc(nx*nz,sizeof(T));
    if(Ks == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Kf); 
    Kf = (T *) calloc(nx*nz,sizeof(T));
    if(Kf == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Mu); 
    Mu = (T *) calloc(nx*nz,sizeof(T));
    if(Mu == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Mob); 
    Mob = (T *) calloc(nx*nz,sizeof(T));
    if(Mob == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");
    free(Psi); 
    Psi = (T *) calloc(nx*nz,sizeof(T));
    if(Psi == NULL) rs_error("ModelPoroelastic2D::createModel: Failed to allocate memory.");

    this->setRealized(true);
}

template<typename T>
void ModelPoroelastic2D<T>::createPaddedmodel() {
    int nx = this->getNx();
    int nz = this->getNz();

    /* Reallocate models */
    free(Lu); 
    Lu = (T *) calloc(nx*nz,sizeof(T));
    if(Lu == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");
    free(LuM); 
    LuM = (T *) calloc(nx*nz,sizeof(T));
    if(LuM == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Alpha); 
    Alpha = (T *) calloc(nx*nz,sizeof(T));
    if(Alpha == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Beta); 
    Beta = (T *) calloc(nx*nz,sizeof(T));
    if(Beta == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(M_xz); 
    M_xz = (T *) calloc(nx*nz,sizeof(T));
    if(M_xz == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Rho_x); 
    Rho_x = (T *) calloc(nx*nz,sizeof(T));
    if(Rho_x == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Rho_z); 
    Rho_z = (T *) calloc(nx*nz,sizeof(T));
    if(Rho_z == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Rhof_x); 
    Rhof_x = (T *) calloc(nx*nz,sizeof(T));
    if(Rhof_x == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Rhof_z); 
    Rhof_z = (T *) calloc(nx*nz,sizeof(T));
    if(Rhof_z == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Mob_x); 
    Mob_x = (T *) calloc(nx*nz,sizeof(T));
    if(Mob_x == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Mob_z); 
    Mob_z = (T *) calloc(nx*nz,sizeof(T));
    if(Mob_z == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Psi_x); 
    Psi_x = (T *) calloc(nx*nz,sizeof(T));
    if(Psi_x == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");

    free(Psi_z); 
    Psi_z = (T *) calloc(nx*nz,sizeof(T));
    if(Psi_z == NULL) rs_error("ModelPoroelastic2D::createPaddedmodel: Failed to allocate memory.");
}

template<typename T>
std::shared_ptr<ModelPoroelastic2D<T>> ModelPoroelastic2D<T>::getLocal(std::shared_ptr<Data2D<T>> data, T aperture, bool map) {

    std::shared_ptr<ModelPoroelastic2D<T>> local;
    T dx = this->getDx();
    T ox = this->getOx();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);

    /* Create local model */
    local = std::make_shared<ModelPoroelastic2D<T>>(size, this->getNz(), this->getLpml(), dx, this->getDz(), (ox + start*dx) , this->getOz(), this->getF0(), this->getFs());

    /*Realizing local model */
    local->createModel();

	/* Copying from big model into local model */
    T *Rho = local->getRho();
    T *Rhof = local->getRhof();
    T *Por = local->getPor();
    T *Kd = local->getKd();
    T *Ks = local->getKs();
    T *Kf = local->getKf();
    T *Mu = local->getMu();
    T *Mob = local->getMob();
    T *Psi = local->getPsi();

    /* Allocate two traces to read models from file */
    T *Rhotrace = (T *) calloc(nx, sizeof(T));
    if(Rhotrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Rhoftrace = (T *) calloc(nx, sizeof(T));
    if(Rhoftrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Portrace = (T *) calloc(nx, sizeof(T));
    if(Portrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Kdtrace = (T *) calloc(nx, sizeof(T));
    if(Kdtrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Kstrace = (T *) calloc(nx, sizeof(T));
    if(Kstrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Kftrace = (T *) calloc(nx, sizeof(T));
    if(Kftrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Mutrace = (T *) calloc(nx, sizeof(T));
    if(Mutrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Mobtrace = (T *) calloc(nx, sizeof(T));
    if(Mobtrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Psitrace = (T *) calloc(nx, sizeof(T));
    if(Psitrace == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    // Open files for reading
    bool status;
    std::shared_ptr<File> FRho (new File());
    status = FRho->input(Rhofile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Rho file.");
    }
    std::shared_ptr<File> FRhof (new File());
    status = FRhof->input(Rhoffile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Rhof file.");
    }
    std::shared_ptr<File> FPor (new File());
    status = FPor->input(Porfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Por file.");
    }
    std::shared_ptr<File> FKd (new File());
    status = FKd->input(Kdfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Kd file.");
    }
    std::shared_ptr<File> FKs (new File());
    status = FKs->input(Ksfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Ks file.");
    }
    std::shared_ptr<File> FKf (new File());
    status = FKf->input(Kffile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Kf file.");
    }
    std::shared_ptr<File> FMu (new File());
    status = FMu->input(Mufile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Mu file.");
    }
    std::shared_ptr<File> FMob (new File());
    status = FMob->input(Mobfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Mob file.");
    }
    std::shared_ptr<File> FPsi (new File());
    status = FPsi->input(Psifile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getLocal : Error reading from Psi file.");
    }

    off_t i = start;
    off_t lpos, fpos;
    Index l2d(size,nz);
    Index f2d(nx,nz);
    for(size_t i1=0; i1<nz; i1++) {
        fpos = f2d(0, i1)*sizeof(T);
        FRho->read(Rhotrace, nx, fpos);
        if(FRho->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Rho file");
        FRhof->read(Rhoftrace, nx, fpos);
        if(FRhof->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Rhof file");
        FPor->read(Portrace, nx, fpos);
        if(FPor->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Por file");
        FKd->read(Kdtrace, nx, fpos);
        if(FKd->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Kd file");
        FKs->read(Kstrace, nx, fpos);
        if(FKs->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Ks file");
        FKf->read(Kftrace, nx, fpos);
        if(FKf->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Kf file");
        FMu->read(Mutrace, nx, fpos);
        if(FMu->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Mu file");
        FMob->read(Mobtrace, nx, fpos);
        if(FMob->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Mob file");
        FPsi->read(Psitrace, nx, fpos);
        if(FPsi->getFail()) rs_error("ModelPoroelastic2D::getLocal: Error reading from Psi file");
        for(size_t i2=0; i2<size; i2++) {
            lpos = i + i2;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Rho[l2d(i2,i1)] = Rhotrace[lpos];
            Rhof[l2d(i2,i1)] = Rhoftrace[lpos];
            Por[l2d(i2,i1)] = Portrace[lpos];
            Kd[l2d(i2,i1)] = Kdtrace[lpos];
            Ks[l2d(i2,i1)] = Kstrace[lpos];
            Kf[l2d(i2,i1)] = Kftrace[lpos];
            Mu[l2d(i2,i1)] = Mutrace[lpos];
            Mob[l2d(i2,i1)] = Mobtrace[lpos];
            Psi[l2d(i2,i1)] = Psitrace[lpos];
        }
    }

    /* Free traces */
    free(Rhotrace);
    free(Rhoftrace);
    free(Portrace);
    free(Kdtrace);
    free(Kstrace);
    free(Kftrace);
    free(Mutrace);
    free(Mobtrace);
    free(Psitrace);

    return local;
}

template<typename T>
std::shared_ptr<ModelPoroelastic2D<T>> ModelPoroelastic2D<T>::getDomainmodel(std::shared_ptr<Data2D<T>> data, T aperture, bool map, const int d, const int nd0, const int nd1, const int order) {
    std::shared_ptr<ModelPoroelastic2D<T>> local;
    T dx = this->getDx();
    T dz = this->getDz();
    T ox = this->getOx();
    T oz = this->getOz();
    size_t nz = this->getNz();
    size_t nx = this->getNx();
    size_t size;
    off_t start;
    int nxd,nzd;
    int ix0,iz0;
    int lpml = this->getLpml();

    /* Determine grid positions and sizes */
    this->getLocalsize2d(data, aperture, map, &start, &size);
    (this->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*this->getLpml(),d,nd0,1,nd1,order);
    nxd = (this->getDomain())->getNx_pad();
    nzd = (this->getDomain())->getNz_pad();
    ix0 = (this->getDomain())->getIx0();
    iz0 = (this->getDomain())->getIz0();


    /* Create domain model */
    local = std::make_shared<ModelPoroelastic2D<T>>(nxd, nzd, lpml, dx, dz, (ox + (start+ix0-lpml)*dx) , (oz + (iz0-lpml)*dz), this->getF0(), this->getFs());
    (local->getDomain())->setupDomain3D(size+2*lpml,1,nz+2*lpml,d,nd0,1,nd1,order);

    /*Realizing local model */
    local->createModel();
    local->createPaddedmodel();

    /* Copying from big model into local model */
    T *Rho = local->getRho();
    T *Rhof = local->getRhof();
    T *Por = local->getPor();
    T *Kd = local->getKd();
    T *Ks = local->getKs();
    T *Kf = local->getKf();
    T *Mu = local->getMu();
    T *Mob = local->getMob();
    T *Psi = local->getPsi();
    T *Lu = local->getLu();
    T *LuM = local->getLuM();
    T *Alpha = local->getAlpha();
    T *Beta = local->getBeta();
    T *M_xz = local->getM_xz();
    T *Rho_x = local->getRho_x();
    T *Rho_z = local->getRho_z();
    T *Rhof_x = local->getRhof_x();
    T *Rhof_z = local->getRhof_z();

    T *Mob_x = local->getMob_x();
    T *Mob_z = local->getMob_z();
    T *Psi_x = local->getPsi_x();
    T *Psi_z = local->getPsi_z();

    /* Allocate two traces to read models from file */
    T *Rhotrace = (T *) calloc(nx, sizeof(T));
    if(Rhotrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Rhoftrace = (T *) calloc(nx, sizeof(T));
    if(Rhoftrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Portrace = (T *) calloc(nx, sizeof(T));
    if(Portrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Kdtrace = (T *) calloc(nx, sizeof(T));
    if(Kdtrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Kstrace = (T *) calloc(nx, sizeof(T));
    if(Kstrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Kftrace = (T *) calloc(nx, sizeof(T));
    if(Kftrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Mutrace = (T *) calloc(nx, sizeof(T));
    if(Mutrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Mobtrace = (T *) calloc(nx, sizeof(T));
    if(Mobtrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Psitrace = (T *) calloc(nx, sizeof(T));
    if(Psitrace == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Rhotrace_adv = (T *) calloc(nx, sizeof(T));
    if(Rhotrace_adv == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Rhoftrace_adv = (T *) calloc(nx, sizeof(T));
    if(Rhoftrace_adv == NULL) rs_error("ModelPoroelastic2d::getDomainmodel: Failed to allocate memory.");
    T *Mutrace_adv = (T *) calloc(nx, sizeof(T));
    if(Mutrace_adv == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Mobtrace_adv = (T *) calloc(nx, sizeof(T));
    if(Mobtrace_adv == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");
    T *Psitrace_adv = (T *) calloc(nx, sizeof(T));
    if(Psitrace_adv == NULL) rs_error("ModelPoroelastic2d::getLocal: Failed to allocate memory.");

    // Open files for reading
    bool status;
    std::shared_ptr<File> FRho (new File());
    status = FRho->input(Rhofile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Rho file.");
    }
    std::shared_ptr<File> FRhof (new File());
    status = FRhof->input(Rhoffile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Rhof file.");
    }
    std::shared_ptr<File> FPor (new File());
    status = FPor->input(Porfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Por file.");
    }
    std::shared_ptr<File> FKd (new File());
    status = FKd->input(Kdfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Kd file.");
    }
    std::shared_ptr<File> FKs (new File());
    status = FKs->input(Ksfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Ks file.");
    }
    std::shared_ptr<File> FKf (new File());
    status = FKf->input(Kffile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Kf file.");
    }
    std::shared_ptr<File> FMu (new File());
    status = FMu->input(Mufile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Mu file.");
    }
    std::shared_ptr<File> FMob (new File());
    status = FMob->input(Mobfile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Mob file.");
    }
    std::shared_ptr<File> FPsi (new File());
    status = FPsi->input(Psifile);
    if(status == FILE_ERR){
	    rs_error("ModelPoroelastic2D::getDomainmodel : Error reading from Psi file.");
    }
    
    off_t i = start;
    off_t lpos, fpos;
    Index l2d(nxd,nzd);
    Index f2d(nx,nz);
    T M1, M2, M3, M4;

    for(size_t i1=0; i1<nzd; i1++) {
        lpos = iz0 + i1 - lpml;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        FRho->read(Rhotrace, nx, fpos);
        if(FRho->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Rho file");
        FRhof->read(Rhoftrace, nx, fpos);
        if(FRhof->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Rhof file");
        FPor->read(Portrace, nx, fpos);
        if(FPor->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Por file");
        FKd->read(Kdtrace, nx, fpos);
        if(FKd->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Kd file");
        FKs->read(Kstrace, nx, fpos);
        if(FKs->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Ks file");
        FKf->read(Kftrace, nx, fpos);
        if(FKf->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Kf file");
        FMu->read(Mutrace, nx, fpos);
        if(FMu->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Mu file");
        FMob->read(Mobtrace, nx, fpos);
        if(FMob->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Mob file");
        FPsi->read(Psitrace, nx, fpos);
        if(FPsi->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Psi file");

        // Read advanced trace
        lpos = iz0 + i1 - lpml + 1;
        if(lpos < 0) lpos = 0;
        if(lpos > (nz-1)) lpos = nz - 1;
        fpos = f2d(0, lpos)*sizeof(T);
        FRho->read(Rhotrace_adv, nx, fpos);
        if(FRho->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Rho file");
        FRhof->read(Rhoftrace_adv, nx, fpos);
        if(FRhof->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Rhof file");
        FMu->read(Mutrace_adv, nx, fpos);
        if(FMu->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Mu file");
        FMob->read(Mobtrace_adv, nx, fpos);
        if(FMob->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Mob file");
        FPsi->read(Psitrace_adv, nx, fpos);
        if(FPsi->getFail()) rs_error("ModelPoroelastic2D::getDomainmodel: Error reading from Psi file");
        for(size_t i2=0; i2<nxd; i2++) {
            lpos = i + i2 + ix0 - lpml;
            if(lpos < 0) lpos = 0;
            if(lpos > (nx-1)) lpos = nx - 1;
            Rho[l2d(i2,i1)] = Rhotrace[lpos];
            Rhof[l2d(i2,i1)] = Rhoftrace[lpos];
            Por[l2d(i2,i1)] = Portrace[lpos];
            Kd[l2d(i2,i1)] = Kdtrace[lpos];
            Ks[l2d(i2,i1)] = Kstrace[lpos];
            Kf[l2d(i2,i1)] = Kftrace[lpos];
            Mu[l2d(i2,i1)] = Mutrace[lpos];
            Mob[l2d(i2,i1)] = Mobtrace[lpos];
            Psi[l2d(i2,i1)] = Psitrace[lpos];
            if(Mobtrace[lpos] <= 0.0) rs_error("ModelPoroelastic2D::getDomainmodel: Zero mobility found.");
            if(Kstrace[lpos] <= 0.0) rs_error("ModelPoroelastic2D::getDomainmodel: Zero Ks found.");
            if(lpos < nx-1){
                Rho_x[l2d(i2,i1)] = (Rhotrace[lpos]+Rhotrace[lpos+1])/2.0;
            }else{
                Rho_x[l2d(i2,i1)] = (Rhotrace[lpos]);
            }
            Rho_z[l2d(i2,i1)] = (Rhotrace[lpos]+Rhotrace_adv[lpos])/2.0;
            if(lpos < nx-1){
                Rhof_x[l2d(i2,i1)] = (Rhoftrace[lpos]+Rhoftrace[lpos+1])/2.0;
            }else{
                Rhof_x[l2d(i2,i1)] = (Rhoftrace[lpos]);
            }
            Rhof_z[l2d(i2,i1)] = (Rhoftrace[lpos]+Rhoftrace_adv[lpos])/2.0;

            if(lpos < nx-1){
                Mob_x[l2d(i2,i1)] = 2.0/(Mobtrace[lpos]+Mobtrace[lpos+1]);
            }else{
                Mob_x[l2d(i2,i1)] = 1.0/(Mobtrace[lpos]);
            }
            Mob_z[l2d(i2,i1)] = 2.0/(Mobtrace[lpos]+Mobtrace_adv[lpos]);
            if(lpos < nx-1){
                Psi_x[l2d(i2,i1)] = (Psitrace[lpos]+Psitrace[lpos+1])/2.0;
            }else{
                Psi_x[l2d(i2,i1)] = (Psitrace[lpos]);
            }
            Psi_z[l2d(i2,i1)] = (Psitrace[lpos]+Psitrace_adv[lpos])/2.0;
            M1 = Mutrace[lpos];
            M3 = Mutrace_adv[lpos];
            if(lpos < nx-1){
               M2 = Mutrace[lpos+1];
               M4 = Mutrace_adv[lpos+1];
            }else{
               M2 = Mutrace[lpos];
               M4 = Mutrace_adv[lpos];
            }
            M_xz[l2d(i2,i1)] = 0.25*(M1+M2+M3+M4);

            Alpha[l2d(i2,i1)] = 1.0 - (Kdtrace[lpos]/Kstrace[lpos]);
            Beta[l2d(i2,i1)] = 1.0/((Portrace[lpos]/Kftrace[lpos]) + ((Alpha[l2d(i2,i1)] - Portrace[lpos])/Kstrace[lpos]));
			Lu[l2d(i2,i1)] = Kdtrace[lpos] + (Alpha[l2d(i2,i1)]*Alpha[l2d(i2,i1)]*Beta[l2d(i2,i1)]) - ((2.0/3.0)*Mutrace[lpos]);
			LuM[l2d(i2,i1)] = Lu[l2d(i2,i1)] + 2.0*Mutrace[lpos];
        }
    }

    // In case of free surface
    if(this->getFs() && ((iz0 <= lpml) && ((iz0+nzd) >= lpml))){
        for(size_t ix=0; ix<nxd; ix++){
            Rho_x[l2d(ix,lpml-iz0)] /= 2.0;
            LuM[l2d(ix,lpml-iz0)] = (LuM[l2d(ix,lpml-iz0)] - Lu[l2d(ix,lpml-iz0)])/2.0;
            Lu[l2d(ix,lpml-iz0)] = 0.0;
            Alpha[l2d(ix,lpml-iz0)] = 0.0;
            Beta[l2d(ix,lpml-iz0)] = 0.0;
        }
    }



    /* Free traces */
    free(Rhotrace);   
    free(Rhoftrace);   
    free(Portrace);   
    free(Kdtrace);   
    free(Kstrace);   
    free(Kftrace);   
    free(Mutrace);   
    free(Mobtrace);   
    free(Psitrace);   
    free(Rhotrace_adv);
    free(Rhoftrace_adv);
    free(Mutrace_adv);   
    free(Mobtrace_adv);
    free(Psitrace_adv);
    
    return local;
}


template<typename T>
ModelPoroelastic2D<T>::~ModelPoroelastic2D() {
    free(Rho);   
    free(Rhof);   
    free(Por);   
    free(Kd);   
    free(Ks);   
    free(Kf);   
    free(Mu);   
    free(Mob);   
    free(Psi);   
    free(Lu);     
    free(LuM);    
    free(Alpha);     
    free(Beta);     
    free(M_xz);   
    free(Rho_x);     
    free(Rho_z);     
    free(Rhof_x);   
    free(Rhof_z);   
    free(Mob_x);  
    free(Mob_z);  
    free(Psi_x);  
    free(Psi_z);  
}



// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Model<float>;
template class Model<double>;

template class ModelEikonal2D<float>;
template class ModelEikonal3D<float>;
template class ModelAcoustic2D<float>;
template class ModelAcoustic3D<float>;
template class ModelElastic2D<float>;
template class ModelElastic3D<float>;
template class ModelPoroelastic2D<float>;

template class ModelEikonal2D<double>;
template class ModelEikonal3D<double>;
template class ModelAcoustic2D<double>;
template class ModelAcoustic3D<double>;
template class ModelElastic2D<double>;
template class ModelElastic3D<double>;
template class ModelPoroelastic2D<double>;

template class ModelViscoelastic2D<float>;
template class ModelViscoelastic2D<double>;

template class ModelViscoelastic3D<float>;
template class ModelViscoelastic3D<double>;

template class ModelVti2D<float>;
template class ModelVti2D<double>;

template class ModelOrtho3D<float>;
template class ModelOrtho3D<double>;
}
