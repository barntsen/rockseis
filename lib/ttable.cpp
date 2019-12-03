#include "ttable.h"

namespace rockseis {
// constructor
template<typename T>
Ttable<T>::Ttable()
{
    ptree = kd_create(3);
    allocated = false;
    geomset = false;
}

template<typename T>
Ttable<T>::Ttable(const int _dim, const int _ntable)
{
    allocated = false;
    geometry = std::make_shared<Geometry3D<T>>(_ntable); 
    geomset = true;
    ptree = kd_create(3);
    this->setDim(_dim);
    lpml = 0; 
    allocated = false;
}

template<typename T>
Ttable<T>::~Ttable() {
    // Free variables
    kd_free(ptree);
}



// =============== 2D TTABLE CLASS =============== //
template<typename T>
Ttable2D<T>::Ttable2D(std::shared_ptr<ModelEikonal2D<T>> model, int _ntable):Ttable<T>(2, _ntable)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _lpml;
    int _dim;

/* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    // TODO implement the possibility to pad the traveltime table
    //_lpml = model->getLpml();
    _lpml = 0;
    _dim = model->getDim();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox-(_lpml*_dx));
    this->setOy(_oy);
    this->setOz(_oz-(_lpml*_dz));
    this->setDim(_dim);
    this->setNtable(_ntable);
    this->setLpml(_lpml);
}

template<typename T>
Ttable2D<T>::Ttable2D(std::string tablefile)
{
    bool status;
    /* Set filename */
    this->setFilename(tablefile);
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(tablefile);
    if(status == FILE_ERR) rs_error("Ttable2D::Ttable2D: Error opening data file for reading: ", tablefile);

    std::string posfile = tablefile+"-pos";
    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    status = Fpos->input(posfile);
    if(status == FILE_ERR) rs_error("Ttable2D::Ttable2D: Error opening position file for reading: ", tablefile);

    size_t ntable = Fdata->getN(4);
    if(ntable < 1) rs_error("Number of tables in file ", tablefile,  " less than 1.");
    if(Fpos->getN(2) != ntable) rs_error("Ttable2D::Ttable2D: Number of gathers mismatch in datafile and positions file.");
    this->makeGeom(ntable);
    this->setNtable(ntable);

    /* Get necessary parameters from file */
    this->setNx(Fdata->getN(1));
    this->setNy(Fdata->getN(2));
    this->setNz(Fdata->getN(3));
    this->setDx(Fdata->getD(1));
    this->setDy(Fdata->getD(2));
    this->setDz(Fdata->getD(3));
    this->setOx(Fdata->getO(1));
    this->setOy(Fdata->getO(2));
    this->setOz(Fdata->getO(3));
    Fdata->close();
    this->setLpml(0);

    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;
    if(n1 == 0) rs_error("Ttable::Ttable2D: Zero dimension array. Check Nx, Ny, Nz values.");
    if(this->getNy() > 1) {
        rs_error("Ttable::Ttable2D: Traveltime table is 3 dimensional.");
    }

    kdtree *ptree = this->getKdtree();
    // Read coordinates and populate tree
    double pos[3];
    T val;
    for (size_t i=0; i < ntable; i++){
        Fpos->read(&val, 1);
        pos[0] = (double) val;
        Fpos->read(&val, 1);
        pos[1] = (double) val;
        Fpos->read(&val, 1);
        pos[2] = (double) val;
        if(kd_insert(ptree, &pos[0], i) != 0){
            rs_error("Ttable::Ttable2D: Error populating tree.");
        }

    }
    Fpos->close();
}

template<typename T>
void Ttable2D<T>::insertSource(std::shared_ptr<rockseis::Data2D<T>> source, bool maptype, int traceno){
    Point2D<T> *pos;
    Point3D<T> *posi;
    size_t ntrace = source->getNtrace();

    if(traceno < 0 || traceno > ntrace-1){
        rs_error("Ttable2D<T>::insertSource: traceno out of bounds.");
    }

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        pos = (source->getGeom())->getScoords();
    }else{
        pos = (source->getGeom())->getGcoords();
    }

    posi = (this->getGeom())->getScoords();
    posi[0].x = pos[traceno].x;
    posi[0].y = 0.0;
    posi[0].z = pos[traceno].y;
}

template<typename T>
void Ttable2D<T>::allocTtable() 
{
    if(this->getAllocated() == false) {
        this->data = (T *) calloc(this->getNx_pml()*this->getNz_pml(), sizeof(T));
        this->wrk = (T *) calloc(this->getNx_pml()*this->getNz_pml(), sizeof(T));
        this->setAllocated(true);
    }else{
        free(this->data); free(this->wrk);
        this->data = (T *) calloc(this->getNx_pml()*this->getNz_pml(), sizeof(T));
        this->wrk = (T *) calloc(this->getNx_pml()*this->getNz_pml(), sizeof(T));
        this->setAllocated(true);
    }
}

template<typename T>
void Ttable2D<T>::createEmptyttable() 
{
    int ntable = this->getNtable();
    if(ntable < 1) rs_error("Ttable2D::createEmptyTtable2D: number of tables is less than 1.");

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable2D::createEmptyTtable2D: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    Fdata->output(datafile);

    std::string posfile = datafile+"-pos";

    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    Fpos->output(posfile);

    size_t nx = this->getNx_pml();
    size_t ny = 1;
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;
    if(n1 == 0) rs_error("Ttable::createEmptydata: Zero dimension array. Check Nx, Ny, Nz values.");
    Fdata->setN(1, nx);
    Fdata->setN(2, ny);
    Fdata->setN(3, nz);
    Fdata->setN(4,ntable);
    Fdata->setD(1,this->getDx());
    Fdata->setD(2,1.0);
    Fdata->setD(3,this->getDz());
    Fdata->setD(4,1);
    Fdata->setO(1,this->getOx());
    Fdata->setO(2,0.0);
    Fdata->setO(3,this->getOz());
    Fdata->setO(4,0);
    Fdata->setData_format(sizeof(T));
    Fdata->setType(rockseis::REGULAR);
    Fdata->writeHeader();
    Fdata->seekp(Fdata->getStartofdata());

    Fpos->setN(1,3);
    Fpos->setD(1,1.0);
    Fpos->setN(2, ntable);
    Fpos->setD(2,1.0);
    Fpos->setData_format(sizeof(T));
    Fpos->setType(rockseis::REGULAR);
    Fpos->writeHeader();
    Fpos->seekp(Fpos->getStartofdata());

    this->allocTtable();
    T val=0.0;
    for (int i=0; i < ntable; i++){
        Fpos->write(&val, 1);
        Fpos->write(&val, 1);
        Fpos->write(&val, 1);
        Fdata->write(this->getData(),n1);
    }
    Fdata->close();
    Fpos->close();
}

template<typename T>
void Ttable2D<T>::fetchTtabledata(std::shared_ptr<RaysAcoustic2D<T>> rays, std::shared_ptr<Data2D<T>> source, const size_t number)
{
    size_t nx = rays->getNx();
    size_t nz = rays->getNz();

    if(number > this->getNtable()-1) rs_error("Ttable2D::fetchTtabledata: Trying to write a table with number that is larger than ntable");
    if(this->getNx() != nx) rs_error("Ttable2D::fetchTtabledata: Nx mismatch.");
    if(this->getNz() != nz) rs_error("Ttable2D::fetchTtabledata: Nz mismatch.");

    if(!this->getAllocated()) {
        this->allocTtable();
    }
    T *data = this->getData();
    T *TT = rays->getTT();

    Point2D<T> *scoords_rays = (source->getGeom())->getScoords();
    Point3D<T> *scoords_tt = (this->getGeom())->getScoords();

    scoords_tt[number].x = scoords_rays[0].x;
    scoords_tt[number].y = 0.0;
    scoords_tt[number].z = scoords_rays[0].y;

    size_t nx_pml = rays->getNx_pml();
    size_t nz_pml = rays->getNz_pml();
    int lpml = rays->getLpml();

    Index Itt(nx,nz);
    Index Irays(nx_pml, nz_pml);
    for(size_t ix=0; ix < nx; ix++){
        for(size_t iz=0; iz < nz; iz++){
            data[Itt(ix,iz)] = TT[Irays(ix+lpml, iz+lpml)];
        }
    }
}

template<typename T>
void Ttable2D<T>::writeTtable(const size_t number)
{
    if(number > this->getNtable()-1) rs_error("Ttable2D::writeTtable: Trying to write a table with number that is larger than ntable");
    if(!this->getAllocated())  rs_error("Ttable2D::writeTtable: Data not allocated. Run fetchTtabledata, before writing the table.");

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable2D::writeTtable: No file assigned. ");
    }

    std::string posfile = this->getFilename()+"-pos";
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(datafile);
    if(status == FILE_ERR) rs_error("Ttable2D::writeTtable: Error opening data file for appending: ", datafile);
    if(Fdata->getN(4) != this->getNtable()) rs_error("Ttable2D::writeTtable: Number of gathers mismatch in datafile.");
    //Get gather size information
    size_t nx = this->getNx_pml();
    size_t ny = 1;
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;

    if(Fdata->getN(1) != nx) rs_error("Ttable2D::writeTtable: Nx mismatch.");
    if(Fdata->getN(2) != ny) rs_error("Ttable2D::writeTtable: Ny mismatch.");
    if(Fdata->getN(3) != nz) rs_error("Ttable2D::writeTtable: Nz mismatch.");

    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    status = Fpos->append(posfile);
    if(status == FILE_ERR) rs_error("Ttable2D::writeTtable: Error opening positions file for appending: ", posfile);
    if(Fpos->getN(2) != this->getNtable()) rs_error("Ttable2D::writeTtable: Number of gathers mismatch in posfile.");
    
    Point3D<T> *scoords = (this->getGeom())->getScoords();
    //Write gather
    Fpos->seekp(Fpos->getStartofdata() + number*3*sizeof(T));
    Fpos->write(&scoords[number].x, 1);
    Fpos->write(&scoords[number].y, 1);
    Fpos->write(&scoords[number].z, 1);
    if(Fpos->getFail()) rs_error("Ttable2D::writeTtable: Error writting table positions to output file");
    Fdata->seekp(Fdata->getStartofdata() + number*n1*sizeof(T));
    Fdata->write(this->getData(), n1);
    if(Fdata->getFail()) rs_error("Ttable2D::writeTtable: Error writting table to output file");
    Fdata->close();
    Fpos->close();
}

template<typename T>
void Ttable2D<T>::interpTtable(std::shared_ptr<Ttable2D<T>> ttablei, T rad) {
    // Variables
    bool status;
    struct kdres *presults;
    double weight, weight_sum;
    int pch;
    double dist;

    // Get sizes of table
    size_t nx = this->getNx_pml();
    size_t ny = 1;
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;

    // Get interpolation position
    double posi[3], tpos[3]; 
    Point3D<T> *scoords = (ttablei->getGeom())->getScoords();
    posi[0] = (double) scoords[0].x;
    posi[1] = (double) scoords[0].y;
    posi[2] = (double) scoords[0].z;

    kdtree *ptree = this->getKdtree();
    // Find neighbors
    presults = kd_nearest_range( ptree, &posi[0], rad );
    int nr=kd_res_size(presults);

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable2D::interpTtable2D: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(datafile);
    if(status == FILE_ERR) rs_error("Ttable2D::interpTtable2D: Error opening data file for input: ", datafile);

    T *data = this->getData();
    T *wrk = this->getWrk();
    weight_sum=0.0;
    if(!this->getAllocated()) rs_error("Ttable2D::interpTtable2D: Data in source ttable is not allocated.");
    if(!ttablei->getAllocated()) rs_error("Ttable2D::interpTtable2D: Data in target ttable is not allocated.");
    // Do the interpolation (Modified Shepard's method)
    for (size_t id=0; id<n1; id++){
    // Initialize data
        data[id]=0.0;
    }
    // TODO: SORT RESULTS BY DISTANCE AND GET THE ONLY CLOSEST 5
    for(int i1=0; i1<nr; i1++){
        /* Fetch index and position of one of the traces in range */
        pch =  kd_res_item( presults, tpos );

        /* Read trace with the corresponding index */
        Fdata->seekp(Fdata->getStartofdata() + pch*n1*sizeof(T));
        Fdata->read(wrk, n1);
        dist=sqrt(SQ(tpos[0]-posi[0]) + SQ(tpos[1]-posi[1])+ SQ(tpos[2]-posi[2]));
        if(dist==0.0){
            for (size_t id=0; id<n1; id++){
                data[id]=wrk[id];
            }
            weight_sum=0.0;
            break;
        }else{
            weight= SQ((rad-dist)/(rad*dist));
            weight_sum+=weight;
            for (size_t id=0; id<n1; id++){
                data[id]+=weight*wrk[id];
            }
        }
        /* go to the next trace */
        kd_res_next( presults );
    }
    if(weight_sum>0){
        for (size_t id=0; id<n1; id++){
            data[id]/=weight_sum;
        }
    }

    // Free results of tree search
    kd_res_free( presults );


    // Get regular coordinates of source and target table and do trilinear interpolation  
    size_t nx_i = ttablei->getNx_pml();
    size_t nz_i = ttablei->getNz_pml();

    T dx = this->getDx();
    T dz = this->getDz();

    T ox = this->getOx();
    T oz = this->getOz();

    T dx_i = ttablei->getDx();
    T dz_i = ttablei->getDz();

    T ox_i = ttablei->getOx();
    T oz_i = ttablei->getOz();
    int ix0,ix1,iz0,iz1;
    T *datai = ttablei->getData();
    T d;

    for (int ix=0; ix < nx_i; ix++){
        ix0 = this->tabfloor(((ix*dx_i + ox_i)-ox)/dx);
        ix1 = ix0 + 1;
        if(ix0<0){
            ix0=0; 
            ix1=0;
        }
        if(ix0 >= nx-1){
            ix0=nx-1;
            ix1=nx-1;
        }
        for (int iz=0; iz < nz_i; iz++){
            iz0 = this->tabfloor(((iz*dz_i + oz_i)-oz)/dz);
            iz1 = iz0 + 1;
            if(iz0<0){
                iz0=0; 
                iz1=0;
            }
            if(iz0 >= nz-1){
                iz0=nz-1;
                iz1=nz-1;
            }
            interpoint2D[0] = data[k2D(ix0,iz0)];
            interpoint2D[1] = data[k2D(ix1,iz0)];
            interpoint2D[2] = data[k2D(ix0,iz1)];
            interpoint2D[3] = data[k2D(ix1,iz1)];

            d = ((iz*dz_i + oz_i) - (iz0*dz + oz))/dz;
            interpoint1D[0] = (1.0-d)*interpoint2D[0] + d*interpoint2D[2];
            interpoint1D[1] = (1.0-d)*interpoint2D[1] + d*interpoint2D[3];

            d = ((ix*dx_i + ox_i) - (ix0*dx + ox))/dx;
            datai[k2Di(ix,iz)] = (1.0-d)*interpoint1D[0] + d*interpoint1D[1];
        }
    }
}

template<typename T>
void Ttable2D<T>::freeTtable()
{
    if(this->getAllocated()) {
        free(this->data);
        free(this->wrk);
    }
}

template<typename T>
Ttable2D<T>::~Ttable2D() {
    // Nothing to do
    if(this->getAllocated()) {
        free(this->data);
        free(this->wrk);
    }
}

// =============== 3D TTABLE CLASS =============== //
template<typename T>
Ttable3D<T>::Ttable3D(std::shared_ptr<ModelEikonal3D<T>> model, int _ntable):Ttable<T>(3, _ntable)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _lpml;
    int _dim;

/* Get necessary parameters from model class */
    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();
    // TODO implement the possibility to pad the traveltime table
    //_lpml = model->getLpml();
    _lpml = 0;
    _dim = model->getDim();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox-(_lpml*_dx));
    this->setOy(_oy-(_lpml*_dy));
    this->setOz(_oz-(_lpml*_dz));
    this->setDim(_dim);
    this->setNtable(_ntable);
    this->setLpml(_lpml);
    this->setData(NULL);

}

template<typename T>
Ttable3D<T>::Ttable3D(std::string tablefile)
{
    bool status;
    /* Set filename */
    this->setFilename(tablefile);
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(tablefile);
    if(status == FILE_ERR) rs_error("Ttable3D::Ttable3D: Error opening data file for reading: ", tablefile);

    std::string posfile = tablefile+"-pos";
    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    status = Fpos->input(posfile);
    if(status == FILE_ERR) rs_error("Ttable3D::Ttable3D: Error opening position file for reading: ", tablefile);

    size_t ntable = Fdata->getN(4);
    if(ntable < 1) rs_error("Number of tables in file ", tablefile,  " less than 1.");
    if(Fpos->getN(2) != ntable) rs_error("Ttable3D::Ttable3D: Number of gathers mismatch in datafile and positions file.");
    this->makeGeom(ntable);
    this->setNtable(ntable);

    /* Get necessary parameters from file */
    this->setNx(Fdata->getN(1));
    this->setNy(Fdata->getN(2));
    this->setNz(Fdata->getN(3));
    this->setDx(Fdata->getD(1));
    this->setDy(Fdata->getD(2));
    this->setDz(Fdata->getD(3));
    this->setOx(Fdata->getO(1));
    this->setOy(Fdata->getO(2));
    this->setOz(Fdata->getO(3));
    this->setLpml(0);
    Fdata->close();

    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;
    if(n1 == 0) rs_error("Ttable::Ttable3D: Zero dimension array. Check Nx, Ny, Nz values.");
    this->setDim(2);
    this->setData(NULL);

    // Read coordinates and populate tree
    kdtree *ptree = this->getKdtree();
    double pos[3];
    T val;
    for (size_t i=0; i < ntable; i++){
        Fpos->read(&val, 1);
        pos[0] = (double) val;
        Fpos->read(&val, 1);
        pos[1] = (double) val;
        Fpos->read(&val, 1);
        pos[2] = (double) val;
        if(kd_insert(ptree, &pos[0], i) != 0){
            rs_error("Ttable::Ttable3D: Error populating tree.");
        }

    }
    Fpos->close();
}

template<typename T>
void Ttable3D<T>::insertSource(std::shared_ptr<rockseis::Data3D<T>> source, bool maptype, int traceno){
    Point3D<T> *pos;
    Point3D<T> *posi;
    size_t ntrace = source->getNtrace();

    if(traceno < 0 || traceno > ntrace-1){
        rs_error("RaysAcoustic3D<T>::insertSource: traceno out of bounds.");
    }

    // Get correct map (source or receiver mapping)
    if(maptype == SMAP) {
        pos = (source->getGeom())->getScoords();
    }else{
        pos = (source->getGeom())->getGcoords();
    }

    posi = (this->getGeom())->getScoords();
    posi[0].x = pos[traceno].x;
    posi[0].y = pos[traceno].y;
    posi[0].z = pos[traceno].z;
}

template<typename T>
void Ttable3D<T>::allocTtable() 
{
    if(this->getAllocated() == false) {
        this->data = (T *) calloc(this->getNx_pml()*this->getNy_pml()*this->getNz_pml(), sizeof(T));
        this->wrk = (T *) calloc(this->getNx_pml()*this->getNy_pml()*this->getNz_pml(), sizeof(T));
        this->setAllocated(true);
    }else{
        free(this->data); free(this->wrk);
        this->data = (T *) calloc(this->getNx_pml()*this->getNy_pml()*this->getNz_pml(), sizeof(T));
        this->wrk = (T *) calloc(this->getNx_pml()*this->getNy_pml()*this->getNz_pml(), sizeof(T));
        this->setAllocated(true);
    }
}

template<typename T>
void Ttable3D<T>::createEmptyttable() 
{
    int ntable = this->getNtable();
    if(ntable < 1) rs_error("Ttable3D::createEmptyTtable3D: number of tables is less than 1.");

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable3D::createEmptyTtable3D: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    Fdata->output(datafile);

    std::string posfile = datafile+"-pos";

    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    Fpos->output(posfile);

    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;
    if(n1 == 0) rs_error("Ttable::createEmptydata: Zero dimension array. Check Nx, Ny, Nz values.");
    Fdata->setN(1, nx);
    Fdata->setN(2, ny);
    Fdata->setN(3, nz);
    Fdata->setN(4,ntable);
    Fdata->setD(1,this->getDx());
    Fdata->setD(2,this->getDy());
    Fdata->setD(3,this->getDz());
    Fdata->setD(4,1);
    Fdata->setO(1,this->getOx());
    Fdata->setO(2,this->getOy());
    Fdata->setO(3,this->getOz());
    Fdata->setO(4,0);
    Fdata->setData_format(sizeof(T));
    Fdata->setType(rockseis::REGULAR);
    Fdata->writeHeader();
    Fdata->seekp(Fdata->getStartofdata());

    Fpos->setN(1,3);
    Fpos->setD(1,1.0);
    Fpos->setN(2, ntable);
    Fpos->setD(2,1.0);
    Fpos->setData_format(sizeof(T));
    Fpos->setType(rockseis::REGULAR);
    Fpos->writeHeader();
    Fpos->seekp(Fpos->getStartofdata());

    this->allocTtable();
    T val=0.0;
    for (int i=0; i < ntable; i++){
        Fpos->write(&val, 1);
        Fpos->write(&val, 1);
        Fpos->write(&val, 1);
        Fdata->write(this->getData(),n1);
    }
    Fdata->close();
    Fpos->close();
}

template<typename T>
void Ttable3D<T>::fetchTtabledata(std::shared_ptr<RaysAcoustic3D<T>> rays, std::shared_ptr<Data3D<T>> source, const size_t number)
{
    size_t nx = rays->getNx();
    size_t ny = rays->getNy();
    size_t nz = rays->getNz();

    if(number > this->getNtable()-1) rs_error("Ttable3D::fetchTtabledata: Trying to write a table with number that is larger than ntable");
    if(this->getNx() != nx) rs_error("Ttable3D::fetchTtabledata: Nx mismatch.");
    if(this->getNy() != ny) rs_error("Ttable3D::fetchTtabledata: Ny mismatch.");
    if(this->getNz() != nz) rs_error("Ttable3D::fetchTtabledata: Nz mismatch.");

    if(!this->getAllocated()) {
        this->allocTtable();
    }
    T *data = this->getData();
    T *TT = rays->getTT();

    Point3D<T> *scoords_rays = (source->getGeom())->getScoords();
    Point3D<T> *scoords_tt = (this->getGeom())->getScoords();

    scoords_tt[number].x = scoords_rays[0].x;
    scoords_tt[number].y = scoords_rays[0].y;
    scoords_tt[number].z = scoords_rays[0].z;

    size_t nx_pml = rays->getNx_pml();
    size_t ny_pml = rays->getNy_pml();
    size_t nz_pml = rays->getNz_pml();
    int lpml = rays->getLpml();

    Index Itt(nx,ny,nz);
    Index Irays(nx_pml, ny_pml, nz_pml);
    for(size_t ix=0; ix < nx; ix++){
        for(size_t iy=0; iy < ny; iy++){
            for(size_t iz=0; iz < nz; iz++){
                data[Itt(ix,iy,iz)] = TT[Irays(ix+lpml,iy+lpml,iz+lpml)];
            }
        }
    }

}

template<typename T>
void Ttable3D<T>::writeTtable(const size_t number)
{
    if(number > this->getNtable()-1) rs_error("Ttable3D::writeTtable: Trying to write a table with number that is larger than ntable");
    if(!this->getAllocated())  rs_error("Ttable3D::writeTtable: Data not allocated. Run fetchTtabledata, before writing the table.");

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable3D::writeTtable: No file assigned. ");
    }

    std::string posfile = this->getFilename()+"-pos";
    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(datafile);
    if(status == FILE_ERR) rs_error("Ttable3D::writeTtable: Error opening data file for appending: ", datafile);
    if(Fdata->getN(4) != this->getNtable()) rs_error("Ttable3D::writeTtable: Number of gathers mismatch in datafile.");
    //Get gather size information
    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;

    if(Fdata->getN(1) != nx) rs_error("Ttable3D::writeTtable: Nx mismatch.");
    if(Fdata->getN(2) != ny) rs_error("Ttable3D::writeTtable: Ny mismatch.");
    if(Fdata->getN(3) != nz) rs_error("Ttable3D::writeTtable: Nz mismatch.");

    std::shared_ptr<rockseis::File> Fpos (new rockseis::File());
    status = Fpos->append(posfile);
    if(status == FILE_ERR) rs_error("Ttable3D::writeTtable: Error opening positions file for appending: ", posfile);
    if(Fpos->getN(2) != this->getNtable()) rs_error("Ttable3D::writeTtable: Number of gathers mismatch in posfile.");
    
    Point3D<T> *scoords = (this->getGeom())->getScoords();
    //Write gather
    Fpos->seekp(Fpos->getStartofdata() + number*3*sizeof(T));
    Fpos->write(&scoords[number].x, 1);
    Fpos->write(&scoords[number].y, 1);
    Fpos->write(&scoords[number].z, 1);
    if(Fpos->getFail()) rs_error("Ttable3D::writeTtable: Error writting table positions to output file");
    Fdata->seekp(Fdata->getStartofdata() + number*n1*sizeof(T));
    Fdata->write(this->getData(), n1);
    if(Fdata->getFail()) rs_error("Ttable3D::writeTtable: Error writting table to output file");
    Fdata->close();
    Fpos->close();
}

template<typename T>
void Ttable3D<T>::interpTtable(std::shared_ptr<Ttable3D<T>> ttablei, T rad) {
    // Variables
    bool status;
    struct kdres *presults;
    double weight, weight_sum;
    int pch;
    double dist;

    // Get sizes of table
    size_t nx = this->getNx_pml();
    size_t ny = this->getNy_pml();
    size_t nz = this->getNz_pml();
    size_t n1 = nx*ny*nz;

    // Get interpolation position
    double posi[3], tpos[3]; 
    Point3D<T> *scoords = (ttablei->getGeom())->getScoords();
    posi[0] = (double) scoords[0].x;
    posi[1] = (double) scoords[0].y;
    posi[2] = (double) scoords[0].z;

    // Find neighbors
    kdtree *ptree = this->getKdtree();
    presults = kd_nearest_range( ptree, &posi[0], rad );
    int nr=kd_res_size(presults);

    std::string datafile = this->getFilename();
    if(datafile.empty()){
        rs_error("Ttable3D::interpTtable: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->input(datafile);
    if(status == FILE_ERR) rs_error("Ttable3D::interpTtable: Error opening data file for input: ", datafile);

    T *data = this->getData();
    T *wrk = this->getWrk();
    weight_sum=0.0;
    if(!this->getAllocated()) rs_error("Ttable3D::interpTtable: Data in source ttable is not allocated.");
    if(!ttablei->getAllocated()) rs_error("Ttable3D::interpTtable: Data in target ttable is not allocated.");
    // Do the interpolation (Modified Shepard's method)
    for (size_t id=0; id<n1; id++){
    // Initialize data
        data[id]=0.0;
    }
    // TODO: SORT RESULTS BY DISTANCE AND GET THE ONLY CLOSEST 5
    for(int i1=0; i1<nr; i1++){
        /* Fetch index and position of one of the traces in range */
        pch =  kd_res_item( presults, tpos );

        /* Read trace with the corresponding index */
        Fdata->seekp(Fdata->getStartofdata() + pch*n1*sizeof(T));
        Fdata->read(wrk, n1);
        dist=sqrt(SQ(tpos[0]-posi[0]) + SQ(tpos[1]-posi[1])+ SQ(tpos[2]-posi[2]));
        if(dist==0.0){
            for (size_t id=0; id<n1; id++){
                data[id]=wrk[id];
            }
            weight_sum=0.0;
            break;
        }else{
            weight= SQ((rad-dist)/(rad*dist));
            weight_sum+=weight;
            for (size_t id=0; id<n1; id++){
                data[id]+=weight*wrk[id];
            }
        }
        /* go to the next trace */
        kd_res_next( presults );
    }
    if(weight_sum>0){
        for (size_t id=0; id<n1; id++){
            data[id]/=weight_sum;
        }
    }

    // Free results of tree search
    kd_res_free( presults );


    // Get regular coordinates of source and target table and do trilinear interpolation  
    size_t nx_i = ttablei->getNx();
    size_t ny_i = ttablei->getNy();
    size_t nz_i = ttablei->getNz();

    T dx = this->getDx();
    T dy = this->getDy();
    T dz = this->getDz();

    T ox = this->getOx();
    T oy = this->getOy();
    T oz = this->getOz();

    T dx_i = ttablei->getDx();
    T dy_i = ttablei->getDy();
    T dz_i = ttablei->getDz();

    T ox_i = ttablei->getOx();
    T oy_i = ttablei->getOy();
    T oz_i = ttablei->getOz();
    int ix0,ix1,iy0,iy1,iz0,iz1;
    T *datai = ttablei->getData();
    T d;

    for (int ix=0; ix < nx_i; ix++){
        ix0 = this->tabfloor(((ix*dx_i + ox_i)-ox)/dx);
        ix1 = ix0 + 1;
        if(ix0<0){
            ix0=0; 
            ix1=0;
        }
        if(ix0 >= nx-1){
            ix0=nx-1;
            ix1=nx-1;
        }
        for (int iy=0; iy < ny_i; iy++){
            iy0 = this->tabfloor(((iy*dy_i + oy_i)-oy)/dy);
            iy1 = iy0 + 1;
            if(iy0<0){
                iy0=0; 
                iy1=0;
            }
            if(iy0 >= ny-1){
                iy0=ny-1;
                iy1=ny-1;
            }
            for (int iz=0; iz < nz_i; iz++){
                iz0 = this->tabfloor(((iz*dz_i + oz_i)-oz)/dz);
                iz1 = iz0 + 1;
                if(iz0<0){
                    iz0=0; 
                    iz1=0;
                }
                if(iz0 >= nz-1){
                    iz0=nz-1;
                    iz1=nz-1;
                }
                interpoint3D[0] = data[k3D(ix0,iy0,iz0)];
                interpoint3D[1] = data[k3D(ix1,iy0,iz0)];
                interpoint3D[2] = data[k3D(ix0,iy1,iz0)];
                interpoint3D[3] = data[k3D(ix1,iy1,iz0)];
                interpoint3D[4] = data[k3D(ix0,iy0,iz1)];
                interpoint3D[5] = data[k3D(ix1,iy0,iz1)];
                interpoint3D[6] = data[k3D(ix0,iy1,iz1)];
                interpoint3D[7] = data[k3D(ix1,iy1,iz1)];

                d = ((iz*dz_i + oz_i) - (iz0*dz + oz))/dz;
                interpoint2D[0] = (1.0-d)*interpoint3D[0] + d*interpoint3D[4];
                interpoint2D[1] = (1.0-d)*interpoint3D[1] + d*interpoint3D[5];
                interpoint2D[2] = (1.0-d)*interpoint3D[2] + d*interpoint3D[6];
                interpoint2D[3] = (1.0-d)*interpoint3D[3] + d*interpoint3D[7];

                d = ((iy*dy_i + oy_i) - (iy0*dy + oy))/dy;
                interpoint1D[0] = (1.0-d)*interpoint2D[0] + d*interpoint2D[2];
                interpoint1D[1] = (1.0-d)*interpoint2D[1] + d*interpoint2D[3];

                d = ((ix*dx_i + ox_i) - (ix0*dx + ox))/dx;
                datai[k3Di(ix,iy,iz)] = (1.0-d)*interpoint1D[0] + d*interpoint1D[1];
            }
        }
    }
}

template<typename T>
void Ttable3D<T>::freeTtable()
{
    if(this->getAllocated()) {
        free(this->data);
        free(this->wrk);
    }
}

template<typename T>
Ttable3D<T>::~Ttable3D() {
    // Nothing to do
    if(this->getAllocated()) {
        free(this->data);
        free(this->wrk);
    }
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Ttable<float>;
template class Ttable<double>;

template class Ttable2D<float>;
template class Ttable2D<double>;

template class Ttable3D<float>;
template class Ttable3D<double>;

}
