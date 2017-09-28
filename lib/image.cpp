#include "image.h"

#define ki2D(i,j,k,l) ((l)*nhx*nz*nx + (k)*nx*nz + (j)*nx +(i))
#define ks2D(i,j) ((j)*nxs + (i))
#define kr2D(i,j) ((j)*nxr + (i))

#define ki3D(i,j,k,l,m,n) ((n)*nhy*nhx*nx*ny*nz + (m)*nhx*nx*ny*nz + (l)*nx*ny*nz + (k)*nx*ny + (j)*nx + (i))
#define ks3D(i,j,k) ((k)*nxs*nys + (j)*nxs + (i))
#define kr3D(i,j,k) ((k)*nxr*nyr + (j)*nxr + (i))

namespace rockseis {


/// =============== 2D IMAGE CLASS =============== //

// constructor
template<typename T>
Image2D<T>::Image2D(std::string _imagefile): Model<T>(2)
{
    bool status;
    long int _nx, _ny, _nz;
    int _nhx, _nhz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

    /* Get necessary parameters from file */

    //Opening file for reading
    imagefile = _imagefile;
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(imagefile.c_str());
    if(status == FILE_ERR){
	    rs_error("Image2D::Error reading from ", imagefile );
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Image2D::Numerical precision in ", imagefile, " mismatch with image class contructor.");
    }

    _nx=Fin->getN(1);
    _ny=Fin->getN(2);
    _nz=Fin->getN(3);
    _dx=Fin->getD(1);
    _dy=Fin->getD(2);
    _dz=Fin->getD(3);
    _ox=Fin->getO(1);
    _oy=Fin->getO(2);
    _oz=Fin->getO(3);

    _nhx=Fin->getN(4);
    _nhz=Fin->getN(6);

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNhx(_nhx);
    this->setNhz(_nhz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    allocated = false;
}

template<typename T>
Image2D<T>::Image2D(std::string _imagefile, std::shared_ptr<ModelAcoustic2D<T>> model, int nhx, int nhz): Model<T>(2)
{
    long int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setNhx(nhx);
    this->setNhz(nhz);
    imagefile = _imagefile;
    allocated = false;
}

template<typename T>
Image2D<T>::Image2D(std::string _imagefile, std::shared_ptr<ModelElastic2D<T>> model, int nhx, int nhz): Model<T>(2)
{
    long int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setNhx(nhx);
    this->setNhz(nhz);
    imagefile = _imagefile;
    allocated = false;
}

template<typename T>
bool Image2D<T>::read()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image2D::read: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(imagefile.c_str());
    if(status == FILE_ERR){
	    rs_error("Image2D::read: Error reading from " , imagefile ,". ");
    }
    long int nx, nz;
    int nhx, nhz;
    nx = this->getNx();
    nz = this->getNz();
    nhx = this->getNhx();
    nhz = this->getNhz();

    // Allocate space for image;
    if(!allocated) this->allocateImage();

    if ( Fin->is_open() )
    {
        //Read image from file
        Fin->read(&imagedata[0], nx*nz*nhx*nhz);
        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Image2D<T>::createEmpty()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image2D::createEmpty: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(imagefile.c_str());

	long int nx, nz;
	int nhx, nhz;
	T dx,dz; 
	T ox,oz;
	nx = this->getNx();
	nz = this->getNz();
	nhx = this->getNhx();
	nhz = this->getNhz();
	dx = this->getDx();
	dz = this->getDz();
	ox = this->getOx();
	oz = this->getOz();

    if ( Fout->is_open())
    {
        // Create geometry
        Fout->setN(1,nx);
        Fout->setD(1,dx);
        Fout->setO(1,ox);
        Fout->setN(3,nz);
        Fout->setD(3,dz);
        Fout->setO(3,oz);
        Fout->setN(4,nhx);
        Fout->setD(4,dx);
        Fout->setO(4,-dx*(nhx-1)/2);
        Fout->setN(6,nhz);
        Fout->setD(6,dz);
        Fout->setO(6,-dz*(nhz-1)/2);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(REGULAR);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        T val = 0.0;
        // Write data
        for (long int i=0; i < nx*nz*nhx*nhz; i++)
        {
            Fout->write(&val, 1);
        }
        if(Fout->getFail()) rs_error("Image2D::createEmpty: Failed to write data to file");
        status = FILE_OK;
    }else{
        rs_error("Image2D::createEmpty: Error opening file. ");
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Image2D<T>::stackImage(std::string infile)
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image2D::stackImage: No file assigned. ");
    }

    if(infile.empty()){
	    rs_error("Image2D::stackImage: Input filename is empty. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->append(imagefile.c_str());
    if ( !Fout->is_open()) rs_error("Image2D::stackImage: Failed to open output image.");

    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    Fin->input(infile.c_str());
    if( !Fin->is_open() ) rs_error("Image2D::stackImage: Failed to open input image.");

    long int nxg, nzg;
    int nhxg, nhzg;
	T dxg, dzg, oxg, ozg;
    nxg = this->getNx();
    nzg = this->getNz();
    nhxg = this->getNhx();
    nhzg = this->getNhz();
    dxg = this->getDx();
    dzg = this->getDz();
    oxg = this->getOx();
    ozg = this->getOz();

    long int nxl, nzl;
    int nhxl, nhzl;
	T dxl, dzl, oxl, ozl;
    nxl = Fin->getN(1);
    nzl = Fin->getN(3);
    nhxl = Fin->getN(4);
    nhzl = Fin->getN(6);
	dxl = Fin->getD(1);
	dzl = Fin->getD(3);
	oxl = Fin->getO(1);
	ozl = Fin->getO(3);

	if(nhxg != nhxl || nhzg != nhzl || dxg != dxl || dzg != dzl){
		rs_error("Image2D::stackImage: Images are not compatible. Cannot stack.");
	}

    T *trcin;
    trcin = (T *) calloc(nxl, sizeof(T));
    if(trcin == NULL) rs_error("Image2D::stackImage: Failed to allocate memory for input trace.");

    T *trcout;
    trcout = (T *) calloc(nxg, sizeof(T));
    if(trcout == NULL) rs_error("Image2D::stackImage: Failed to allocate memory for input trace.");

	// Stack data
	int ix, iz, ihx, ihz;
    int ix_start, iz_start;
	ix_start = (int) (oxl/dxl - oxg/dxg);
	iz_start = (int) (ozl/dzl - ozg/dzg);
	Index Iin(nxl, nzl, nhxl, nhzl);
	Index Iout(nxg, nzg, nhxl, nhzl);
	for(ihz=0; ihz<nhzl; ihz++) {
		for(ihx=0; ihx<nhxl; ihx++) {
			for(iz=0; iz<nzl; iz++) {
				if((iz + iz_start) >= 0 && (iz + iz_start) < nzg){
					// Read traces 
					Fin->read(trcin, nxl, Iin(0,iz,ihx,ihz)*sizeof(T));
					Fout->read(trcout, nxg, Iout(0,(iz+iz_start),ihx,ihz)*sizeof(T));
					if(Fin->getFail()) rs_error("Image2D::stackImage: Failed to read data from file");
					if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to write data to file");
					for(ix = 0; ix < nxl; ix++) {
						if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
							trcout[ix + ix_start] += trcin[ix];
						}
					}
					// Write trc 
					Fout->write(trcout, nxg, Iout(0, (iz+iz_start),ihx,ihz)*sizeof(T));
					if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to write data to file");
				}
			}
		}
	}
	Fin->close();
	Fout->close();
	free(trcin);
	free(trcout);
	status = FILE_OK;
    return status;
}

template<typename T>
bool Image2D<T>::write()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image2D::writeImage: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    if(strcmp(imagefile.c_str(), "stdout")) {
        Fout->output(imagefile.c_str());
    }else{
        Fout->output();
    }

    long int nx, nz;
    int nhx, nhz;
    T dx,dz; 
    T ox,oz;
    nx = this->getNx();
    nz = this->getNz();
    nhx = this->getNhx();
    nhz = this->getNhz();
    dx = this->getDx();
    dz = this->getDz();
    ox = this->getOx();
    oz = this->getOz();

    if ( Fout->is_open() && allocated)
    {
        // Create geometry
        Fout->setN(1,nx);
        Fout->setD(1,dx);
        Fout->setO(1,ox);
        Fout->setN(3,nz);
        Fout->setD(3,dz);
        Fout->setO(3,oz);
        Fout->setN(4,nhx);
        Fout->setD(4,dx);
        Fout->setO(4,-dx*(nhx-1)/2);
        Fout->setN(6,nhz);
        Fout->setD(6,dz);
        Fout->setO(6,-dz*(nhz-1)/2);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(REGULAR);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        // Write image
        Fout->write(&imagedata[0], nx*nz*nhx*nhz);

        status = FILE_OK;
    }else{
	    rs_error("Image2D::writeImage: Image is not allocated. ");
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
void Image2D<T>::allocateImage()
{
    // Allocate the memory for the image
    imagedata = (T *) calloc(this->getNhx()*this->getNhz()*this->getNx()*this->getNz(), sizeof(T));
    if(imagedata == NULL) rs_error("Image2D::allocate: Failed to allocate memory for imagedata");
    allocated = true;
}

template<typename T>
void Image2D<T>::freeImage()
{
    if(allocated){
        free(imagedata);
    }
    allocated = false;
}




// destructor
template<typename T>
Image2D<T>::~Image2D() {
    // Freeing all variables
    if(allocated){
        free(imagedata);
    }
}

/// =============== 3D IMAGE CLASS =============== //

// constructor
template<typename T>
Image3D<T>::Image3D(std::string imagefile): Model<T>(3)
{
    bool status;
    long int _nx, _ny, _nz;
    int _nhx, _nhy, _nhz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

/* Get necessary parameters from file */

    //Opening file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(imagefile.c_str());
    if(status == FILE_ERR){
	    rs_error("Image3D::Error reading from ", imagefile );
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Image3D::Numerical precision in ", imagefile, " mismatch with image class contructor.");
    }

    _nx=Fin->getN(1);
    _ny=Fin->getN(2);
    _nz=Fin->getN(3);
    _dx=Fin->getD(1);
    _dy=Fin->getD(2);
    _dz=Fin->getD(3);
    _ox=Fin->getO(1);
    _oy=Fin->getO(2);
    _oz=Fin->getO(3);

    _nhx=Fin->getN(4);
    _nhy=Fin->getN(5);
    _nhz=Fin->getN(6);

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNhx(_nhx);
    this->setNhy(_nhy);
    this->setNhz(_nhz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    allocated = false;
}

template<typename T>
Image3D<T>::Image3D(std::string _imagefile, std::shared_ptr<ModelAcoustic3D<T>> model, int nhx, int nhy, int nhz): Model<T>(3)
{
    long int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setNhx(nhx);
    this->setNhy(nhy);
    this->setNhz(nhz);
    imagefile = _imagefile;
    allocated = false;
}

template<typename T>
Image3D<T>::Image3D(std::string _imagefile, std::shared_ptr<ModelElastic3D<T>> model, int nhx, int nhy, int nhz): Model<T>(3)
{
    long int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 

    _nx=model->getNx();
    _ny=model->getNy();
    _nz=model->getNz();
    _dx=model->getDx();
    _dy=model->getDy();
    _dz=model->getDz();
    _ox=model->getOx();
    _oy=model->getOy();
    _oz=model->getOz();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setNhx(nhx);
    this->setNhy(nhy);
    this->setNhz(nhz);
    imagefile = _imagefile;
    allocated = false;
}

template<typename T>
bool Image3D<T>::read()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image3D::read: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(imagefile.c_str());
    if(status == FILE_ERR){
	    rs_error("Image3D::read: Error reading from " , imagefile ,". ");
    }
    long int nx, ny, nz;
    int nhx, nhy, nhz;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    nhx = this->getNhx();
    nhy = this->getNhy();
    nhz = this->getNhz();

    // Allocate space for image;
    if(!allocated) this->allocateImage();

    if ( Fin->is_open() )
    {
        //Read image from file
        Fin->read(&imagedata[0], nx*ny*nz*nhx*nhy*nhz);
        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Image3D<T>::write()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image3D::writeImage: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    if(strcmp(imagefile.c_str(), "stdout")) {
        Fout->output(imagefile.c_str());
    }else{
        Fout->output();
    }

    long int nx, ny, nz;
    int nhx, nhy, nhz;
    T dx,dy,dz; 
    T ox,oy,oz;
    nx = this->getNx();
    ny = this->getNy();
    nz = this->getNz();
    nhx = this->getNhx();
    nhy = this->getNhy();
    nhz = this->getNhz();
    dx = this->getDx();
    dy = this->getDy();
    dz = this->getDz();
    ox = this->getOx();
    oy = this->getOy();
    oz = this->getOz();

    if ( Fout->is_open() && allocated)
    {
        // Create geometry
        Fout->setN(1,nx);
        Fout->setD(1,dx);
        Fout->setO(1,ox);
        Fout->setN(2,ny);
        Fout->setD(2,dy);
        Fout->setO(2,oy);
        Fout->setN(3,nz);
        Fout->setD(3,dz);
        Fout->setO(3,oz);
        Fout->setN(4,nhx);
        Fout->setD(4,dx);
        Fout->setO(4,-dx*(nhx-1)/2);
        Fout->setN(5,nhy);
        Fout->setD(5,dy);
        Fout->setO(5,-dy*(nhy-1)/2);
        Fout->setN(6,nhz);
        Fout->setD(6,dz);
        Fout->setO(6,-dz*(nhz-1)/2);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(REGULAR);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        // Write image
        Fout->write(&imagedata[0], nx*ny*nz*nhx*nhy*nhz);

        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
void Image3D<T>::allocateImage()
{
    // Allocate the memory for the image
    imagedata = (T *) calloc(this->getNhx()*this->getNhy()*this->getNhz()*this->getNx()*this->getNy()*this->getNz(), sizeof(T));
    if(imagedata == NULL) rs_error("Image3D::allocate: Failed to allocate memory for imagedata");
    allocated = true;
}

template<typename T>
bool Image3D<T>::createEmpty()
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image3D::createEmpty: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(imagefile.c_str());

	long int nx, ny, nz;
	int nhx, nhy, nhz;
	T dx,dy,dz; 
	T ox,oy,oz;
	nx = this->getNx();
	ny = this->getNy();
	nz = this->getNz();
	nhx = this->getNhx();
	nhy = this->getNhy();
	nhz = this->getNhz();
	dx = this->getDx();
	dy = this->getDy();
	dz = this->getDz();
	ox = this->getOx();
	oy = this->getOy();
	oz = this->getOz();

    if ( Fout->is_open())
    {
        // Create geometry
        Fout->setN(1,nx);
        Fout->setD(1,dx);
        Fout->setO(1,ox);
        Fout->setN(2,ny);
        Fout->setD(2,dy);
        Fout->setO(2,oy);
        Fout->setN(3,nz);
        Fout->setD(3,dz);
        Fout->setO(3,oz);
        Fout->setN(4,nhx);
        Fout->setD(4,dx);
        Fout->setO(4,-dx*(nhx-1)/2);
        Fout->setN(5,nhy);
        Fout->setD(5,dx);
        Fout->setO(5,-dy*(nhy-1)/2);
        Fout->setN(6,nhz);
        Fout->setD(6,dz);
        Fout->setO(6,-dz*(nhz-1)/2);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(REGULAR);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        T val = 0.0;
        // Write data
        for (long int i=0; i < nx*ny*nz*nhx*nhy*nhz; i++)
        {
            Fout->write(&val, 1);
        }
        if(Fout->getFail()) rs_error("Image3D::createEmpty: Failed to write data to file");
        status = FILE_OK;
    }else{
        rs_error("Image3D::createEmpty: Error opening file. ");
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Image3D<T>::stackImage(std::string infile)
{
    bool status;
    if(imagefile.empty()){
	    rs_error("Image3D::stackImage: No file assigned. ");
    }

    if(infile.empty()){
	    rs_error("Image3D::stackImage: Input filename is empty. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->append(imagefile.c_str());
    if ( !Fout->is_open()) rs_error("Image3D::stackImage: Failed to open output image.");

    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    Fin->input(infile.c_str());
    if( !Fin->is_open() ) rs_error("Image3D::stackImage: Failed to open input image.");

    long int nxg, nyg,nzg;
    int nhxg, nhyg, nhzg;
	T dxg, dyg, dzg, oxg, oyg, ozg;
    nxg = this->getNx();
    nyg = this->getNy();
    nzg = this->getNz();
    nhxg = this->getNhx();
    nhyg = this->getNhy();
    nhzg = this->getNhz();
    dxg = this->getDx();
    dyg = this->getDy();
    dzg = this->getDz();
    oxg = this->getOx();
    oyg = this->getOy();
    ozg = this->getOz();

    long int nxl, nyl, nzl;
    int nhxl, nhyl, nhzl;
	T dxl, dyl, dzl, oxl, oyl, ozl;
    nxl = Fin->getN(1);
    nyl = Fin->getN(2);
    nzl = Fin->getN(3);
    nhxl = Fin->getN(4);
    nhyl = Fin->getN(5);
    nhzl = Fin->getN(6);
	dxl = Fin->getD(1);
	dyl = Fin->getD(2);
	dzl = Fin->getD(3);
	oxl = Fin->getO(1);
	oyl = Fin->getO(2);
	ozl = Fin->getO(3);

	if(nhxg != nhxl || nhyg != nhyl || nhzg != nhzl || dxg != dxl || dyg != dyl || dzg != dzl){
		rs_error("Image3D::stackImage: Images are not compatible. Cannot stack.");
	}

    T *trcin;
    trcin = (T *) calloc(nxl, sizeof(T));
    if(trcin == NULL) rs_error("Image3D::stackImage: Failed to allocate memory for input trace.");

    T *trcout;
    trcout = (T *) calloc(nxg, sizeof(T));
    if(trcout == NULL) rs_error("Image3D::stackImage: Failed to allocate memory for input trace.");

	// Stack data
	int ix, iy, iz, ihx, ihy, ihz;
    int ix_start, iy_start, iz_start;
	ix_start = (int) (oxl/dxl - oxg/dxg);
	iy_start = (int) (oyl/dyl - oyg/dyg);
	iz_start = (int) (ozl/dzl - ozg/dzg);
	Index Iin(nxl, nyl, nzl, nhxl, nhyl, nhzl);
	Index Iout(nxg, nyg, nzg, nhxl, nhyl, nhzl);
    for(ihz=0; ihz<nhzl; ihz++) {
        for(ihy=0; ihy<nhyl; ihy++) {
            for(ihx=0; ihx<nhxl; ihx++) {
                for(iz=0; iz<nzl; iz++) {
                    if((iz + iz_start) >= 0 && (iz + iz_start) < nzg){
                        for(iy=0; iy<nyl; iy++) {
                            if((iy + iy_start) >= 0 && (iy + iy_start) < nyg){
                                // Read traces 
                                Fin->read(trcin, nxl, Iin(0,iy,iz,ihx,ihy,ihz)*sizeof(T));
                                Fout->read(trcout, nxg, Iout(0, (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)*sizeof(T));
                                if(Fin->getFail()) rs_error("Image3D::stackImage: Failed to read data from file");
                                if(Fout->getFail()) rs_error("Image3D::stackImage: Failed to write data to file");
                                for(ix = 0; ix < nxl; ix++) {
                                    if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
                                        trcout[ix + ix_start] += trcin[ix];
                                    }
                                }
                                // Write trc 
                                Fout->write(trcout, nxg, Iout(0, (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)*sizeof(T));
                                if(Fout->getFail()) rs_error("Image3D::stackImage: Failed to write data to file");
                            }
                        }
                    }
                }
            }
        }
    }
	Fin->close();
	Fout->close();
	free(trcin);
	free(trcout);
	status = FILE_OK;
    return status;
}


template<typename T>
void Image3D<T>::freeImage()
{
    if(allocated){
        free(imagedata);
    }
    allocated = false;
}

// destructor
template<typename T>
Image3D<T>::~Image3D() {
    // Freeing all variables
    if(allocated){
        free(imagedata);
    }
    allocated = false;
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Image2D<float>;
template class Image2D<double>;

template class Image3D<float>;
template class Image3D<double>;

}


