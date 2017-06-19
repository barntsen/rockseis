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


template<typename T>
void Image2D<T>::crossCorr(T *ws, int pads, T* wr, int padr)
{
    if(!allocated) this->allocateImage();
	int ix, iz, ihx, ihz;
	int nhx = this->getNhx();
	int nhz = this->getNhz();
	int nx = this->getNx();
	int nxs = nx+2*pads;
	int nxr = nx+2*padr;
	int nz = this->getNz();
	int hx, hz;
	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihz=0; ihz<nhz; ihz++){
			hz= -(nhz-1)/2 + ihz;
			for (ix=0; ix<nx; ix++){
				if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
				{
					for (iz=0; iz<nz; iz++){
						if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz))
							imagedata[ki2D(ix,iz,ihx,ihz)] += ws[ks2D(ix+pads, iz+pads)]*wr[kr2D(ix+padr, iz+padr)];
					}	
				}
			}
		}
	}
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
void Image3D<T>::freeImage()
{
    if(allocated){
        free(imagedata);
    }
    allocated = false;
}

template<typename T>
void Image3D<T>::crossCorr(T *ws, int pads, T* wr, int padr)
{
    if(!allocated) this->allocateImage();
	int ix, iy, iz, ihx, ihy, ihz;
	int nhx = this->getNhx();
	int nhy = this->getNhy();
	int nhz = this->getNhz();
	int nx = this->getNx();
	int nxs = nx + 2*pads;
	int nxr = nx + 2*padr;
	int ny = this->getNy();
	int nys = ny + 2*pads;
	int nyr = ny + 2*padr;
	int nz = this->getNz();
	int hx, hy, hz;
	for (ihx=0; ihx<nhx; ihx++){
		hx= -(nhx-1)/2 + ihx;
		for (ihy=0; ihy<nhy; ihy++){
			hy= -(nhy-1)/2 + ihy;
			for (ihz=0; ihz<nhz; ihz++){
				hz= -(nhz-1)/2 + ihz;
				for (ix=0; ix<nx; ix++){
					if( ((ix-hx) >= 0) && ((ix-hx) < nx) && ((ix+hx) >= 0) && ((ix+hx) < nx))
						for (iy=0; iy<ny; iy++){
							if( ((iy-hy) >= 0) && ((iy-hy) < ny) && ((iy+hy) >= 0) && ((iy+hy) < ny))
								for (iz=0; iz<nz; iz++){
									if( ((iz-hz) >= 0) && ((iz-hz) < nz) && ((iz+hz) >= 0) && ((iz+hz) < nz))
										imagedata[ki3D(ix,iy,iz,ihx,ihy,ihz)] += ws[ks3D(ix+pads,iy+pads,iz+pads)]*wr[kr3D(ix+padr,iy+padr,iz+padr)];
								}	
						}
				}
			}
		}
	}
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


