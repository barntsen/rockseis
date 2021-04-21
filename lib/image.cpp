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

   if(_nhx == 0) _nhx = 1;
   if(_nhz == 0) _nhz = 1;

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
   incore = false;
}

template<typename T>
Image2D<T>::Image2D(std::string _imagefile, std::shared_ptr<ModelEikonal2D<T>> model, int nhx, int nhz): Model<T>(2)
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
   incore = false;
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
   incore = false;
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
   incore = false;
}

template<typename T>
Image2D<T>::Image2D(const int _nx, const int _nz,  const int _nhx, const int _nhz, const T _dx, const T _dz, const T _ox, const T _oz)
{
   this->setNx(_nx);
   this->setNy(1);
   this->setNz(_nz);
   this->setDx(_dx);
   this->setDy(1.0);
   this->setDz(_dz);
   this->setOx(_ox);
   this->setOy(0.0);
   this->setOz(_oz);
   this->setNhx(_nhx);
   this->setNhz(_nhz);
   allocated = false;
   incore = false;
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
   trcin = (T *) calloc(nxl*nzl, sizeof(T));
   if(trcin == NULL) rs_error("Image2D::stackImage: Failed to allocate memory for input trace.");

   T *trcout;
   trcout = (T *) calloc(nxg*nzg, sizeof(T));
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
         // Read traces 
         Fin->read(trcin, nxl*nzl, Iin(0,0,ihx,ihz)*sizeof(T));
         Fout->read(trcout, nxg*nzg, Iout(0,0,ihx,ihz)*sizeof(T));
         if(Fin->getFail()) rs_error("Image2D::stackImage: Failed to read data from file");
         if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to write data to file");
         // Stack
         for(iz=0; iz<nzl; iz++) {
            if((iz + iz_start) >= 0 && (iz + iz_start) < nzg){
               for(ix = 0; ix < nxl; ix++) {
                  if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
                     trcout[(iz + iz_start)*nxg + (ix + ix_start)] += trcin[iz*nxl + ix];
                  }
               }
            }
         }
         // Write traces
         Fout->write(trcout, nxg*nzg, Iout(0, 0, ihx, ihz)*sizeof(T));
         if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to write data to file");
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
bool Image2D<T>::stackImage_parallel(std::string infile, int padlx, int padhx, int padlz, int padhz)
{
   bool status;
   if(infile.empty()){
      rs_error("Image2D::stackImage: Input filename is empty. ");
   }

   std::shared_ptr<rockseis::File> Fout (new rockseis::File());
   status = Fout->append(infile);
   if(status == FILE_ERR) rs_error("Image2D::stackImage: Failed to open output image.");

   long int nxg, nzg;
   int nhxg, nhzg;
   T dxg, dzg, oxg, ozg;
   nxg = Fout->getN(1);
   nzg = Fout->getN(3);
   nhxg = Fout->getN(4);
   nhzg = Fout->getN(6);
   dxg = Fout->getD(1);
   dzg = Fout->getD(3);
   oxg = Fout->getO(1);
   ozg = Fout->getO(3);

   long int nxl, nzl;
   int nhxl, nhzl;
   T dxl, dzl, oxl, ozl;
   nxl = this->getNx();
   nzl = this->getNz();
   nhxl = this->getNhx();
   nhzl = this->getNhz();
   dxl = this->getDx();
   dzl = this->getDz();
   oxl = this->getOx();
   ozl = this->getOz();


   if(nhxg != nhxl || nhzg != nhzl || dxg != dxl || dzg != dzl){
      rs_error("Image2D::stackImage_parallel: Images are not compatible. Cannot stack.");
   }

   T *trcin = this->getImagedata();
   T trcout;

   // Stack data
   int ix, iz, ihx, ihz;
   int ix_start, iz_start;
   ix_start = (int) (oxl/dxl - oxg/dxg);
   iz_start = (int) (ozl/dzl - ozg/dzg);

   Index Iin(nxl, nzl, nhxl, nhzl);
   Index Iout(nxg, nzg, nhxl, nhzl);
   for(ihz=0; ihz<nhzl; ihz++) {
      for(ihx=0; ihx<nhxl; ihx++) {
         // Stack
         for(iz = padlz; iz < nzl-padhz; iz++) {
            if((iz + iz_start) >= 0 && (iz + iz_start) < nzg){
               for(ix = padlx; ix < nxl-padhx; ix++) {
                  if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
                     // Read traces 
                     Fout->read(&trcout, 1, Iout(ix + ix_start,iz + iz_start,ihx,ihz)*sizeof(T));
                     if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to read data from file");
                     trcout += trcin[Iin(ix,iz,ihx,ihz)];
                     // Write traces
                     Fout->write(&trcout, 1, Iout(ix + ix_start,iz + iz_start,ihx,ihz)*sizeof(T));
                     if(Fout->getFail()) rs_error("Image2D::stackImage: Failed to write data to file");
                  }
               }
            }
         }
      }
   }
   Fout->close();
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

   if(!this->getAllocated()){
      rs_error("Image2D::writeImage: Image is not allocated. ");
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
   if(nhx == 0) nhx = 1;
   if(nhz == 0) nhz = 1;

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
std::shared_ptr<rockseis::Image2D<T>> Image2D<T>::getLocal(std::shared_ptr<rockseis::Data2D<T>> data, T aperture, bool map) 
{
   std::shared_ptr<rockseis::Image2D<T>> local;
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
   size_t nz = this->getNz();
   size_t nx = this->getNx();
   size_t size;
   off_t start;

   size_t nhx = this->getNhx();
   size_t nhz = this->getNhz();

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
         size = (size_t) (rintf((max-min + aperture)/dx) + 1);
         if( size % 2 == 0 ) size++; // Get odd size due to symmetry
      }else{
         size = (size_t) (2*rintf(off/dx) + 1);
      }
      start = (off_t) (rintf((min - ox)/dx) - (size - 1)/2); 
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
      size = (size_t) (rintf((max-min + 2*fabs(daperture))/dx) + 2);
      start = (off_t) (rintf((min - ox)/dx) - rintf(fabs(daperture/dx))) - 1; 
   }

   /* Create local model */
   local = std::make_shared<rockseis::Image2D<T>>(size, this->getNz(), this->getNhx(), this->getNhz(), dx, this->getDz(), (ox + start*dx), this->getOz());

   /*Realizing local model */
   local->allocateImage();

   /* Copying from big model into local model */
   T *Im = local->getImagedata();

   /* Allocate two traces to read models from file */
   T *imtrace = (T *) calloc(nx, sizeof(T));
   if(imtrace == NULL) rs_error("Image2D::getlocal: Failed to allocate memory.");

   // Open files for reading
   bool status;
   std::shared_ptr<rockseis::File> Fim (new rockseis::File());
   status = Fim->input(imagefile);
   if(status == FILE_ERR){
      rs_error("Image2D::getLocal : Error reading from image file.");
   }

   off_t i = start;
   off_t lpos, fpos;
   rockseis::Index l2d(size,nz,nhx,nhz);
   rockseis::Index f2d(nx,nz,nhx,nhz);
   for(long int ih1=0; ih1<nhz; ih1++) {
      for(long int ih2=0; ih2<nhx; ih2++) {
         for(long int i1=0; i1<nz; i1++) {
            fpos = f2d(0, i1, ih2, ih1)*sizeof(T);
            Fim->read(imtrace, nx, fpos);
            if(Fim->getFail()) rs_error("Image2D::getLocal: Error reading from image file");
            for(long int i2=0; i2<size; i2++) {
               lpos = i + i2;
               if(lpos < 0) lpos = 0;
               if(lpos > (nx-1)) lpos = nx - 1;
               Im[l2d(i2,i1,ih2,ih1)] = imtrace[lpos];
            }
         }
      }
   }

   /* Free traces */
   free(imtrace);

   return local;
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


   if(_nhx == 0) _nhx = 1;
   if(_nhy == 0) _nhy = 1;
   if(_nhz == 0) _nhz = 1;

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
   incore = false;
}

template<typename T>
Image3D<T>::Image3D(std::string _imagefile, std::shared_ptr<ModelEikonal3D<T>> model, int nhx, int nhy, int nhz): Model<T>(3)
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
   incore = false;
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
   incore = false;
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
   incore = false;
}

template<typename T>
Image3D<T>::Image3D(const int _nx, const int _ny, const int _nz, const int _nhx, const int _nhy, const int _nhz, const T _dx, const T _dy, const T _dz, const T _ox, const T _oy, const T _oz)
{
   this->setNx(_nx);
   this->setNy(_ny);
   this->setNz(_nz);
   this->setDx(_dx);
   this->setDy(_dy);
   this->setDz(_dz);
   this->setOx(_ox);
   this->setOy(_oy);
   this->setOz(_oz);
   this->setNhx(_nhx);
   this->setNhy(_nhy);
   this->setNhz(_nhz);
   allocated = false;
   incore = false;
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

   if(!this->getAllocated()){
      rs_error("Image3D::writeImage: Image is not allocated. ");
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

   if(nhx == 0) nhx = 1;
   if(nhy == 0) nhy = 1;
   if(nhz == 0) nhz = 1;

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
bool Image3D<T>::stackImage(std::shared_ptr<Image3D<T>>  imagein)
{
   if(!imagein->getAllocated()) rs_error("Image3D<T>::stackImage: Input image data not allocated.");
   if(!this->getAllocated()) rs_error("Image3D<T>::stackImage: Output image data not allocated.");
   bool status;
   long int nxg, nyg,nzg;
   int nhxg, nhyg, nhzg;
   T *dataout=this->getImagedata();
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

   T *datain=imagein->getImagedata();
   nxl = imagein->getNx();
   nyl = imagein->getNy();
   nzl = imagein->getNz();
   nhxl = imagein->getNhx();
   nhyl = imagein->getNhy();
   nhzl = imagein->getNhz();
   dxl = imagein->getDx();
   dyl = imagein->getDy();
   dzl = imagein->getDz();
   oxl = imagein->getOx();
   oyl = imagein->getOy();
   ozl = imagein->getOz();

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
                        for(ix=0; ix < nxl; ix++){
                           trcin[ix] = datain[Iin(ix,iy,iz,ihx,ihy,ihz)];
                        }

                        for(ix=0; ix < nxg; ix++){
                           trcout[ix] = dataout[Iout(ix, (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)];
                        }
                        for(ix = 0; ix < nxl; ix++) {
                           if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
                              trcout[ix + ix_start] += trcin[ix];
                           }
                        }
                        // Write trc 
                        for(ix=0; ix < nxg; ix++){
                           dataout[Iout(ix, (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)] = trcout[ix];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   free(trcin);
   free(trcout);
   status = FILE_OK;
   return status;
}


template<typename T>
bool Image3D<T>::stackImage_parallel(std::string infile,int padlx, int padhx, int padly, int padhy, int padlz, int padhz) 
{
   bool status;
   if(infile.empty()){
      rs_error("Image3D::stackImage_parallel: Input filename is empty. ");
   }

   std::shared_ptr<rockseis::File> Fout (new rockseis::File());
   Fout->append(infile.c_str());
   if ( !Fout->is_open()) rs_error("Image3D::stackImage_parallel: Failed to open output image: ", infile);

   long int nxg, nyg,nzg;
   int nhxg, nhyg, nhzg;
   T dxg, dyg, dzg, oxg, oyg, ozg;
   nxg = Fout->getN(1);
   nyg = Fout->getN(2);
   nzg = Fout->getN(3);
   nhxg = Fout->getN(4);
   nhyg = Fout->getN(5);
   nhzg = Fout->getN(6);
   dxg = Fout->getD(1);
   dyg = Fout->getD(2);
   dzg = Fout->getD(3);
   oxg = Fout->getO(1);
   oyg = Fout->getO(2);
   ozg = Fout->getO(3);

   long int nxl, nyl, nzl;
   int nhxl, nhyl, nhzl;
   T dxl, dyl, dzl, oxl, oyl, ozl;
   nxl = this->getNx();
   nyl = this->getNy();
   nzl = this->getNz();
   nhxl = this->getNhx();
   nhyl = this->getNhy();
   nhzl = this->getNhz();
   dxl = this->getDx();
   dyl = this->getDy();
   dzl = this->getDz();
   oxl = this->getOx();
   oyl = this->getOy();
   ozl = this->getOz();

   if(nhxg != nhxl || nhyg != nhyl || nhzg != nhzl || dxg != dxl || dyg != dyl || dzg != dzl){
      rs_error("Image3D::stackImage_parallel: Images are not compatible. Cannot stack.");
   }

   T *trcin = this->getImagedata();
   T trcout;

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
            for(iz=padlz; iz<nzl-padhz; iz++) {
               if((iz + iz_start) >= 0 && (iz + iz_start) < nzg){
                  for(iy=padly; iy<nyl-padhy; iy++) {
                     if((iy + iy_start) >= 0 && (iy + iy_start) < nyg){
                        for(ix = padlx; ix < nxl-padhx; ix++) {
                           if((ix + ix_start) >= 0 && (ix + ix_start) < nxg){
                              // Read traces 
                              Fout->read(&trcout, 1, Iout((ix+ix_start), (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)*sizeof(T));
                              if(Fout->getFail()) rs_error("Image3D::stackImage_parallel: Failed to write data to file");
                              trcout += trcin[Iin(ix,iy,iz,ihx,ihy,ihz)];
                              // Write trc 
                              Fout->write(&trcout, 1, Iout((ix+ix_start), (iy+iy_start), (iz+iz_start),ihx,ihy,ihz)*sizeof(T));
                              if(Fout->getFail()) rs_error("Image3D::stackImage_parallel: Failed to write data to file");
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   Fout->close();
   status = FILE_OK;
   return status;
}

template<typename T>
std::shared_ptr<rockseis::Image3D<T>> Image3D<T>::getLocal(std::shared_ptr<rockseis::Data3D<T>> data, T aperture_x, T aperture_y, bool map) {

   std::shared_ptr<rockseis::Image3D<T>> local;
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
   size_t nx = this->getNx();
   size_t ny = this->getNy();
   size_t nz = this->getNz();

   size_t size_x;
   off_t start_x;
   size_t size_y;
   off_t start_y;

   size_t nhx = this->getNhx();
   size_t nhy = this->getNhy();
   size_t nhz = this->getNhz();

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
         size_x = (size_t) (rintf((max_x-min_x + aperture_x)/dx) + 1);
         if( size_x % 2 == 0 ) size_x++; // Get odd size due to symmetry
      }else{
         size_x = (size_t) (2*rintf(off_x/dx) + 1);
      }
      start_x = (off_t) (rintf((min_x - ox)/dx) - (size_x - 1)/2); 
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
      size_x = (size_t) (rintf((max_x-min_x + 2*fabs(daperture_x))/dx) + 2);
      start_x = (off_t) (rintf((min_x - ox)/dx) - rintf(fabs(daperture_x/dx))) - 1; 
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
         size_y = (size_t) (rintf((max_y-min_y + aperture_y)/dy) + 1);
         if( size_y % 2 == 0 ) size_y++; // Get odd size due to symmetry
      }else{
         size_y = (size_t) (2*rintf(off_y/dy) + 1);
      }
      start_y = (off_t) (rintf((min_y - oy)/dy) - (size_y - 1)/2); 
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
      size_y = (size_t) (rintf((max_y-min_y + 2*fabs(daperture_y))/dy) + 2);
      start_y = (off_t) (rintf((min_y - oy)/dy) - rintf(fabs(daperture_y/dy))) - 1; 
   }

   double oxl, oyl; 
   oxl = (ox + start_x*dx);
   oyl = (oy + start_y*dy);

   /* Create local model */
   local = std::make_shared<rockseis::Image3D<T>>(size_x, size_y, nz, nhx, nhy, nhz, dx, dy, this->getDz(), oxl, oyl, this->getOz());

   /*Realizing local model */
   local->allocateImage();

   /* Copying from big model into local model */
   T *Im = local->getImagedata();

   /* Allocate two traces to read models from file */
   T *imtrace = (T *) calloc(nx*ny, sizeof(T));
   if(imtrace == NULL) rs_error("Image3D::getlocal: Failed to allocate memory.");

   // Open files for reading
   bool status;
   std::shared_ptr<rockseis::File> Fim (new rockseis::File());
   status = Fim->input(imagefile);
   if(status == FILE_ERR){
      rs_error("Image3D::getLocal : Error reading from image file.");
   }

   off_t i = start_x;
   off_t j = start_y;
   off_t lpos_x, lpos_y, fpos;
   rockseis::Index l3d(size_x, size_y, nz, nhx, nhy, nhz);
   rockseis::Index f3d(nx, ny, nz, nhx, nhy, nhz);
   rockseis::Index l2d(nx, ny);
   for(size_t ih1=0; ih1<nhz; ih1++) {
      for(size_t ih2=0; ih2<nhy; ih2++) {
         for(size_t ih3=0; ih3<nhx; ih3++) {
            for(size_t i1=0; i1<nz; i1++) {
               fpos = f3d(0, 0, i1, ih3, ih2, ih1)*sizeof(T);
               Fim->read(imtrace, nx*ny, fpos);
               if(Fim->getFail()) rs_error("Image3D::getLocal: Error reading from image file");
               for(size_t i2=0; i2<size_y; i2++) {
                  lpos_y = j + i2;
                  if(lpos_y < 0) lpos_y = 0;
                  if(lpos_y > (ny-1)) lpos_y = ny - 1;
                  for(size_t i3=0; i3<size_x; i3++) {
                     lpos_x = i + i3;
                     if(lpos_x < 0) lpos_x = 0;
                     if(lpos_x > (nx-1)) lpos_x = nx - 1;
                     Im[l3d(i3,i2,i1,ih3,ih2,ih1)] = imtrace[l2d(lpos_x, lpos_y)];
                  }
               }
            }
         }
      }
   }

   /* Free traces */
   free(imtrace);

   return local;
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


