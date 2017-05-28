#include "snap.h"

namespace rockseis {
// constructor
template<typename T>
Snapshot<T>::Snapshot()
{
    for(int i=0; i < NPTR; i++) allocated[i] = false;
    open = false;
    field = PRESSURE;
    geometry = std::make_shared<Geometry<T>>(); 
    Fp = std::make_shared<File>();
    snapit = 0;
}

//Open Snap
template<typename T>
bool Snapshot<T>::openSnap(std::string filename, char flag) {
    if(!filename.empty()){
        switch(flag){
            case 'w':
                this->filename = filename;
                this->Fp->output(this->filename);
                this->open = true;
                this->Fp->setN(1,this->getNx());
                this->Fp->setN(2,this->getNy());
                this->Fp->setN(3,this->getNz());
                this->Fp->setN(4,this->getSnapnt());
                this->Fp->setD(1,this->getDx());
                this->Fp->setD(2,this->getDy());
                this->Fp->setD(3,this->getDz());
                this->Fp->setD(4,this->getSnapdt());
                this->Fp->setO(1,this->getOx());
                this->Fp->setO(2,this->getOy());
                this->Fp->setO(3,this->getOz());
                this->Fp->setO(4,this->getSnapot());
                this->Fp->setData_format(sizeof(T));
                this->Fp->setType(rockseis::SNAPSHOT);
                this->Fp->writeHeader();
                this->Fp->seekp(this->Fp->getStartofdata());
                this->setSnapit(0);
                break;
            case 'r':
                if(this->open) rs_error("Snapshot cannot be opened two times.");
                this->filename = filename;
                if(this->Fp->input(this->filename) == FILE_ERR)
                {
                    rs_error("Snap::openSnap: Error opening snapshot file for reading.");
                }
                if(this->getNx() != this->Fp->getN(1)) rs_error("Snap::openSnap: Mismatch in nx size of snaps");
                if(this->getNy() != this->Fp->getN(2)) rs_error("Snap::openSnap: Mismatch in ny size of snaps");
                if(this->getNz() != this->Fp->getN(3)) rs_error("Snap::openSnap: Mismatch in nz size of snaps");
                if(sizeof(T) != this->Fp->getData_format()) rs_error("Snap::openSnap: Mismatch in precision of snaps");
                if(this->getSnapnt() != this->Fp->getN(4)) rs_error("Snap::openSnap: Mismatch in number of snaps");
                if(this->Fp->getType() != rockseis::SNAPSHOT) rs_error("Snap::openSnap: Mismatch in file type");
                this->open = true;
                this->setSnapit(this->getSnapnt() - 1);
                break;
            case 'a':
                if(this->open) rs_error("Snapshot cannot be opened two times.");
                if(this->Fp->append(this->filename) == FILE_ERR)
                {
                    rs_error("Snap::openSnap: Error opening snapshot file for reading.");
                }
                if(this->getNx() != this->Fp->getN(1)) rs_error("Snap::openSnap: Mismatch in nx size of snaps");
                if(this->getNy() != this->Fp->getN(2)) rs_error("Snap::openSnap: Mismatch in ny size of snaps");
                if(this->getNz() != this->Fp->getN(3)) rs_error("Snap::openSnap: Mismatch in nz size of snaps");
                if(sizeof(T) != this->Fp->getData_format()) rs_error("Snap::openSnap: Mismatch in precision of snaps");
                if(this->getSnapnt() != this->Fp->getN(4)) rs_error("Snap::openSnap: Mismatch in number of snaps");
                if(this->Fp->getType() != rockseis::SNAPSHOT) rs_error("Snap::openSnap: Mismatch in file type");
                this->open = true;
                this->setSnapit(this->getSnapnt() - 1);
                break;
            default:
                rs_error("Snap::openSnap: Invalid flag.");
                break;
        }
    }else{
        rs_error("Snap::openSnap: No filename set.");
    }
    return SNAP_OK;
}

//Open Edge
template<typename T>
bool Snapshot<T>::openEdge(std::string filename, char flag) {
    if(!filename.empty()){
    int count = 0;
        switch(flag){
            case 'w':
                if(this->open) rs_error("Snapshot cannot be opened two times.");
                this->filename = filename;
                this->Fp->output(this->filename);
                this->open = true;
                for(int i=0; i<NPTR; i++){
                    if(data[i] != NULL) count++;
                }
                if( count==0 ) rs_error("Snapshot::openEdge: No field set for snapping");
                if(this->getDim() ==2){
                    this->Fp->setN(1,2*this->getNx() + 2*this->getNz());
                }else{
                    this->Fp->setN(1,2*this->getNx()*this->getNz() + 2*this->getNy()*this->getNx() + 2*this->getNz()*this->getNy());
                }
                this->Fp->setN(2,count);
                this->Fp->setN(3,1);
                this->Fp->setN(4,this->getSnapnt());
                this->Fp->setD(1,this->getDx());
                this->Fp->setD(2,this->getDy());
                this->Fp->setD(3,this->getDz());
                this->Fp->setD(4,this->getSnapdt());
                this->Fp->setO(1,this->getOx());
                this->Fp->setO(2,this->getOy());
                this->Fp->setO(3,this->getOz());
                this->Fp->setO(4,this->getSnapot());
                this->Fp->setData_format(sizeof(T));
                this->Fp->setType(rockseis::EDGESNAP);
                this->Fp->writeHeader();
                this->Fp->seekp(this->Fp->getStartofdata());
                this->setSnapit(0);
                break;
            case 'r':
                if(this->open) rs_error("Snapshot cannot be opened two times.");
                this->filename = filename;
                if(this->Fp->input(this->filename) == FILE_ERR)
                {
                    rs_error("Snap::openEdge: Error opening snapshot file for reading.");
                }
                if( this->getDim() == 2 ) {
                if((2*this->getNx() + 2*this->getNz()) != this->Fp->getN(1)) 
                    rs_error("Snap::openEdge: Mismatch in size of edge snaps");
                }else{
                if((2*this->getNx()*this->getNz() + 2*this->getNy()*this->getNx() + 2*this->getNz()*this->getNy()) != this->Fp->getN(1)) 
                    rs_error("Snap::openEdge: Mismatch in nx size of snaps");
                }
                if(sizeof(T) != this->Fp->getData_format()) rs_error("Snap::openEdge: Mismatch in precision of snaps");
                if(this->getSnapnt() != this->Fp->getN(4)) rs_error("Snap::openEdge: Mismatch in number of snaps");
                if(this->Fp->getType() != rockseis::EDGESNAP) rs_error("Snap::openEdge: Mismatch in file type");
                this->open = true;
                this->setSnapit(this->getSnapnt() - 1);
                break;
            case 'a':
                if(this->open) rs_error("Snapshot cannot be opened two times.");
                this->filename = filename;
                if(this->Fp->append(this->filename) == FILE_ERR)
                {
                    rs_error("Snap::openEdge: Error opening snapshot file for reading.");
                }
                if( this->getDim() == 2 ) {
                    if((2*this->getNx() + 2*this->getNz()) != this->Fp->getN(1)) 
                        rs_error("Snap::openEdge: Mismatch in size of edge snaps");
                }else{
                    if((2*this->getNx()*this->getNz() + 2*this->getNy()*this->getNx() + 2*this->getNz()*this->getNy()) != this->Fp->getN(1)) 
                        rs_error("Snap::openEdge: Mismatch in nx size of snaps");
                }
                if(sizeof(T) != this->Fp->getData_format()) rs_error("Snap::openEdge: Mismatch in precision of snaps");
                if(this->getSnapnt() != this->Fp->getN(4)) rs_error("Snap::openEdge: Mismatch in number of snaps");
                if(this->Fp->getType() != rockseis::EDGESNAP) rs_error("Snap::openEdge: Mismatch in file type");
                this->open = true;
                this->setSnapit(this->getSnapnt() - 1);
                break;
            default:
                rs_error("Snap::openEdge: Invalid flag.");
                break;
        }
    }else{
        rs_error("Snap::openEdge: No filename set.");
    }
    return SNAP_OK;
}

template<typename T>
void Snapshot<T>::closeSnap() {
    if(this->open) {
        this->Fp->close();
        this->open = false;
    }
}

template<typename T>
void Snapshot<T>::removeSnap() {
    if(this->open) {
        this->Fp->close();
        this->open = false;
    }
	if(!this->filename.empty()){
		if( remove( filename.c_str() ) != 0 ){
			rs_error( "Snapshot::removeSnap: Error deleting file: ", filename);
		}
	}
}

template<typename T>
void Snapshot<T>::allocSnap(int i) 
{
    if(i >= NPTR || i < 0) rs_error("Snap::allocSnap: Trying to allocate data out of bounds.");
    if(this->allocated[i] == false) {
        this->data[i] = (T *) calloc(this->getNx()*this->getNy()*this->getNz(), sizeof(T));
        this->allocated[i] = true;
    }else{
        free(this->data[i]);
        this->data[i] = (T *) calloc(this->getNx()*this->getNy()*this->getNz(), sizeof(T));
        this->allocated[i] = true;
    }
}

template<typename T>
void Snapshot<T>::freeSnaps()
{
    for(int i=0; i<NPTR; i++){
        if(this->allocated[i]) free(this->data[i]);
    }
}


template<typename T>
Snapshot<T>::~Snapshot() {
    this->closeSnap();
    for(int i=0; i<NPTR; i++){
        if(this->allocated[i]) free(this->data[i]);
    }
}


// =============== 2D SNAPSHOT CLASS =============== //

template<typename T>
Snapshot2D<T>::Snapshot2D(std::shared_ptr<WavesAcoustic2D<T>> waves, int snapinc)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _nt;
    T _dt, _ot;
    int _lpml;
    int _dim;
    int snapnt, enddiff;
    T snapdt, snapot;

/* Get necessary parameters from waves class */
    _nx=waves->getNx();
    _ny=waves->getNy();
    _nz=waves->getNz();
    _dx=waves->getDx();
    _dy=waves->getDy();
    _dz=waves->getDz();
    _ox=waves->getOx();
    _oy=waves->getOy();
    _oz=waves->getOz();
    _lpml = waves->getLpml();
    _dim = waves->getDim();
    _nt = waves->getNt();
    _dt = waves->getDt();
    _ot = waves->getOt();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setLpml(_lpml);
    this->setDim(_dim);
    this->setSnapinc(snapinc);

    snapnt = (_nt-1)/snapinc + 1;
    snapdt = _dt*snapinc;
    snapot = _ot;
    enddiff = (int) rintf(((_nt-1)*_dt - (snapnt-1)*snapdt)/_dt);
    this->setSnapnt(snapnt);
    this->setSnapdt(snapdt);
    this->setSnapot(snapot);
    this->setEnddiff(enddiff);
    for(int i=0; i<NPTR; i++){
        this->setData(NULL, i);
    }
}

template<typename T>
Snapshot2D<T>::Snapshot2D(std::shared_ptr<WavesElastic2D<T>> waves, int snapinc)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _nt;
    T _dt, _ot;
    int _lpml;
    int _dim;
    int snapnt, enddiff;
    T snapdt, snapot;

/* Get necessary parameters from waves class */
    _nx=waves->getNx();
    _ny=waves->getNy();
    _nz=waves->getNz();
    _dx=waves->getDx();
    _dy=waves->getDy();
    _dz=waves->getDz();
    _ox=waves->getOx();
    _oy=waves->getOy();
    _oz=waves->getOz();
    _lpml = waves->getLpml();
    _dim = waves->getDim();
    _nt = waves->getNt();
    _dt = waves->getDt();
    _ot = waves->getOt();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setLpml(_lpml);
    this->setDim(_dim);
    this->setSnapinc(snapinc);

    snapnt = (_nt-1)/snapinc + 1;
    snapdt = _dt*snapinc;
    snapot = _ot;
    enddiff = (int) rintf(((_nt-1)*_dt - (snapnt-1)*snapdt)/_dt);
    this->setSnapnt(snapnt);
    this->setSnapdt(snapdt);
    this->setSnapot(snapot);
    this->setEnddiff(enddiff);
    for(int i=0; i<NPTR; i++){
        this->setData(NULL, i);
    }
}

// Write Snapshots
template<typename T>
void Snapshot2D<T>::writeSnap(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    int nx_pml = this->getNx_pml();
    int nz_pml = this->getNz_pml();
    int lpml = this->getLpml();
    Index I(nx_pml,nz_pml);
    int i,j;
    int snapit = this->getSnapit();
    T *data1 = this->getData(0);
    T *data2 = this->getData(1);
    T *data3 = this->getData(2);
    std::shared_ptr<rockseis::File> Fp = this->getFp();
    T val = 0;
    if(this->getOpen()){
        if((it % this->getSnapinc()) == 0){
            this->setSnapit(snapit + 1); // Increment snap counter
            //Write snapshot
            for(j=0; j<nz; j++){
                for(i=0; i<nx; i++){
                    if(data1 != NULL) val = data1[I(i+lpml,j+lpml)];
                    if(data2 != NULL) val += data2[I(i+lpml,j+lpml)];
                    if(data3 != NULL) val += data3[I(i+lpml,j+lpml)];
                    Fp->write(&val,1); 
                }
            }
        }
    }
    if(Fp->getFail()) rs_error("Snapshot2D::writeSnap: Error writting to file.");
}

// read Snapshots
template<typename T>
void Snapshot2D<T>::readSnap(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    int snapit = this->getSnapit();
    std::shared_ptr<rockseis::File> Fp = this->getFp();
    off_t pos;

    if(this->getOpen() && this->getAllocated(0)){
        if(((it-this->getEnddiff()) % this->getSnapinc()) == 0){
            this->setSnapit(snapit - 1); // Increment snap counter
            //Read snapshot
            pos  = nx*nz*snapit*sizeof(T);
            Fp->read(this->getData(0), nz*nx, pos); 
        }
    }
    if(Fp->getFail()) rs_error("Snapshot2D::readSnap: Error reading from file.");
}

// Write Edges
template<typename T>
void Snapshot2D<T>::writeEdge(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    int nx_pml = this->getNx_pml();
    int nz_pml = this->getNz_pml();
    int lpml = this->getLpml();
    Index I(nx_pml,nz_pml);
    int i,j;
    int snapit = this->getSnapit();
    T *data[NPTR];
    for(i=0; i<NPTR; i++) data[i] = this->getData(i);
    std::shared_ptr<rockseis::File> Fp = this->getFp();
    T val = 0;
    if(this->getOpen()){
        if((it % this->getSnapinc()) == 0){
            this->setSnapit(snapit + 1); // Increment snap counter
            //Write snapshot
            for(i=0; i<NPTR; i++)
            {
                if(data[i] != NULL) 
                {
                    //Left and right
                    for(j=0; j<nz; j++)
                    {
                        val = data[i][I(lpml,j+lpml)]; // Left
                        Fp->write(&val,1); 
                        val = data[i][I(nx+lpml,j+lpml)]; //Right
                        Fp->write(&val,1); 
                    }
                    //Top and bottom
                    for(j=0; j<nx; j++)
                    {
                        val = data[i][I(j+lpml,lpml)]; // Top
                        Fp->write(&val,1); 
                        val = data[i][I(j+lpml,nz+lpml)]; //Bottom
                        Fp->write(&val,1); 
                    }

                }
            }
        }
    }
    if(Fp->getFail()) rs_error("Snapshot2D::writeEdge: Error writting to file.");
}

// read Edges
template<typename T>
void Snapshot2D<T>::readEdge(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    int nx_pml = this->getNx_pml();
    int nz_pml = this->getNz_pml();
    int lpml = this->getLpml();
    int snapit = this->getSnapit();
    std::shared_ptr<rockseis::File> Fp = this->getFp();
    off_t pos;

    int i,j;
    Index I(nx_pml,nz_pml);
    T *data[NPTR];
    for(i=0; i<NPTR; i++) data[i] = this->getData(i);

    T val = 0;
    if(this->getOpen()){
        if(((it-this->getEnddiff()) % this->getSnapinc()) == 0){
            this->setSnapit(snapit - 1); // Increment snap counter
            pos  = (2*nx + 2*nz)*snapit*sizeof(T);
            Fp->seekg(Fp->getStartofdata() + pos);
            for(i=0; i<NPTR; i++)
            {
                if(data[i] != NULL) 
                {
                    //Left and right
                    for(j=0; j<nz; j++)
                    {
                        Fp->read(&val,1); 
                        data[i][I(lpml,j+lpml)] = val; // Left
                        Fp->read(&val,1); 
                        data[i][I(nx+lpml,j+lpml)] = val; //Right
                    }
                    //Top and bottom
                    for(j=0; j<nx; j++)
                    {
                        Fp->read(&val,1); 
                        data[i][I(j+lpml,lpml)] = val; // Top
                        Fp->read(&val,1); 
                        data[i][I(j+lpml,nz+lpml)] = val; //Bottom
                    }
                }
            }
        }
    }
    if(Fp->getFail()) rs_error("Snapshot2D::readEdge: Error reading from file.");
}

template<typename T>
Snapshot2D<T>::~Snapshot2D() 
{
    // Do nothing.
}

// =============== 3D SNAPSHOT CLASS =============== //

template<typename T>
Snapshot3D<T>::Snapshot3D(std::shared_ptr<WavesAcoustic3D<T>> waves, int snapinc)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _nt;
    T _dt, _ot;
    int _lpml;
    int _dim;
    int snapnt, enddiff;
    T snapdt, snapot;

/* Get necessary parameters from waves class */
    _nx=waves->getNx();
    _ny=waves->getNy();
    _nz=waves->getNz();
    _dx=waves->getDx();
    _dy=waves->getDy();
    _dz=waves->getDz();
    _ox=waves->getOx();
    _oy=waves->getOy();
    _oz=waves->getOz();
    _lpml = waves->getLpml();
    _nt = waves->getNt();
    _dt = waves->getDt();
    _ot = waves->getOt();
    _dim = waves->getDim();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setLpml(_lpml);
    this->setDim(_dim);
    this->setSnapinc(snapinc);

    snapnt = (_nt-1)/snapinc + 1;
    snapdt = _dt*snapinc;
    snapot = _ot;
    enddiff = (int) rintf(((_nt-1)*_dt - (snapnt-1)*snapdt)/_dt);
    this->setSnapnt(snapnt);
    this->setSnapdt(snapdt);
    this->setSnapot(snapot);
    this->setEnddiff(enddiff);
    for(int i=0; i<NPTR; i++){
        this->setData(NULL, i);
    }
}

template<typename T>
Snapshot3D<T>::Snapshot3D(std::shared_ptr<WavesElastic3D<T>> waves, int snapinc)
{
    int _nx, _ny, _nz;
    T _dx, _dy, _dz; 
    T _ox, _oy, _oz; 
    int _nt;
    T _dt, _ot;
    int _lpml;
    int _dim;
    int snapnt, enddiff;
    T snapdt, snapot;

/* Get necessary parameters from waves class */
    _nx=waves->getNx();
    _ny=waves->getNy();
    _nz=waves->getNz();
    _dx=waves->getDx();
    _dy=waves->getDy();
    _dz=waves->getDz();
    _ox=waves->getOx();
    _oy=waves->getOy();
    _oz=waves->getOz();
    _lpml = waves->getLpml();
    _nt = waves->getNt();
    _dt = waves->getDt();
    _ot = waves->getOt();
    _dim = waves->getDim();

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setLpml(_lpml);
    this->setSnapinc(snapinc);
    this->setDim(_dim);

    snapnt = (_nt-1)/snapinc + 1;
    snapdt = _dt*snapinc;
    snapot = _ot;
    enddiff = (int) rintf(((_nt-1)*_dt - (snapnt-1)*snapdt)/_dt);
    this->setSnapnt(snapnt);
    this->setSnapdt(snapdt);
    this->setSnapot(snapot);
    this->setEnddiff(enddiff);
    for(int i=0; i<NPTR; i++){
        this->setData(NULL, i);
    }
}

    // Write Snapshots
template<typename T>
void Snapshot3D<T>::writeSnap(const int it){
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();
    int nx_pml = this->getNx_pml();
    int ny_pml = this->getNy_pml();
    int nz_pml = this->getNz_pml();
    int lpml = this->getLpml();
    Index I(nx_pml,ny_pml,nz_pml);
    int i,j,k;
    int snapit = this->getSnapit();
    T *data1 = this->getData(0);
    T *data2 = this->getData(1);
    T *data3 = this->getData(2);
    std::shared_ptr<rockseis::File> Fp = this->getFp();
    T val = 0;
    if(this->getOpen()){
        if((it % this->getSnapinc()) == 0){
            this->setSnapit(snapit + 1); // Increment snap counter
            //Write snapshot
            for(k=0; k<nz; k++){
                for(j=0; j<ny; j++){
                    for(i=0; i<nx; i++){
                        if(data1 != NULL) val = data1[I(i+lpml,j+lpml,k+lpml)];
                        if(data2 != NULL) val += data2[I(i+lpml,j+lpml,k+lpml)];
                        if(data3 != NULL) val += data3[I(i+lpml,j+lpml,k+lpml)];
                        Fp->write(&val,1); 
                    }
                }
            }
        }
    }
    if(Fp->getFail()) rs_error("Snapshot3D::writeSnap: Error writting to file.");
}

// read Snapshots
template<typename T>
void Snapshot3D<T>::readSnap(const int it){
    int nx = this->getNx();
    int ny = this->getNy();
    int nz = this->getNz();
    int snapit = this->getSnapit();
    off_t pos;
    std::shared_ptr<rockseis::File> Fp = this->getFp();

    if(this->getOpen() && this->getAllocated(0)){
        if(((it-this->getEnddiff()) % this->getSnapinc()) == 0){
            this->setSnapit(snapit - 1); // Increment snap counter
            //Read snapshot
            pos = nx*ny*nz*snapit*sizeof(T);
            Fp->read(this->getData(0),nz*ny*nx, pos); 
        }
    }
}

template<typename T>
Snapshot3D<T>::~Snapshot3D() 
{
    // Do nothing.
}




// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Snapshot<float>;
template class Snapshot<double>;
template class Snapshot2D<float>;
template class Snapshot2D<double>;
template class Snapshot3D<float>;
template class Snapshot3D<double>;

}
