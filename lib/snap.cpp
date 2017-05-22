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
    snapinc=1;
    enddiff=0;
    
}

// Make constructors for different Wave classes 
// fix Write snap to distinguish between 2D and 3D dimension 
// And also between Pressure and other Snaps

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
                this->Fp->setN(3,this->getNz());
                this->Fp->setN(4,this->getSnapnt());
                this->Fp->setD(1,this->getDx());
                this->Fp->setD(3,this->getDz());
                this->Fp->setD(4,this->getSnapdt());
                this->Fp->setO(1,this->getOx());
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

// Write Snapshots
template<typename T>
void Snapshot<T>::writeSnap(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    int nx_pml = this->getNx_pml();
    int nz_pml = this->getNz_pml();
    int lpml = this->getLpml();
    Index I(nx_pml,nz_pml);
    int i,j;
    int snapit = this->getSnapit();

    if(this->open){
       if((it % this->getSnapinc()) == 0){
           this->setSnapit(snapit + 1); // Increment snap counter
           //Write snapshot
           for(j=0; j<nz; j++){
               for(i=0; i<nx; i++){
                   Fp->write(&data[0][I(i+lpml,j+lpml)],1); 
               }
           }
        }
    }
}

// Write Snapshots
template<typename T>
void Snapshot<T>::readSnap(const int it){
    int nx = this->getNx();
    int nz = this->getNz();
    Index I(nx,nz);
    int snapit = this->getSnapit();

    if(this->open && this->allocated[0]){
       if(((it-this->getEnddiff()) % this->getSnapinc()) == 0){
           this->setSnapit(snapit - 1); // Increment snap counter
           //Read snapshot
           this->Fp->seekp(nx*nz*snapit*sizeof(T));
           this->Fp->read(this->data[0],nz*nx); 
       }
    }
}


template<typename T>
void Snapshot<T>::closeSnap() {
    if(this->open) {
        this->Fp->close();
        this->open = false;
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
    for(int i=0; i<NPTR; i++){
        if(this->allocated[i]) free(this->data[i]);
    }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Snapshot<float>;
template class Snapshot<double>;

}
