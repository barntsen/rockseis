// Include statements
#include "snap.h"

namespace rockseis {

// =============== ABSTRACT WAVES CLASS =============== //
template<typename T>
Snap<T>::Snap() {
    geometry = std::make_shared<Geometry<T>>(); 
    geometry->setN(1, 1);
    geometry->setN(2, 1);
    geometry->setN(3, 1);
    geometry->setN(4, 1);
    geometry->setD(1, 1.);
    geometry->setD(2, 1.);
    geometry->setD(3, 1.);
    geometry->setD(4, 1.);
    geometry->setO(1, 0.);
    geometry->setO(2, 0.);
    geometry->setO(3, 0.);
    geometry->setO(4, 0.);
    lpml = 10;
}

template<typename T>
Snap<T>::~Snap() {
    // Nothing here
}

// =============== 2D ACOUSTIC MODEL CLASS =============== //
template<typename T>
SnapAcoustic2D<T>::SnapAcoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, const int _nt_mod, T _dt_mod, T _ot_mod, const int snap_inc): Snap<T>(){
    int _nx, _ny, _nz, _nt;
    T _dx, _dy, _dz, _dt;
    T _ox, _oy, _oz, _ot;
    int _lpml;
    int end_diff;

    /* Get necessary parameters from model class */
    _nx = model->getNx();
    _ny = model->getNy();
    _nz = model->getNz();
    _dx = model->getDx();
    _dy = model->getDy();
    _dz = model->getDz();
    _ox = model->getOx();
    _oy = model->getOy();
    _oz = model->getOz();
    _lpml = model->getLpml();

    if(snap_inc < 1) {
	    /* Error */ 
	    std::cerr << "SnapAcoustic2D::snap_inc must be an integer larger than 0.\n";
	    exit(1);
    }
    if(snap_inc > _nt_mod) {
	    /* Error */ 
	    std::cerr << "SnapAcoustic2D::snap_inc is larger than modelling nt.\n";
	    exit(1);
    }

    // Setting snap time axis
    _nt = _nt_mod/snap_inc;
    _dt = _dt_mod*snap_inc;
    _ot = _ot_mod;
    end_diff = (_nt_mod*_dt_mod - _nt*_dt)/_dt_mod;
    this->setEnd_diff(end_diff);

    this->setNx(_nx);
    this->setNy(_ny);
    this->setNz(_nz);
    this->setNt(_nt);
    this->setNt_mod(_nt_mod);
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setDt_mod(_dt_mod);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setOt_mod(_ot_mod);
    this->setLpml(_lpml);

    // Setting bools
    Ax = 0;
    Az = 0;
    P = 0;
}

template<typename T>
void SnapAcoustic2D<T>::putPfile(std::string filename) {
    if(!filename.empty()){
        Pfile = filename;
        FP = std::make_shared<File>();
        FP->output(Pfile);
        P = true;
        FP->setN(1,this->getNx());
        FP->setN(3,this->getNz());
        FP->setN(4,this->getNt());
        FP->setD(1,this->getDx());
        FP->setD(3,this->getDz());
        FP->setD(4,this->getDt());
        FP->setO(1,this->getOx());
        FP->setO(3,this->getOz());
        FP->setO(4,this->getOt());
        FP->setData_format(sizeof(T));
        FP->writeHeader();
    }else{
        std::cerr << "SnapAcoustic2D::PutPfile: No filename set.\n";
        exit(1);
    }
}

template<typename T>
void SnapAcoustic2D<T>::putAxfile(std::string filename) {
    if(!filename.empty()){
        Axfile = filename;
        FAx = std::make_shared<File>();
        FAx->output(Axfile);
        Ax = true;
        FAx->setN(1,this->getNx());
        FAx->setN(3,this->getNz());
        FAx->setN(4,this->getNt());
        FAx->setD(1,this->getDx());
        FAx->setD(3,this->getDz());
        FAx->setD(4,this->getDt());
        FAx->setO(1,this->getOx());
        FAx->setO(3,this->getOz());
        FAx->setO(4,this->getOt());
        FAx->setData_format(sizeof(T));
        FAx->writeHeader();
    }else{
        std::cerr << "SnapAcoustic2D::PutAxfile: No filename set.\n";
        exit(1);
    }
}

template<typename T>
void SnapAcoustic2D<T>::putAzfile(std::string filename) {
    if(!filename.empty()){
        Azfile = filename;
        FAz = std::make_shared<File>();
        FAz->output(Azfile);
        Az = true;
        FAz->setN(1,this->getNx());
        FAz->setN(3,this->getNz());
        FAz->setN(4,this->getNt());
        FAz->setD(1,this->getDx());
        FAz->setD(3,this->getDz());
        FAz->setD(4,this->getDt());
        FAz->setO(1,this->getOx());
        FAz->setO(3,this->getOz());
        FAz->setO(4,this->getOt());
        FAz->setData_format(sizeof(T));
        FAz->writeHeader();
    }else{
        std::cerr << "SnapAcoustic2D::PutAzfile: No filename set.\n";
        exit(1);
    }
}

template<typename T>
SnapAcoustic2D<T>::~SnapAcoustic2D() {
	/* Do nothing */
}

// =============== 3D ACOUSTIC WAVES CLASS =============== //
template<typename T>
SnapAcoustic3D<T>::SnapAcoustic3D(std::shared_ptr<ModelAcoustic3D<T>> model, const int snap_inc){
    // Nothing here
}

template<typename T>
SnapAcoustic3D<T>::~SnapAcoustic3D() {
}

// =============== 2D ELASTIC WAVES CLASS =============== //
/** The 2D Elastic WAVES model class
 *
 */

template<typename T>
SnapElastic2D<T>::SnapElastic2D(std::shared_ptr<ModelElastic2D<T>> model, const int snap_inc)	///< Constructor
{
    // Do nothing
}

template<typename T>
SnapElastic2D<T>::~SnapElastic2D() {
}

// =============== 3D ELASTIC WAVES CLASS =============== //
/** The 3D Elastic WAVES model class
 *
 */

template<typename T>
SnapElastic3D<T>::SnapElastic3D(std::shared_ptr<ModelElastic3D<T>> model, const int snap_inc){
    // Nothing here
}

template<typename T>
SnapElastic3D<T>::~SnapElastic3D() {
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class SnapAcoustic2D<float>;
template class SnapAcoustic3D<float>;
template class SnapElastic2D<float>;
template class SnapElastic3D<float>;

}


