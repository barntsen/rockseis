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
SnapAcoustic2D<T>::SnapAcoustic2D(std::shared_ptr<WavesAcoustic2D<T>> waves, const int snap_inc): Snap<T>(){
    int _nx, _ny, _nz, _nt, _nt_mod;
    T _dx, _dy, _dz, _dt, _dt_mod; 
    T _ox, _oy, _oz, _ot, _ot_mod; 
    int _lpml;
    int end_diff;

    /* Get necessary parameters from model class */
    _nx = waves->getNx();
    _ny = waves->getNy();
    _nz = waves->getNz();
    _nt_mod = waves->getNt();
    _dx = waves->getDx();
    _dy = waves->getDy();
    _dz = waves->getDz();
    _dt_mod = waves->getDt();
    _ox = waves->getOx();
    _oy = waves->getOy();
    _oz = waves->getOz();
    _ot_mod = waves->getOt();
    _lpml = waves->getLpml();

    if(snap_inc < 1) {
	    /* Error */ 
	    std::cerr << "SnapAcoustic2D::snap_inc must be an integer larger than 0\n";
	    exit(1);
    }
    if(snap_inc > _nt_mod) {
	    /* Error */ 
	    std::cerr << "SnapAcoustic2D::snap_inc is larger than nt. \n";
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
    this->setDx(_dx);
    this->setDy(_dy);
    this->setDz(_dz);
    this->setDt(_dt);
    this->setOx(_ox);
    this->setOy(_oy);
    this->setOz(_oz);
    this->setOt(_ot);
    this->setLpml(_lpml);
}

template<typename T>
SnapAcoustic2D<T>::~SnapAcoustic2D() {
	/* Do nothing */
}

// =============== 3D ACOUSTIC WAVES CLASS =============== //
template<typename T>
SnapAcoustic3D<T>::SnapAcoustic3D(std::shared_ptr<WavesAcoustic3D<T>> waves, const int snap_inc){
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
SnapElastic2D<T>::SnapElastic2D(std::shared_ptr<WavesElastic2D<T>> waves, const int snap_inc)	///< Constructor
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
SnapElastic3D<T>::SnapElastic3D(std::shared_ptr<WavesElastic3D<T>> waves, const int snap_inc){
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


