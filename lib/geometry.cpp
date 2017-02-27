#include "geometry.h"

namespace rockseis {

/// =============== ABSTRACT GEOMETRY CLASS =============== //
// constructor
template<typename T>
Geometry<T>::Geometry()
{
    int i;
    for(i = 0; i < MAXDIMS; i++)
    {
        n[i] = 0;
        d[i] = 0.0;
        o[i] = 0.0;
    } 
}

// Clear geometry
template<typename T>
void Geometry<T>::clear(){
    int i;
    for(i = 0; i < MAXDIMS; i++)
    {    
        n[i] = 0;
        d[i] = 0.0;
        o[i] = 0.0;
    }
}

// Print geometry info
template<typename T>
void Geometry<T>::print(){
    std::cerr << "==================================================\n";
    std::cerr << "==================== Geometry ====================\n";
    int i;
    for(i = 0; i < MAXDIMS; i++)
    {    
        std::cerr << "n[" << i << "] = " << n[i] << "    d[" << i << "] = " << d[i] << "    o[" << i << "] = " << o[i] << "\n"; 
    }
    std::cerr << "==================================================\n";
}

template<typename T>
bool Geometry<T>::compare(std::shared_ptr<Geometry<T>> other){
    int i;
    for(i = 0; i < MAXDIMS; i++)
    {    
	    if(n[i] > 0) 
	    {
		    if(n[i] != other->getN(i+1)) return 1; 
		    if(d[i] != other->getD(i+1)) return 1; 
		    if(o[i] != other->getO(i+1)) return 1; 
	    }
    }
    return 0;
}



template<typename T>
void Geometry<T>::setN(int dim, size_t val) 
{
    if( (dim <= MAXDIMS) && (dim > 0) ) {
        n[dim-1] = val;  
    }
    if(val > 0){
        // If subdimensions are 0, set them to be 1 (Non-singletons)
        for (int i = 0; i < dim-1; i++)
        {
            if(n[i] == 0) n[i] = 1;
        }
    }
}


// destructor
template<typename T>
Geometry<T>::~Geometry(){
    /* Do nothing*/
}


/// =============== 2D DATA GEOMETRY CLASS =============== //

// constructor
template<typename T>
Geometry2D<T>::Geometry2D(size_t ntrace): Geometry<T>()
{
    this->setN(1, ntrace);
    scoords = (rockseis::Point2D<T> *) calloc(ntrace,sizeof(rockseis::Point2D<T>));
    gcoords = (rockseis::Point2D<T> *) calloc(ntrace,sizeof(rockseis::Point2D<T>));
    smap = (rockseis::Point2D<int> *) calloc(ntrace,sizeof(rockseis::Point2D<int>));
    gmap = (rockseis::Point2D<int> *) calloc(ntrace,sizeof(rockseis::Point2D<int>));
}

// destructor
template<typename T>
Geometry2D<T>::~Geometry2D() {
    // Freeing all variables
    free(scoords);
    free(gcoords);
    free(smap);
    free(gmap);
}


// create map
template<typename T>
void Geometry2D<T>::makeMap(std::shared_ptr<Geometry<T>> _geom) {
	size_t n = this->getN(1);  //Get number of traces 
	// Get regular model parameters
	int nx = _geom->getN(1);
	int ny = _geom->getN(3);
	T dx = _geom->getD(1);
	T dy = _geom->getD(3);
	T ox = _geom->getO(1);
	T oy = _geom->getO(3);
	// Compute index smap
	int pos;
	for (size_t i = 0; i < n ; i++){
		pos =  (int) round(scoords[i].x - ox)/dx;
		if(pos >=0 && pos < nx)
		{
			smap[i].x  = pos; // index is within bounds
		}else
		{
			smap[i].x  = -1;  // index is off bounds
		}
		pos =  (int) round(scoords[i].y - oy)/dy;
		if(pos >=0 && pos < ny)
		{
			smap[i].y  = pos; // indey is within bounds
		}else
		{
			smap[i].y  = -1;  // indey is off bounds
		}
	}
	// Compute index gmap
	for (size_t i = 0; i < n ; i++){
		pos =  (int) round(gcoords[i].x - ox)/dx;
		if(pos >=0 && pos < nx)
		{
			gmap[i].x  = pos; // index is within bounds
		}else
		{
			gmap[i].x  = -1;  // index is off bounds
		}
		pos =  (int) round(gcoords[i].y - oy)/dy;
		if(pos >=0 && pos < ny)
		{
			gmap[i].y  = pos; // indey is within bounds
		}else
		{
			gmap[i].y  = -1;  // indey is off bounds
		}
	}
}


/// =============== 3D DATA GEOMETRY CLASS =============== //

// constructor
template<typename T>
Geometry3D<T>::Geometry3D(int ntrace): Geometry<T>()
{
    this->setN(1, ntrace);
    scoords = (rockseis::Point3D<T> *) calloc(ntrace,sizeof(rockseis::Point3D<T>));
    gcoords = (rockseis::Point3D<T> *) calloc(ntrace,sizeof(rockseis::Point3D<T>));
    smap = (rockseis::Point3D<int> *) calloc(ntrace,sizeof(rockseis::Point3D<int>));
    gmap = (rockseis::Point3D<int> *) calloc(ntrace,sizeof(rockseis::Point3D<int>));
}

template<typename T>
Geometry3D<T>::~Geometry3D() {
    // Freeing all variables
    free(scoords);
    free(gcoords);
    free(smap);
    free(gmap);
}

// create map
template<typename T>
void Geometry3D<T>::makeMap(std::shared_ptr<Geometry<T>> _geom) {
	size_t n = this->getN(1);  //Get number of traces 
	// Get regular model parameters
	int nx = _geom->getN(1);
	int ny = _geom->getN(2);
	int nz = _geom->getN(3);
	T dx = _geom->getD(1);
	T dy = _geom->getD(2);
	T dz = _geom->getD(3);
	T ox = _geom->getO(1);
	T oy = _geom->getO(2);
	T oz = _geom->getO(3);
	// Compute index map
	int pos;
	for (size_t i = 0; i < n ; i++){
		pos =  (int) round(scoords[i].x - ox)/dx;
		if(pos >=0 && pos < nx)
		{
			smap[i].x  = pos; // index is within bounds
		}else
		{
			smap[i].x  = -1;  // index is off bounds
		}
		pos =  (int) round(scoords[i].y - oy)/dy;
		if(pos >=0 && pos < ny)
		{
			smap[i].y  = pos; // indey is within bounds
		}else
		{
			smap[i].y  = -1;  // indey is off bounds
		}
		pos =  (int) round(scoords[i].z - oz)/dz;
		if(pos >=0 && pos < nz)
		{
			smap[i].z  = pos; // indez is within bounds
		}else
		{
			smap[i].z  = -1;  // indez is off bounds
		}
	}
	for (size_t i = 0; i < n ; i++){
		pos =  (int) round(gcoords[i].x - ox)/dx;
		if(pos >=0 && pos < nx)
		{
			gmap[i].x  = pos; // index is within bounds
		}else
		{
			gmap[i].x  = -1;  // index is off bounds
		}
		pos =  (int) round(gcoords[i].y - oy)/dy;
		if(pos >=0 && pos < ny)
		{
			gmap[i].y  = pos; // indey is within bounds
		}else
		{
			gmap[i].y  = -1;  // indey is off bounds
		}
		pos =  (int) round(gcoords[i].z - oz)/dz;
		if(pos >=0 && pos < nz)
		{
			gmap[i].z  = pos; // indez is within bounds
		}else
		{
			gmap[i].z  = -1;  // indez is off bounds
		}
	}

}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Geometry<float>;
template class Geometry<double>;
template class Geometry2D<float>;
template class Geometry3D<float>;
}


