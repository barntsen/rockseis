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
    verbose = false;
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
        if(n[i] != 0) std::cerr << "n[" << i << "] = " << n[i] << "    d[" << i << "] = " << d[i] << "    o[" << i << "] = " << o[i] << "\n"; 
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
size_t Geometry<T>::getNtot(){
    int i;
    size_t Ntot = 1;
    for(i = 0; i < MAXDIMS; i++)
    {    
	    if(n[i] > 0) 
	    {
            Ntot *= n[i];
	    }
    }
    return Ntot;
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
    scoords = (Point2D<T> *) calloc(ntrace,sizeof(Point2D<T>));
    gcoords = (Point2D<T> *) calloc(ntrace,sizeof(Point2D<T>));
    smap = (Point2D<int> *) calloc(ntrace,sizeof(Point2D<int>));
    gmap = (Point2D<int> *) calloc(ntrace,sizeof(Point2D<int>));
    sshift = (Point2D<T> *) calloc(ntrace,sizeof(Point2D<T>));
    gshift = (Point2D<T> *) calloc(ntrace,sizeof(Point2D<T>));
    for (size_t i=0; i< ntrace; i++) 
    {
        scoords[i].x = 0;
        scoords[i].y = 0;
        gcoords[i].x = 0;
        gcoords[i].y = 0;
        smap[i].x = -1;
        smap[i].y = -1;
        gmap[i].x = -1;
        gmap[i].y = -1;
    }
}

// destructor
template<typename T>
Geometry2D<T>::~Geometry2D() {
    // Freeing all variables
    free(scoords);
    free(gcoords);
    free(smap);
    free(gmap);
    free(sshift);
    free(gshift);
}

// create map
template<typename T>
void Geometry2D<T>::makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padhx, int padhy, T shiftx, T shifty) {
	size_t n = this->getN(1);  //Get number of traces 
	// Get regular model parameters
	int nx = _geom->getN(1);
	int ny = _geom->getN(3);
	double dx = (double) _geom->getD(1);
    double dy = (double) _geom->getD(3);
    double ox = (double) _geom->getO(1);
    double oy = (double) _geom->getO(3);

    int pos;
    double res;
    if(map == SMAP)
    {
        // Compute index smap
        for (size_t i = 0; i < n ; i++){
            res = ((double) scoords[i].x - ox + shiftx*dx)/dx;
            pos = (int) std::floor(res);
            //Avoiding roundoff errors
            if(CLOSE(res-pos,1.0,dx*1e-6))
               pos += 1;
            if(pos >=padlx && pos < nx-padhx)
            {
                smap[i].x  = pos; // index is within bounds
                sshift[i].x = (((double) scoords[i].x - ox + shiftx*dx)/dx) - pos;
            }else
            {
                smap[i].x  = -1;  // index is off bounds
            }
            res = ((double) scoords[i].y - oy + shifty*dy)/dy;
            pos = (int) std::floor(res);
            //Avoiding roundoff errors
            if(CLOSE(res-pos,1.0,dy*1e-6))
               pos += 1;
            if(pos >=padly && pos < ny-padhy)
            {
                smap[i].y  = pos; // index is within bounds
                sshift[i].y = (((double) scoords[i].y - oy + shifty*dy)/dy) - pos;
            }else
            {
                smap[i].y  = -1;  // index is off bounds
            }
        }
        // Issue a warning if all coordinates are out of bounds
        bool s_inbound = false;
        for (size_t i =0; i < n; i++) 
        {
            if ((smap[i].x >= 0)  && (smap[i].y >= 0)) s_inbound = true;
        }
        if (!s_inbound && this->getVerbose()) rs_warning("All source positions out of bounds, modelling might produce only zero output.");

    }else{
        // Compute index gmap
        for (size_t i = 0; i < n ; i++){
            res = ((double) gcoords[i].x - ox + shiftx*dx)/dx;
            pos = (int) std::floor(res);
            //Avoiding roundoff errors
            if(CLOSE(res-pos,1.0,dx*1e-6))
               pos += 1;
            if(pos >=padlx && pos < nx-padhx)
            {
                gmap[i].x  = pos; // index is within bounds
                gshift[i].x = (((double) gcoords[i].x - ox + shiftx*dx)/dx) - pos;
            }else
            {
                gmap[i].x  = -1;  // index is off bounds
            }
            res = ((double) gcoords[i].y - oy + shifty*dy)/dy;
            pos = (int) std::floor(res);
            //Avoiding roundoff errors
            if(CLOSE(res-pos,1.0,dy*1e-6))
               pos += 1;
            if(pos >=padly && pos < ny-padhy)
            {
                gmap[i].y  = pos; // index is within bounds
                gshift[i].y = (((double) gcoords[i].y - oy + shifty*dy)/dy) - pos;
            }else
            {
                gmap[i].y  = -1;  // index is off bounds
            }
        }
        // Issue a warning if all coordinates are out of bounds
        bool g_inbound = false;
        for (size_t i =0; i < n; i++) 
        {
            if ((gmap[i].x >= 0)  && (gmap[i].y >= 0)) g_inbound = true;
        }
        if (!g_inbound && this->getVerbose()) rs_warning("All receiver positions are out of bounds, modelling might produce only zero output.");
    }

}

// copy map
template<typename T>
void Geometry2D<T>::copyGmap2Smap() {
    size_t n = this->getN(1);  //Get number of traces 
    for (size_t i = 0; i < n; i++) 
    {
        smap[i].x = gmap[i].x;
        smap[i].y = gmap[i].y;
    }
}

// copy map
template<typename T>
void Geometry2D<T>::copySmap2Gmap() {
    size_t n = this->getN(1);  //Get number of traces 
    for (size_t i = 0; i < n; i++) 
    {
        gmap[i].x = smap[i].x;
        gmap[i].y = smap[i].y;
    }
}

/// =============== 3D DATA GEOMETRY CLASS =============== //

// constructor
template<typename T>
Geometry3D<T>::Geometry3D(int ntrace): Geometry<T>()
{
    this->setN(1, ntrace);
    scoords = (Point3D<T> *) calloc(ntrace,sizeof(Point3D<T>));
    gcoords = (Point3D<T> *) calloc(ntrace,sizeof(Point3D<T>));
    smap = (Point3D<int> *) calloc(ntrace,sizeof(Point3D<int>));
    gmap = (Point3D<int> *) calloc(ntrace,sizeof(Point3D<int>));
    sshift = (Point3D<T> *) calloc(ntrace,sizeof(Point3D<T>));
    gshift = (Point3D<T> *) calloc(ntrace,sizeof(Point3D<T>));
    for (size_t i=0; i< ntrace; i++) 
    {
        scoords[i].x = 0;
        scoords[i].y = 0;
        scoords[i].z = 0;
        gcoords[i].x = 0;
        gcoords[i].y = 0;
        gcoords[i].z = 0;
        smap[i].x = -1;
        smap[i].y = -1;
        smap[i].z = -1;
        gmap[i].x = -1;
        gmap[i].y = -1;
        gmap[i].z = -1;
    }

}

template<typename T>
Geometry3D<T>::~Geometry3D() {
    // Freeing all variables
    free(scoords);
    free(gcoords);
    free(smap);
    free(gmap);
    free(sshift);
    free(gshift);
}

// create map
template<typename T>
void Geometry3D<T>::makeMap(std::shared_ptr<Geometry<T>> _geom, bool map, int padlx, int padly, int padlz, int padhx, int padhy, int padhz, T shiftx, T shifty, T shiftz) {
    size_t n = this->getN(1);  //Get number of traces 
    // Get regular model parameters
    int nx = _geom->getN(1);
    int ny = _geom->getN(2);
    int nz = _geom->getN(3);
    double dx = (double) _geom->getD(1);
    double dy = (double) _geom->getD(2);
    double dz = (double) _geom->getD(3);
    double ox = (double) _geom->getO(1);
    double oy = (double) _geom->getO(2);
    double oz = (double) _geom->getO(3);
    // Compute index map
    int pos;
    if(map == SMAP){
        for (size_t i = 0; i < n ; i++){
            pos = (int) std::floor((scoords[i].x - ox + shiftx*dx)/dx);
            if(pos >=padlx && pos < nx-padhx)
            {
                smap[i].x  = pos; // index is within bounds
                sshift[i].x = ((scoords[i].x - ox + shiftx*dx)/dx) - pos;
            }else
            {
                smap[i].x  = -1;  // index is off bounds
            }
            pos = (int) std::floor((scoords[i].y - oy + shifty*dy)/dy);
            if(pos >=padly && pos < ny-padhy)
            {
                smap[i].y  = pos; // index is within bounds
                sshift[i].y = ((scoords[i].y - oy + shifty*dy)/dy) - pos;
            }else
            {
                smap[i].y  = -1;  // index is off bounds
            }
            pos = (int) std::floor((scoords[i].z - oz + shiftz*dz)/dz);
            if(pos >=padlz && pos < nz-padhz)
            {
                smap[i].z  = pos; // index is within bounds
                sshift[i].z = ((scoords[i].z - oz + shiftz*dz)/dz) - pos;
            }else
            {
                smap[i].z  = -1;  // index is off bounds
            }
        }

        // Issue a warning if all coordinates are out of bounds
        bool s_inbound = false;

        for (size_t i = 0; i < n; i++) 
        {
            if ((smap[i].x >= 0)  && (smap[i].y >= 0)  && (smap[i].z >= 0)) s_inbound = true;
        }
        if (!s_inbound && this->getVerbose()) rs_warning("All source positions out of bounds, modelling might produce only zero output.");
    }else{
        for (size_t i = 0; i < n ; i++){
            pos = (int) std::floor((gcoords[i].x - ox + shiftx*dx)/dx);
            if(pos >=padlx && pos < nx-padhx)
            {
                gmap[i].x  = pos; // index is within bounds
                gshift[i].x = ((gcoords[i].x - ox + shiftx*dx)/dx) - pos;
            }else
            {
                gmap[i].x  = -1;  // index is off bounds
            }
            pos = (int) std::floor((gcoords[i].y - oy + shifty*dy)/dy);
            if(pos >=padly && pos < ny-padhy)
            {
                gmap[i].y  = pos; // index is within bounds
                gshift[i].y = ((gcoords[i].y - oy + shifty*dy)/dy) - pos;
            }else
            {
                gmap[i].y  = -1;  // index is off bounds
            }
            pos = (int) std::floor((gcoords[i].z - oz + shiftz*dz)/dz);
            if(pos >=padlz && pos < nz-padhz)
            {
                gmap[i].z  = pos; // index is within bounds
                gshift[i].z = ((gcoords[i].z - oz + shiftz*dz)/dz) - pos;
            }else
            {
                gmap[i].z  = -1;  // index is off bounds
            }
        }

        // Issue a warning if all coordinates are out of bounds
        bool g_inbound = false;

        for (size_t i = 0; i < n; i++) 
        {
            if ((gmap[i].x >= 0)  && (gmap[i].y >= 0)  && (gmap[i].z >= 0)) g_inbound = true;
        }
        if (!g_inbound && this->getVerbose()) rs_warning("All receiver positions are out of bounds, modelling might produce only zero output.");
    }
}


// copy map
template<typename T>
void Geometry3D<T>::copyGmap2Smap() {
    size_t n = this->getN(1);  //Get number of traces 
    for (size_t i = 0; i < n; i++) 
    {
        smap[i].x = gmap[i].x;
        smap[i].y = gmap[i].y;
        smap[i].z = gmap[i].z;
    }
}

// copy map
template<typename T>
void Geometry3D<T>::copySmap2Gmap() {
    size_t n = this->getN(1);  //Get number of traces 
    for (size_t i = 0; i < n; i++) 
    {
        gmap[i].x = smap[i].x;
        gmap[i].y = smap[i].y;
        gmap[i].z = smap[i].z;
    }
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Geometry<float>;
template class Geometry<double>;
template class Geometry2D<float>;
template class Geometry2D<double>;
template class Geometry3D<float>;
template class Geometry3D<double>;
}


