#include "data.h"

namespace rockseis {
// constructor
template<typename T>
Data<T>::Data(const std::string file, const int _ntrace, const int _nt, const T _dt, const T _ot)
{
    filename=file;
    nt = _nt;
    dt = _dt;
    ot = _ot;
    ntrace = _ntrace;
}

template<typename T>
Data<T>::Data(const int _ntrace, const int _nt, const T _dt)
{
    ntrace = _ntrace;
    nt = _nt;
    dt = _dt;
    ot = 0.0;
}

template<typename T>
Data<T>::Data(const int _ntrace, const int _nt, const T _dt, const T _ot)
{
    ntrace = _ntrace;
    nt = _nt;
    dt = _dt;
    ot = _ot;
}

// destructor
template<typename T>
Data<T>::~Data(){
    /* Do nothing*/
}

/// =============== 2D DATA CLASS =============== //

// constructor
template<typename T>
Data2D<T>::Data2D(const int _ntrace, const int _nt, const T _dt): Data<T>(_ntrace, _nt, _dt)
{
    // Setting variables
    this->setNtrace(_ntrace);
    this->setNt(_nt);
    this->setDt(_dt);

    // Create a 2D data geometry
    geometry = std::make_shared<Geometry2D<T>>(_ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(_ntrace*_nt, sizeof(T));
}

template<typename T>
bool Data2D<T>::readfloatData(std::shared_ptr<File> Fin)
{
    size_t nt = Fin->getN(1);
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    bool status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == 4) )
    {
        // Create geometry
        int i;
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->floatread((float *) &scoords[i].x, 1);
            Fin->floatread((float *) &scoords[i].y, 1);
            Fin->floatread((float *) &gcoords[i].x, 1);
            Fin->floatread((float *) &gcoords[i].y, 1);
            // Read trace data
            Fin->floatread((float *) &data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::readfloatCoords(std::shared_ptr<File> Fin)
{
    bool status;
    size_t nt = Fin->getN(1);
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == 4) )
    {
        // Create geometry
        int i;
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Jump over trace data;
            Fin->seekg(Fin->getStartofdata() + i*(nt+4)*sizeof(float));
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->floatread((float *) &scoords[i].x, 1);
            Fin->floatread((float *) &scoords[i].y, 1);
            Fin->floatread((float *) &gcoords[i].x, 1);
            Fin->floatread((float *) &gcoords[i].y, 1);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::writefloatData(std::shared_ptr<File> Fout)
{
    bool status;
    size_t nt = this->getNt();
    double dt = this->getDt();
    int ntrace = this->getNtrace();

    if ( Fout->is_open() )
    {
        // Create geometry
        Fout->setN(1,nt);
        Fout->setD(1,dt);
        Fout->setO(1,0.0);
        Fout->setN(2,ntrace);
        Fout->setNheader(4);
        Fout->setData_format(4);
        Fout->setHeader_format(4);
        Fout->setType(DATA2D);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        int i;
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Write coordinates to data trace (First source x and y and then receiver x and y
            Fout->floatwrite((float *) &scoords[i].x, 1);
            Fout->floatwrite((float *) &scoords[i].y, 1);
            Fout->floatwrite((float *) &gcoords[i].x, 1);
            Fout->floatwrite((float *) &gcoords[i].y, 1);
            // Write trace data
            Fout->floatwrite((float *) &data[i*nt], nt);
        }	

        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

// destructor
template<typename T>
Data2D<T>::~Data2D() {
    // Freeing all variables
    free(data);
}

/// =============== 3D DATA CLASS =============== //

// constructor
template<typename T>
Data3D<T>::Data3D(const int _ntrace, const int _nt, const T _dt): Data<T>(_ntrace, _nt, _dt)
{
    // Setting variables
    this->setNtrace(_ntrace);
    this->setNt(_nt);
    this->setDt(_dt);

    // Create a 3D data geometry
    geometry =std::make_shared<Geometry3D<T>>(_ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(_ntrace*_nt, sizeof(T));
}

template<typename T>
bool Data3D<T>::readfloatData(std::shared_ptr<File> Fin)
{
    bool status;
    int nt = Fin->getN(1);
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == 6) )
    {
        // Create geometry
        int i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            //read coordinates from data trace (First source x and y and then receiver x and y
            Fin->floatread((float *) &scoords[i].x, 1);
            Fin->floatread((float *) &scoords[i].y, 1);
            Fin->floatread((float *) &scoords[i].z, 1);
            Fin->floatread((float *) &gcoords[i].x, 1);
            Fin->floatread((float *) &gcoords[i].y, 1);
            Fin->floatread((float *) &gcoords[i].z, 1);
            // Write trace data
            Fin->floatread((float *) &data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data3D<T>::readfloatCoords(std::shared_ptr<File> Fin)
{
    bool status;
    size_t nt = Fin->getN(1);
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == 4) )
    {
        // Create geometry
        int i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Jump over trace data;
            Fin->seekg(Fin->getStartofdata() + i*(nt+6)*sizeof(float));
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->floatread((float *) &scoords[i].x, 1);
            Fin->floatread((float *) &scoords[i].y, 1);
            Fin->floatread((float *) &scoords[i].z, 1);
            Fin->floatread((float *) &gcoords[i].x, 1);
            Fin->floatread((float *) &gcoords[i].y, 1);
            Fin->floatread((float *) &gcoords[i].z, 1);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data3D<T>::writefloatData(std::shared_ptr<File> Fout)
{
    bool status;
    int nt = this->getNt();
    double dt = this->getDt();
    int ntrace = this->getNtrace();

    if ( Fout->is_open() )
    {
        // Create geometry
        Fout->setN(1,nt);
        Fout->setD(1,dt);
        Fout->setO(1,0.0);
        Fout->setN(2,ntrace);
        Fout->setNheader(6);
        Fout->setData_format(4);
        Fout->setHeader_format(4);
        Fout->setType(DATA3D);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        int i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            //Write coordinates to data trace (First source x and y and then receiver x and y
            Fout->floatwrite((float *) &scoords[i].x, 1);
            Fout->floatwrite((float *) &scoords[i].y, 1);
            Fout->floatwrite((float *) &scoords[i].z, 1);
            Fout->floatwrite((float *) &gcoords[i].x, 1);
            Fout->floatwrite((float *) &gcoords[i].y, 1);
            Fout->floatwrite((float *) &gcoords[i].z, 1);
            // Write trace data
            Fout->floatwrite((float *) &data[i*nt], nt);
        }	

        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

// destructor
template<typename T>
Data3D<T>::~Data3D() {
    // Freeing all variables
    free(data);
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Data<float>;
template class Data2D<float>;
template class Data3D<float>;

}


