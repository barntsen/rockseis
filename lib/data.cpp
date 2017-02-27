#include "data.h"

namespace rockseis {
// constructor
template<typename T>
Data<T>::Data(const std::string file)
{
    datafile=file;
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
Data2D<T>::Data2D(std::string datafile): Data<T>(datafile)
{
    bool status;

    size_t ntrace;
    size_t nt;
    T dt;
    T ot;
    //Opeing file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data2D::Error reading from input file. \n";
	    exit(1);
    }
    if(Fin->getData_format() != sizeof(T))
    {
        std::cerr << "Data2D::Numerical precision in " << datafile << " mismatch with data class contructor.\n";
        exit(1);
    }
    if(Fin->getNheader() != 4)
    {
        std::cerr << "Data2D:: " << datafile << " is not a 2d data file.\n";
        exit(1);
    }

    ntrace = Fin->getN(2);
    nt = Fin->getN(1);
    dt = Fin->getD(1);
    ot = Fin->getO(1);

    // Setting variables
    this->setNtrace(ntrace);
    this->setNt(nt);
    this->setDt(dt);
    this->setOt(ot);

    // Create a 2D data geometry
    geometry = std::make_shared<Geometry2D<T>>(ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(ntrace*nt, sizeof(T));
}

template<typename T>
bool Data2D<T>::readData()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    std::cerr << "Data2D::readData no file assigned. \n";
	    exit(1);
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data2D::Error reading from " << datafile <<". \n";
	    exit(1);
    }
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
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->read(&scoords[i].x, 1);
            Fin->read(&scoords[i].y, 1);
            Fin->read(&gcoords[i].x, 1);
            Fin->read(&gcoords[i].y, 1);
            // Read trace data
            Fin->read(&data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::readCoords()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    std::cerr << "Data2D::readCoords: No file assigned. \n";
	    exit(1);
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data2D::Error reading from " << datafile <<". \n";
	    exit(1);
    }

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
            Fin->read((float *) &scoords[i].x, 1);
            Fin->read((float *) &scoords[i].y, 1);
            Fin->read((float *) &gcoords[i].x, 1);
            Fin->read((float *) &gcoords[i].y, 1);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::write()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    std::cerr << "Data2D::writeData: No file assigned. \n";
	    exit(1);
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(datafile.c_str());
    size_t nt = this->getNt();
    double dt = this->getDt();
    double ot = this->getOt();
    int ntrace = this->getNtrace();

    if ( Fout->is_open() )
    {
        // Create geometry
        Fout->setN(1,nt);
        Fout->setD(1,dt);
        Fout->setO(1,ot);
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
            Fout->write((float *) &scoords[i].x, 1);
            Fout->write((float *) &scoords[i].y, 1);
            Fout->write((float *) &gcoords[i].x, 1);
            Fout->write((float *) &gcoords[i].y, 1);
            // Write trace data
            Fout->write((float *) &data[i*nt], nt);
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
Data3D<T>::Data3D(std::string datafile): Data<T>(datafile)
{
    bool status;
    this->setFile(datafile);

    size_t ntrace;
    size_t nt;
    T dt;
    T ot;
    //Opeing file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data3D::Error reading from input file. \n";
	    exit(1);
    }
    if(Fin->getData_format() != sizeof(T))
    {
        std::cerr << "Data3D::Numerical precision in " << datafile << " mismatch with data class contructor.\n";
        exit(1);
    }
    if(Fin->getNheader() != 6)
    {
        std::cerr << "Data3D:: " << datafile << " is not a 3d data file.\n";
        exit(1);
    }

    ntrace = Fin->getN(2);
    nt = Fin->getN(1);
    dt = Fin->getD(1);
    ot = Fin->getO(1);

    // Setting variables
    this->setNtrace(ntrace);
    this->setNt(nt);
    this->setDt(dt);
    this->setOt(ot);

    // Create a 3D data geometry
    geometry = std::make_shared<Geometry3D<T>>(ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(ntrace*nt, sizeof(T));
}

template<typename T>
bool Data3D<T>::readData()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    std::cerr << "Data3D::readData no file assigned. \n";
	    exit(1);
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data3D::Error reading from " << datafile <<". \n";
	    exit(1);
    }

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
            Fin->read(&scoords[i].x, 1);
            Fin->read(&scoords[i].y, 1);
            Fin->read(&scoords[i].z, 1);
            Fin->read(&gcoords[i].x, 1);
            Fin->read(&gcoords[i].y, 1);
            Fin->read(&gcoords[i].z, 1);
            // Write trace data
            Fin->read(&data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data3D<T>::readCoords()
{
    bool status;
    std::string datafile = this->getFile();
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    std::cerr << "Data3D::readCoords : Error reading from file. \n";
	    exit(1);
    }
    size_t nt = Fin->getN(1);
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();
    int data_format = Fin->getData_format();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == 6) && (data_format == sizeof(T)))
    {
        // Create geometry
        int i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Jump over trace data;
            Fin->seekg(Fin->getStartofdata() + i*(nt+6)*sizeof(T));
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->read(&scoords[i].x, 1);
            Fin->read(&scoords[i].y, 1);
            Fin->read(&scoords[i].z, 1);
            Fin->read(&gcoords[i].x, 1);
            Fin->read(&gcoords[i].y, 1);
            Fin->read(&gcoords[i].z, 1);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data3D<T>::write()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    std::cerr << "Data2D::writeData: No file assigned. \n";
	    exit(1);
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    Fout->output(datafile.c_str());

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
            Fout->write((float *) &scoords[i].x, 1);
            Fout->write((float *) &scoords[i].y, 1);
            Fout->write((float *) &scoords[i].z, 1);
            Fout->write((float *) &gcoords[i].x, 1);
            Fout->write((float *) &gcoords[i].y, 1);
            Fout->write((float *) &gcoords[i].z, 1);
            // Write trace data
            Fout->write((float *) &data[i*nt], nt);
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


