#include "data.h"

namespace rockseis {
// constructor
template<typename T>
Data<T>::Data(const std::string file)
{
    datafile=file;
    field = PRESSURE;
    Fdata = std::make_shared<File>(); 
}

template<typename T>
Data<T>::Data(const int _ntrace, const int _nt, const T _dt)
{
    ntrace = _ntrace;
    nt = _nt;
    dt = _dt;
    ot = 0.0;
    field = PRESSURE;
    Fdata = std::make_shared<File>(); 
}

template<typename T>
Data<T>::Data(const int _ntrace, const int _nt, const T _dt, const T _ot)
{
    ntrace = _ntrace;
    nt = _nt;
    dt = _dt;
    ot = _ot;
    field = PRESSURE;
    Fdata = std::make_shared<File>(); 
}

template<typename T>
bool Data<T>::open(std::string flag)
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
        rs_error("Data::open: No file assigned.");
    }

    switch ((flag.c_str())[0])
    {
        case 'i':
            status = Fdata->input(datafile);
            if(status == FILE_ERR) rs_error("Error opening file");
            Fdata->seekg(Fdata->getStartofdata());
            break;
        case 'o':
            Fdata->output(datafile);
            Fdata->writeHeader();
            Fdata->seekp(Fdata->getStartofdata());
            status = FILE_OK;
            break;
        case 'a':
            status = Fdata->append(datafile);
            break;
        default:
            status = FILE_ERR;
            break;
    }
    return status;
}

template<typename T>
void Data<T>::close()
{
    Fdata->close();
}

// destructor
template<typename T>
Data<T>::~Data(){
    close();
}

/// =============== 2D DATA CLASS =============== //

// constructor
template<typename T>
Data2D<T>::Data2D(const int _ntrace, const int _nt, const T _dt, const T _ot): Data<T>(_ntrace, _nt, _dt, _ot)
{
    // Create a 2D data geometry
    geometry = std::make_shared<Geometry2D<T>>(_ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(_ntrace*_nt, sizeof(T));

    std::shared_ptr<File> Fdata = this->getFdata();
    Fdata->setN(1,_nt);
    Fdata->setD(1,_dt);
    Fdata->setO(1,_ot);
    Fdata->setN(2,0);
    Fdata->setD(2,1);
    Fdata->setO(2,0);
    Fdata->setNheader(NHEAD2D);
    Fdata->setData_format(sizeof(T));
    Fdata->setHeader_format(sizeof(T));
    Fdata->setType(rockseis::DATA2D);
}

template<typename T>
Data2D<T>::Data2D(std::string datafile): Data<T>(datafile)
{
    bool status;

    size_t ntrace;
    size_t nt;
    T dt;
    T ot;
    //Opening file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data2D::Error reading from ", datafile );
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Data2D::Numerical precision in ", datafile, " mismatch with data class contructor.");
    }
    if(Fin->getNheader() != NHEAD2D)
    {
        rs_error("Data2D:: " , datafile , " is not a 2d data file.");
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
Data2D<T>::Data2D(std::string datafile, const int _nt, const T _dt, const T _ot): Data<T>(datafile)
{
    bool status;
    size_t ntrace;
    //Opeing file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data2D::Error reading from " , datafile ,". ");
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Data2D::Numerical precision in " , datafile , " mismatch with data class contructor.");
    }
    if(Fin->getNheader() != NHEAD2D)
    {
        rs_error("Data2D:: " , datafile , " is not a 2d data file.");
    }

    ntrace = Fin->getN(2);

    // Setting variables
    this->setNtrace(ntrace);
    this->setNt(_nt);
    this->setDt(_dt);
    this->setOt(_ot);

    // Create a 2D data geometry
    geometry = std::make_shared<Geometry2D<T>>(ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(ntrace*_nt, sizeof(T));
}

template<typename T>
bool Data2D<T>::read()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    rs_error("Data2D::read: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data2D::read: Error reading from " , datafile ,". ");
    }
    size_t nt = Fin->getN(1);
    if(nt != this->getNt())
    {
	    rs_error("Data2D::read: Nt in " , datafile ," is different from that allocated during class construction. ");
	    exit(1);
    }
    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == NHEAD2D) )
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
	    rs_error("Data2D::readCoords: No file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data2D::readCoords: Error reading from " , datafile ,". ");
    }

    size_t nt = Fin->getN(1);
    size_t ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == NHEAD2D) )
    {
        // Create geometry
        size_t i;
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Jump over trace data;
            Fin->seekg(Fin->getStartofdata() + i*(nt+NHEAD2D)*sizeof(T));
            //Read coordinates from data trace (First source x and y and then receiver x and y
            Fin->read(&scoords[i].x, 1);
            Fin->read(&scoords[i].y, 1);
            Fin->read(&gcoords[i].x, 1);
            Fin->read(&gcoords[i].y, 1);
        }	
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::copyCoords(std::shared_ptr<Data2D<T>> Data)
{
    bool status;
    size_t ntr1 = this->getNtrace();
    size_t ntr2 = Data->getNtrace();

    if(ntr1 <= ntr2)
    {
    // Copy geometry
        size_t i;
        Point2D<T> *scoords1 = (this->getGeom())->getScoords();
        Point2D<T> *gcoords1 = (this->getGeom())->getGcoords();

        Point2D<T> *scoords2 = (Data->getGeom())->getScoords();
        Point2D<T> *gcoords2 = (Data->getGeom())->getGcoords();
        for (i=0; i < ntr1; i++)
        {
            scoords1[i].x = scoords2[i].x;
            scoords1[i].y = scoords2[i].y;
            gcoords1[i].x = gcoords2[i].x;
            gcoords1[i].y = gcoords2[i].y;
        }	
        status = FILE_OK;	
    }else{
        rs_error("Data2D<T>::copyCoords : Mismatch in number of traces of Data2D objects");
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
	    rs_error("Data2D::writeData: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    if(strcmp(datafile.c_str(), "stdout")) {
        Fout->output(datafile.c_str());
    }else{
        Fout->output();
    }

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
        Fout->setD(2,1.0);
        Fout->setO(2,0.0);
        Fout->setNheader(NHEAD2D);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(DATA2D);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        int i;
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Write coordinates to data trace (First source x and y and then receiver x and y
            Fout->write(&scoords[i].x, 1);
            Fout->write(&scoords[i].y, 1);
            Fout->write(&gcoords[i].x, 1);
            Fout->write(&gcoords[i].y, 1);
            // Write trace data
            Fout->write(&data[i*nt], nt);
        }	

        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data2D<T>::readTraces()
{
    bool status;
    int nt = this->getNt();
    int ntrace = this->getNtrace();
    std::shared_ptr<rockseis::File> Fin;
    Fin = this->getFdata();

    if ( Fin->is_open() )
    {
        // Update geometry
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (int i=0; i < ntrace; i++)
        {
            //Read coordinates to data trace (First source x and y and then receiver x and y
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

    if(Fin->getFail()){
        status = FILE_ERR;
    }else{
        status = FILE_OK;
    }
    return status;
}

template<typename T>
bool Data2D<T>::writeTraces()
{
    bool status;
    int nt = this->getNt();
    int ntrace = this->getNtrace();
    std::shared_ptr<rockseis::File> Fout;
    Fout = this->getFdata();

    if ( Fout->is_open() )
    {
        // Update geometry
        size_t n2 = Fout->getN(2);
        Fout->setN(2, n2 + ntrace);
        Point2D<T> *scoords = (this->getGeom())->getScoords();
        Point2D<T> *gcoords = (this->getGeom())->getGcoords();
        for (int i=0; i < ntrace; i++)
        {
            //Write coordinates to data trace (First source x and y and then receiver x and y
            Fout->write(&scoords[i].x, 1);
            Fout->write(&scoords[i].y, 1);
            Fout->write(&gcoords[i].x, 1);
            Fout->write(&gcoords[i].y, 1);
            // Write trace data
            Fout->write(&data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }
    if(Fout->getFail()){
        status = FILE_ERR;
    }else{
        status = FILE_OK;
    }
    return status;
}

template<typename T>
void Data2D<T>::createEmpty(size_t ntr)
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
        rs_error("Data2D::createEmpty: No file assigned. ");
    }

    status = this->open("o");
    if(status == FILE_ERR) rockseis::rs_error("Data2D::createEmpty: Error opening file for writting");

    for (size_t i=0; i < ntr; i++){
        this->writeTraces();
    }
    this->close();
}

template<typename T>
void Data2D<T>::putTrace(std::string filename, size_t number)
{
    size_t traceno;
    traceno = number;

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(filename);
    if(status == FILE_ERR) rs_error("Data2D::putTrace: Error opening output data file.");
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA2D) rs_error("Data2D::putTrace: Datafile must be of type Data2D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = 1; // Only write one trace!
    if(n1 != this->getNt()) rs_error("Data2D::putTrace: Number of samples in trace and datafile mismatch.");
    if(d1 != this->getDt()) rs_error("Data2D::putTrace: Sampling interval in trace and datafile mismatch.");
    if(o1 != this->getOt()) rs_error("Data2D::putTrace: Origin in trace and datafile mismatch.");

    if(traceno > (Fdata->getN(2)-1)) rs_error("Data2D::putTrace: Trying to put a trace with number that is larger than number of traces in file. ", std::to_string(traceno));

    //Write gather
    Point2D<T> *scoords = (this->getGeom())->getScoords();
    Point2D<T> *gcoords = (this->getGeom())->getGcoords();
    T *tracedata = this->getData();
    for (size_t j=0; j < n2; j++){
        Fdata->seekp(Fdata->getStartofdata() + traceno*(n1+NHEAD2D)*sizeof(T));
        Fdata->write(&scoords[j].x, 1);
        Fdata->write(&scoords[j].y, 1);
        Fdata->write(&gcoords[j].x, 1);
        Fdata->write(&gcoords[j].y, 1);
        Fdata->write(&tracedata[j*n1], n1);
    }

    if(Fdata->getFail()) rs_error("Data2D::PutTrace: Error writting gather to output file");
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
Data3D<T>::Data3D(const int _ntrace, const int _nt, const T _dt, const T _ot): Data<T>(_ntrace, _nt, _dt, _ot)
{
    // Create a 3D data geometry
    geometry = std::make_shared<Geometry3D<T>>(_ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(_ntrace*_nt, sizeof(T));

    // Create a geometry for the header
    std::shared_ptr<File> Fdata = this->getFdata();
    Fdata->setN(1,_nt);
    Fdata->setD(1,_dt);
    Fdata->setO(1,_ot);
    Fdata->setN(2,0);
    Fdata->setD(2,1.0);
    Fdata->setO(2,0.0);
    Fdata->setNheader(NHEAD3D);
    Fdata->setData_format(sizeof(T));
    Fdata->setHeader_format(sizeof(T));
    Fdata->setType(rockseis::DATA3D);
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

    //Opening file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data3D::Error reading from " , datafile ,". ");
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Data3D::Numerical precision in " , datafile , " mismatch with data class contructor.");
    }
    if(Fin->getNheader() != NHEAD3D)
    {
        rs_error("Data3D:: " , datafile , " is not a 3d data file.");
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
Data3D<T>::Data3D(std::string datafile, const int _nt, const T _dt, const T _ot): Data<T>(datafile)
{
    bool status;
    size_t ntrace;
    //Opeing file for reading
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data3D::Error reading from " , datafile ,". ");
    }
    if(Fin->getData_format() != sizeof(T))
    {
        rs_error("Data3D::Numerical precision in " , datafile , " mismatch with data class contructor.");
    }
    if(Fin->getNheader() != NHEAD3D)
    {
        rs_error("Data3D:: " , datafile , " is not a 3d data file.");
    }

    ntrace = Fin->getN(2);

    // Setting variables
    this->setNtrace(ntrace);
    this->setNt(_nt);
    this->setDt(_dt);
    this->setOt(_ot);

    // Create a 3D data geometry
    geometry = std::make_shared<Geometry3D<T>>(ntrace); 

    // Allocate the memory for the data
    data = (T *) calloc(ntrace*_nt, sizeof(T));
}

template<typename T>
bool Data3D<T>::read()
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
	    rs_error("Data3D::read no file assigned. ");
    }
    std::shared_ptr<rockseis::File> Fin (new rockseis::File());
    status = Fin->input(datafile.c_str());
    if(status == FILE_ERR){
	    rs_error("Data3D::Error reading from " , datafile ,". ");
    }

    int nt = Fin->getN(1);
    if(nt != this->getNt())
    {
	    rs_error("Data3D::read: Nt in " , datafile ," is different from that allocated during class construction. ");
    }

    int ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == NHEAD3D) )
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
            // Read trace data
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
	    rs_error("Data3D::readCoords : Error reading from " , datafile ,". ");
    }
    size_t nt = Fin->getN(1);
    size_t ntrace = Fin->getN(2);
    size_t Nheader = Fin->getNheader();
    int data_format = Fin->getData_format();

    status = FILE_OK;
    if ( Fin->is_open()  && (Nheader == NHEAD3D) && (data_format == sizeof(T)))
    {
        // Create geometry
        size_t i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            // Jump over trace data;
            Fin->seekg(Fin->getStartofdata() + i*(nt+NHEAD3D)*sizeof(T));
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
bool Data3D<T>::copyCoords(std::shared_ptr<Data3D<T>> Data)
{
    bool status;
    size_t ntr1 = this->getNtrace();
    size_t ntr2 = Data->getNtrace();

    if(ntr1 <= ntr2)
    {
    // Copy geometry
        size_t i;
        Point3D<T> *scoords1 = (this->getGeom())->getScoords();
        Point3D<T> *gcoords1 = (this->getGeom())->getGcoords();

        Point3D<T> *scoords2 = (Data->getGeom())->getScoords();
        Point3D<T> *gcoords2 = (Data->getGeom())->getGcoords();
        for (i=0; i < ntr1; i++)
        {
            scoords1[i].x = scoords2[i].x;
            scoords1[i].y = scoords2[i].y;
            scoords1[i].z = scoords2[i].z;
            gcoords1[i].x = gcoords2[i].x;
            gcoords1[i].y = gcoords2[i].y;
            gcoords1[i].z = gcoords2[i].z;
        }	
        status = FILE_OK;	
    }else{
        rs_error("Data3D<T>::copyCoords : Mismatch in number of traces of Data3D objects");
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
	    rs_error("Data3D::writeData: No file assigned. ");
    }

    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    if(strcmp(datafile.c_str(), "stdout")) {
        Fout->output(datafile.c_str());
    }else{
        Fout->output();
    }

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
        Fout->setD(2,1.0);
        Fout->setO(2,0.0);
        Fout->setNheader(NHEAD3D);
        Fout->setData_format(sizeof(T));
        Fout->setHeader_format(sizeof(T));
        Fout->setType(DATA3D);
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        int i;
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (i=0; i < ntrace; i++)
        {
            //Write coordinates to data trace (First source x and y and then receiver x and y
            Fout->write(&scoords[i].x, 1);
            Fout->write(&scoords[i].y, 1);
            Fout->write(&scoords[i].z, 1);
            Fout->write(&gcoords[i].x, 1);
            Fout->write(&gcoords[i].y, 1);
            Fout->write(&gcoords[i].z, 1);
            // Write trace data
            Fout->write(&data[i*nt], nt);
        }	
        status = FILE_OK;
    }else{
        status = FILE_ERR;	
    }
    return status;
}

template<typename T>
bool Data3D<T>::readTraces()
{
    bool status;
    int nt = this->getNt();
    int ntrace = this->getNtrace();
    std::shared_ptr<rockseis::File> Fin;
    Fin = this->getFdata();

    if ( Fin->is_open() )
    {
        // Update geometry
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (int i=0; i < ntrace; i++)
        {
            //Read coordinates to data trace (First source x and y and then receiver x and y
            Fin->read(&scoords[i].x, 1);
            Fin->read(&scoords[i].y, 1);
            Fin->read(&scoords[i].z, 1);
            Fin->read(&gcoords[i].x, 1);
            Fin->read(&gcoords[i].y, 1);
            Fin->read(&gcoords[i].z, 1);
            // Read trace data
            Fin->read(&data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }

    if(Fin->getFail()){
        status = FILE_ERR;
    }else{
        status = FILE_OK;
    }
    return status;
}

template<typename T>
bool Data3D<T>::writeTraces()
{
    bool status;
    int nt = this->getNt();
    int ntrace = this->getNtrace();
    std::shared_ptr<rockseis::File> Fout;
    Fout = this->getFdata();

    if ( Fout->is_open() )
    {
        // Update geometry
        size_t n2 = Fout->getN(2);
        Fout->setN(2, n2 + ntrace);
        Point3D<T> *scoords = (this->getGeom())->getScoords();
        Point3D<T> *gcoords = (this->getGeom())->getGcoords();
        for (int i=0; i < ntrace; i++)
        {
            //Write coordinates to data trace (First source x, y and z, and then receiver x, y and z)
            Fout->write(&scoords[i].x, 1);
            Fout->write(&scoords[i].y, 1);
            Fout->write(&scoords[i].z, 1);
            Fout->write(&gcoords[i].x, 1);
            Fout->write(&gcoords[i].y, 1);
            Fout->write(&gcoords[i].z, 1);
            // Write trace data
            Fout->write(&data[i*nt], nt);
        }	
    }else{
        status = FILE_ERR;	
    }

    if(Fout->getFail()){
        status = FILE_ERR;
    }else{
        status = FILE_OK;
    }
    return status;
}

template<typename T>
void Data3D<T>::createEmpty(size_t ntr)
{
    bool status;
    std::string datafile = this->getFile();
    if(datafile.empty()){
        rs_error("Data3D::createEmpty: No file assigned. ");
    }

    status = this->open("o");
    if(status == FILE_ERR) rockseis::rs_error("Data3D::createEmpty: Error opening file for writting");

    for (size_t i=0; i < ntr; i++){
        this->writeTraces();
    }
    this->close();
}

template<typename T>
void Data3D<T>::putTrace(std::string filename, size_t number)
{
    size_t traceno;
    traceno = number;

    bool status;
    std::shared_ptr<rockseis::File> Fdata (new rockseis::File());
    status = Fdata->append(filename);
    if(status == FILE_ERR) rs_error("Data3D::putTrace: Error opening output data file.");
    rs_datatype datatype = Fdata->getType(); 
    if(datatype != DATA3D) rs_error("Data3D::putTrace: Datafile must be of type Data3D.");
    //Get gather size information
    size_t n1 = Fdata->getN(1);
    T d1 = Fdata->getD(1);
    T o1 = Fdata->getO(1);
    size_t n2 = 1; // Only write one trace!
    if(n1 != this->getNt()) rs_error("Data3D::putTrace: Number of samples in trace and datafile mismatch.");
    if(d1 != this->getDt()) rs_error("Data3D::putTrace: Sampling interval in trace and datafile mismatch.");
    if(o1 != this->getOt()) rs_error("Data3D::putTrace: Origin in trace and datafile mismatch.");

    if(traceno > (Fdata->getN(2)-1)) rs_error("Data3D::putTrace: Trying to put a trace with number that is larger than number of traces in file. ", std::to_string(traceno));

    //Write gather
    Point3D<T> *scoords = (this->getGeom())->getScoords();
    Point3D<T> *gcoords = (this->getGeom())->getGcoords();
    T *tracedata = this->getData();
    for (size_t j=0; j < n2; j++){
        Fdata->seekp(Fdata->getStartofdata() + traceno*(n1+NHEAD3D)*sizeof(T));
        Fdata->write(&scoords[j].x, 1);
        Fdata->write(&scoords[j].y, 1);
        Fdata->write(&scoords[j].z, 1);
        Fdata->write(&gcoords[j].x, 1);
        Fdata->write(&gcoords[j].y, 1);
        Fdata->write(&gcoords[j].z, 1);
        Fdata->write(&tracedata[j*n1], n1);
    }

    if(Fdata->getFail()) rs_error("Data3D::PutTrace: Error writting gather to output file");
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

template class Data<double>;
template class Data2D<double>;
template class Data3D<double>;
}


