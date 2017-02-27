#ifndef DATA_H
#define DATA_H

// Include statements
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "file.h"


namespace rockseis {
// =============== ABSTRACT DATA CLASS =============== //
/** The data class
 *
 */
template<typename T>
class Data{
public:
    Data(const std::string file); ///< Constructor
    Data(const int _ntrace, const int _nt, const T _dt); ///< Constructor
    Data(const int _ntrace, const int _nt, const T _dt, const T _ot); ///< Constructor
    ~Data(); ///< Destructor


    // Get functions
    int getNt() { return nt; }		///< Get nt
    T getDt() { return dt; }		///< Get dt
    T getOt() { return ot; }		///< Get ot
    int getNtrace() { return ntrace; }	///< Get number of traces
    std::string getFile() { return datafile; }

    // Set functions
    void setNt(const int _nt) { nt = _nt; }	///< Set Nx
    void setNtrace(const int _ntrace) { ntrace = _ntrace; }	///< Set number of traces
    void setDt(const T _dt) { dt = _dt; }	///< Set sampling interval (dt)
    void setOt(const T _ot) { ot = _ot; }	///< Set origin (ot)
    void setFile(std::string _datafile) { datafile=_datafile; }

private:
    std::string datafile;
    int ntrace;
    int nt;
    T dt;
    T ot;
};

// =============== 2D DATA CLASS =============== //
/** The 2D data class
 *
 */
template<typename T>
class Data2D: public Data<T> {
public:
    Data2D(const int _ntrace, const int _nt, const T _dt); 	///< Constructor
    Data2D(std::string datafile); ///< Constructor
    virtual ~Data2D();       	///< Destructor

    // Get functions
    std::shared_ptr<Geometry2D<T>> getGeom() { return geometry; } ///< Get data geometry
    T *getData() { return data; }  ///< Get data array

    // Make map function
    /** Makes an integer map of the data coordinates over a given model geometry
     * */
    void makeMap(std::shared_ptr<Geometry<T>> _geom) { geometry->makeMap(_geom); } 

    // File input/output
    /** Read and write data files */
    bool readData(); /// Read data and coordinates
    bool readCoords(); /// Read only coordinates
    bool write(); /// Writes data and coordinates


private:
    std::shared_ptr<Geometry2D<T>> geometry; // Data geometry 
    T *data; // Data array
};

// =============== 3D DATA GEOMETRY CLASS =============== //
/** The 3D data class
 *
 */
template<typename T>
class Data3D: public Data<T> {
public:
    Data3D(const int _ntrace, const int _nt, const T _dt); 	///< Constructor
    Data3D(std::string datafile); 	///< Constructor
    virtual ~Data3D();       	///< Destructor

    // Get functions
    std::shared_ptr<Geometry3D<T>> getGeom() { return geometry; } ///< Get data geometry
    T *getData() { return data; } ///< Get data array

    // Make map function
    /** Makes an integer map of the data coordinates over a given model geometry
     * */
    void makeMap(std::shared_ptr<Geometry<T>> _geom) { geometry->makeMap(_geom); }  

    // File input/output
    /** Read and write data files */
    bool readData(); /// Read data and coordinates
    bool readCoords(); /// Read only coordinates
    bool write(); /// Writes data and coordinates

private:
    std::shared_ptr<Geometry3D<T>> geometry;  // Data geometry 
    T *data; // Data array 
};

}

#endif //DATA_H
