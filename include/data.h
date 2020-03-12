#ifndef DATA_H
#define DATA_H

// Include statements
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "file.h"
#include "utils.h"
#include "fft.h"


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
    std::shared_ptr<File> getFdata() { return Fdata; } ///< Get Fdata handle
    int getNtrace() { return ntrace; }	///< Get number of traces
    int getRecinc() { return recinc; }	///< Get recording interval
    bool getAlloc() { return allocated; }	///< Get status if data is allocated
    std::string getFile() { return datafile; } ///< Get filename
    rs_field getField() { return field; } ///< Get data type
    bool open(std::string flag); ///< Open file. Flags can be "i", "o" or "a".
    T gauss(int n, int m); ///<Fourier transform of Gaussian
    void close();  ///< Close file

    // Set functions
    void setNt(const int _nt) { nt = _nt; }	///< Set Nx
    void setNtrace(const int _ntrace) { ntrace = _ntrace; }	///< Set number of traces
    void setDt(const T _dt) { dt = _dt; }	///< Set sampling interval (dt)
    void setOt(const T _ot) { ot = _ot; }	///< Set origin (ot)
    void setFile(std::string _datafile) { datafile=_datafile; } ///< Set filename
    void setFdata(std::shared_ptr<File> Fd) { Fdata = Fd; } ///< Set file handle
    void setField(rs_field _field) { field = _field; } ///< Set data type 
    void setRecinc(const int inc) { recinc = inc; }	///< Set Nx
    void setAlloc(bool val) {allocated = val; }	///< Set status if data is allocated

    // Filter
    void Filter1D(T f0, T f1, T f2, T f3, T df, unsigned long nf, T* W, T *cdata);

    // Hilbert transform
    void Hilbert1D(T *pulse, unsigned long nt);

    // S-transform
    void St1D(std::shared_ptr<rockseis::Fft<T>> fft1d, std::shared_ptr<rockseis::Fft<T>> ifft1d, T *data, int lo, int hi, T *result);
    void iSt1D(T *data, int lo, int hi, T *result);

private:
    std::string datafile;
    int ntrace;
    int nt;
    T dt;
    T ot;
    rs_field field;
    rs_datatype datatype;
    std::shared_ptr<File> Fdata;
    int recinc;
    bool allocated;
};

// =============== 2D DATA CLASS =============== //
/** The 2D data class
 *
 */
template<typename T>
class Data2D: public Data<T> {
public:
    Data2D(const int _ntrace, const int _nt, const T _dt, const T _ot); 	///< Constructor
    Data2D(std::string datafile, const int _nt, const T _dt, const T _ot); 	///< Constructor
    Data2D(std::string datafile); ///< Constructor
    Data2D(std::string datafile, bool allocate); ///< Constructor
    ~Data2D();       	///< Destructor

    // Get functions
    std::shared_ptr<Geometry2D<T>> getGeom() { return geometry; } ///< Get data geometry
    T *getData() { return data; }  ///< Get data array

    // Make map function
    /** Makes an integer map of the data coordinates over a given model geometry
     * */
    void makeMap(std::shared_ptr<Geometry<T>> geom) { geometry->makeMap(geom); } 
    void makeMap(std::shared_ptr<Geometry<T>> geom, bool map) { geometry->makeMap(geom, map); } 
    
    void copyGmap2Smap() { geometry->copyGmap2Smap(); } 
    void copySmap2Gmap() { geometry->copySmap2Gmap(); } 

    // File input/output
    /** Read and write data files */
    bool read(); ///< Read data and coordinates
    bool readCoords(); ///< Read only coordinates
    bool readTraces(); ///< Read a number of traces
    bool copyCoords(std::shared_ptr<Data2D<T>> Data); ///< copy coordinates from other data
    bool write(); ///< Writes data and coordinates
    bool writeTraces(); ///< Writes data and coordinates to the end of the file
    void createEmpty(size_t ntr); ///< Create empty dataset with ntr traces
    void getTrace(std::string filename, size_t number); ///< Writes a trace with a specific number to a file
    void putTrace(std::string filename, size_t number); ///< Writes a trace with a specific number to a file
    void putImage(std::string imagefile); ///< Output an image to the end of the datafile
    
    //Filter
    void apply_filter(T *freqs); ///< Apply a 4 point filter to the data

    //Hilbert transform
    void applyHilbert(); ///< Apply a Hilbert transform to the data

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
    Data3D(const int _ntrace, const int _nt, const T _dt, const T _ot); 	///< Constructor
    Data3D(std::string datafile, const int _nt, const T _dt, const T _ot); 	///< Constructor
    Data3D(std::string datafile); 	///< Constructor
    Data3D(std::string datafile, bool allocate); 	///< Constructor
    virtual ~Data3D();       	///< Destructor

    // Get functions
    std::shared_ptr<Geometry3D<T>> getGeom() { return geometry; } ///< Get data geometry
    T *getData() { return data; } ///< Get data array

    // Make map function
    /** Makes an integer map of the data coordinates over a given model geometry
     * */
    void makeMap(std::shared_ptr<Geometry<T>> geom) { geometry->makeMap(geom); }  
    void makeMap(std::shared_ptr<Geometry<T>> geom, bool map) { geometry->makeMap(geom, map); } 

    void copyGmap2Smap() { geometry->copyGmap2Smap(); } 
    void copySmap2Gmap() { geometry->copySmap2Gmap(); } 

    // File input/output
    /** Read and write data files */
    bool read(); /// Read data and coordinates
    bool readCoords(); /// Read only coordinates
    bool copyCoords(std::shared_ptr<Data3D<T>> Data); ///< copy coordinates from other data
    bool write(); /// Writes data and coordinates
    bool readTraces(); ///< Read a number of traces
    bool writeTraces(); /// Writes data and coordinates to the end of the file
    void createEmpty(size_t ntr); ///< Create empty dataset with ntr traces
    void putTrace(std::string filename, size_t number); ///< Writes a trace with a specific number to a file
    void putImage(std::string imagefile); ///< Output an image to the end of the datafile

    //Filter
    void apply_filter(T *freqs); ///< Apply a 4 point filter to the data

    //Hilbert transform
    void applyHilbert(); ///< Apply a Hilbert transform to the data

private:
    std::shared_ptr<Geometry3D<T>> geometry;  // Data geometry 
    T *data; // Data array 
};

}

#endif //DATA_H
