#ifndef FILE_H
#define FILE_H

// Include statements
#include <string.h>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include "geometry.h"

// Define statements
#define MAGIC_NUMBER "R0CKS"
#define MAGIC_NUMBER_LENGTH 5
#define FILEFORMATVER 0
#define BUFFER 1000000

#define FILE_OK 1 
#define FILE_ERR 0


namespace rockseis {

// =============== ABSTRACT FILE CLASS =============== //
/** The file class
 *
 */

class File{
public:
	File(); ///< Default constructor
	~File(); ///< Destructor
	bool input(); ///< Open stdin for reading
	bool input(std::string filename); ///< Open a file for reading
	void output();  ///< Open standard output for writing 
	void output(std::string filename);  ///< Open a file for writing
	void append(); ///< Open standard output for writting without overwritting
	bool append(std::string filename); ///< Open a file for writting without overwritting
	void close();  ///< Close file
	
	// Write functions
	void writeHeader(); ///< Output header to file
	void seekp(size_t pos) { fstream.seekp(pos); }
	void seekg(size_t pos) { fstream.seekg(pos); }
	void seekp(size_t pos, std::ios_base::seekdir way) { fstream.seekp(pos, way); }
	void seekg(size_t pos, std::ios_base::seekdir way) { fstream.seekg(pos, way); }

	void write(char *buffer, size_t n); ///< Writes n chars to file
	void write(char *buffer, size_t n, size_t pos); ///< Writes n chars to file starting from pos
	void write(float *buffer, size_t n); ///< Writes n floats to file
	void write(float *buffer, size_t n, size_t pos); ///< Writes n floats to file starting from pos
	void write(double *buffer, size_t n); ///< Writes n doubles to file
	void write(double *buffer, size_t n, size_t pos); ///< Writes n doubles to file
	void write(int *buffer, size_t n); ///< Writes n ints to file
	void write(int *buffer, size_t n, size_t pos); ///< Writes n ints to file starting from pos
	void write(size_t *buffer, size_t n); ///< Writes n sze_t to file
	void write(size_t *buffer, size_t n, size_t pos); ///< Writes n size_t to file starting from pos
	bool is_open() { return fstream.is_open(); } ///< Checks if file is open
    bool checksum(); ///< Check if file size matches header and geometry information

	// Read functions
	void readHeader(); ///< Read header from file
	void read(char *buffer, size_t n); ///< Reads n chars from file
	void read(char *buffer, size_t n, size_t pos); ///< Reads n chars from file starting from pos
	void read(float *buffer, size_t n); ///< Reads n floats from file
	void read(float *buffer, size_t n, size_t pos); ///< Reads n floats from file starting from pos
	void read(double *buffer, size_t n); ///< Reads n doubles from file
	void read(double *buffer, size_t n, size_t pos); ///< Reads n doubles from file starting from pos
	void read(int *buffer, size_t n); ///< Reads n ints from file
	void read(int *buffer, size_t n, size_t pos); ///< Reads n ints from file starting from pos
	void read(size_t *buffer, size_t n); ///< Reads n size_ts from file
	void read(size_t *buffer, size_t n, size_t pos); ///< Reads n size_ts from file starting from pos

	// Get geometry
	std::shared_ptr<Geometry<double>> getGeom() { return geometry; }
	size_t getN(int dim) { return geometry->getN(dim); }  ///< Get dimension size
	double getD(int dim) { return geometry->getD(dim); }	///< Get sampling interval
	double getO(int dim) { return geometry->getO(dim); }	///< Get origin

    bool getHeaderstat() { return headerstat; } ///< Check if we need to write out header

	// Other geometry functions
	void printGeometry() { geometry->print(); }
	void clearGeometry() { geometry->clear(); }
	bool compareGeometry(std::shared_ptr<File> other) {return geometry->compare(other->getGeom()); }

	// Set geometry
        void setN(int dim, size_t val) { geometry->setN(dim, val); }	///< Set dimension size
        void setD(int dim, double val) { geometry->setD(dim, val); }	///< Set sampling interval
        void setO(int dim, double val) { geometry->setO(dim, val); }	///< Set origin

	// Set functions
	void setNheader(size_t val) { Nheader = val; } ///< set number of header values
	void setData_format(int val) { data_format = val; } ///< set data format 
	void setHeader_format(int val) { header_format = val; } ///< set header format 
	void setType(rs_datatype val) { type = val; } ///< set datatype
    void setHeaderstat(bool stat) { headerstat = stat; } ///< Set if header has been modified

	// Get functions
	size_t getNheader() { return Nheader; } ///< get number of header values
	int getData_format() { return data_format; } ///< get data format 
	int getHeader_format() { return header_format; } ///< get header format 
	rs_datatype getType() { return type; } ///< get datatype
	size_t getStartofdata() { return startofdata; } 

private:
	std::string magicnumber; ///< File format identifier
	std::string filename; ///< File to Read and/or Write to
	std::string mode; ///< Read, write, ASCII, binary 
	size_t version; ///< File format version
	std::fstream fstream;  ///< File Stream
	size_t Ndims; ///< Geometry information
	size_t Nheader; ///< Geometry information
	int data_format; ///< Size information on data (for binary data)
	int header_format; ///< Size information on header (for binary and irregular data)
    rs_datatype type;  ///< Information about the file content (Regular model, 2D data, 3D data, etc.)
	std::shared_ptr<Geometry<double>> geometry;  ///< regular geometry information from file (n1,n2, ...; d1, d2, ...; o1, o2, ...);
	size_t startofdata; ///<Position for seeking to start of data
    bool headerstat;

};

}

#endif //FILE_H
