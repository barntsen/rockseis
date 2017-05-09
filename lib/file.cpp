#include "file.h"

namespace rockseis {
// Dafault constructor
File::File()
{
    magicnumber = MAGIC_NUMBER;
    version = FILEFORMATVER;
    filename = "";
    Ndims = 1;
    Nheader = 0;
    data_format = -1;
    header_format = -1;
    type=GENERIC;
    geometry = std::make_shared<Geometry<double>>(); 
    headerstat = false;


    //Compute start of data
	startofdata=0;
	startofdata += magicnumber.size();
	startofdata += sizeof(int);
	startofdata += sizeof(int);
	startofdata += sizeof(int);
	startofdata += sizeof(size_t);
	startofdata += sizeof(size_t);
	for(int i=0; i<MAXDIMS; i++){
		startofdata += sizeof(size_t);
		startofdata += sizeof(double);
		startofdata += sizeof(double);
    }
}

bool File::input()
{
    bool status;
    char buffer[MAGIC_NUMBER_LENGTH+1];
    memset(buffer, 0, sizeof(buffer));
    fstream.open("/dev/stdin", std::ios::in | std::ios::binary);
    if(fstream.fail()){
        status = FILE_ERR; 
    }else{
        //Read header and check if file is of Rockseis format
        fstream.read(&buffer[0], MAGIC_NUMBER_LENGTH*sizeof(char));
        if(strcmp(buffer, MAGIC_NUMBER))
        {
            // Fail 
            status = FILE_ERR;
        }else{
            //Success
            geometry->clear();
            readHeader();

            status = checksum();
            if(status == FILE_ERR) rs_warning("Input file failed checksum check.");
        }
    }
    return status;
}

bool File::input(std::string filename)
{
    bool status;
    char buffer[MAGIC_NUMBER_LENGTH+1];
    memset(buffer, 0, sizeof(buffer));
    if(strcmp(filename.c_str(), "stdin")){
        fstream.open(filename, std::ios::in | std::ios::binary);
    }else{
        fstream.open("/dev/stdin", std::ios::in | std::ios::binary);
    }
    if(fstream.fail()){
        status = FILE_ERR; 
    }else{
        //Read header and check if file is of Rockseis format
        fstream.read(&buffer[0], MAGIC_NUMBER_LENGTH*sizeof(char));
        if(strcmp(buffer, MAGIC_NUMBER))
        {
            // Fail 
            status = FILE_ERR;
        }else{
            //Success
            geometry->clear();
            readHeader();
            status = checksum();
            if(status == FILE_ERR) rs_warning("File ", filename, " failed checksum check.");
        }
    }
    return status;
}

void File::output()
{
	fstream.open("/dev/stdout", std::ios::out | std::ios::binary);
    if(!fstream.is_open() || fstream.fail()) rs_error("Error opening standard output (stdout) for writting.");
}

void File::output(std::string filename)
{
    if(strcmp(filename.c_str(), "stdout")){
        fstream.open(filename, std::ios::out | std::ios::binary);
        if(!fstream.is_open() || fstream.fail()) rs_error("Error opening ", filename, " for writting.");
    }else{
        fstream.open("/dev/stdout", std::ios::out | std::ios::binary);
        if(!fstream.is_open() || fstream.fail()) rs_error("Error opening standard output (stdout) for writting.");
    }

}

void File::append()
{
    fstream.open("/dev/stdout", std::ios::app | std::ios::binary);
    if(!fstream.is_open() || fstream.fail()) rs_error("Error opening standard output (stdout) for appending.");
}

bool File::append(std::string filename)
{
    bool status;
    char buffer[MAGIC_NUMBER_LENGTH+1];
    size_t pos;
    memset(buffer, 0, sizeof(buffer));
	fstream.open(filename, std::ios::app | std::ios::binary);
    if(fstream.fail() || !fstream.is_open()){
        status = FILE_ERR; 
    }else{
        //Read header and check if file is of Rockseis format
        pos = fstream.tellp();
        fstream.seekp(0);
        fstream.read(&buffer[0], MAGIC_NUMBER_LENGTH*sizeof(char));
        if(strcmp(buffer, MAGIC_NUMBER))
        {
            // Fail 
            status = FILE_ERR;
        }else{
            //Success
            geometry->clear();
            readHeader();
            status = checksum();
            if(status == FILE_ERR) rs_warning("File ", filename, " failed checksum check.");
            fstream.seekp(pos);
        }
    }
    return status;
}

void File::close()
{
	fstream.close();
}

void File::writeHeader()
{
	//Get current position 
	size_t current_pos;
	current_pos=fstream.tellp();

	// Rewind file
	fstream.seekp(0);

	//Write header 
	fstream.write(magicnumber.c_str(), magicnumber.size());

	fstream.write(reinterpret_cast<char *> (&data_format), 1*sizeof(int));
	fstream.write(reinterpret_cast<char *> (&header_format), 1*sizeof(int));
	fstream.write(reinterpret_cast<char *> (&type), 1*sizeof(int));
	fstream.write(reinterpret_cast<char *> (&Nheader), 1*sizeof(size_t));
	fstream.write(reinterpret_cast<char *> (&Ndims), 1*sizeof(size_t));
	size_t i;
	size_t N;
	for(i=0; i<MAXDIMS; i++){
		N = geometry->getN(i+1);
		fstream.write(reinterpret_cast<char *> (&N), 1*sizeof(size_t));
	}
	double val;
	for(i=0; i<MAXDIMS; i++){
		val = geometry->getD(i+1);
		fstream.write(reinterpret_cast<char *> (&val), 1*sizeof(double));
	}
	for(i=0; i<MAXDIMS; i++){
		val = geometry->getO(i+1);
		fstream.write(reinterpret_cast<char *> (&val), 1*sizeof(double));
	}

	// Seek back to current position
	fstream.seekp(current_pos);
    // Reset header modification status
    this->setHeaderstat(false);
}

void File::readHeader()
{
    // Record number of bytes before startofdata

    // Rewind file
    fstream.seekg(magicnumber.size());

    //Read header 
    fstream.read(reinterpret_cast<char *> (&data_format), 1*sizeof(int));
    fstream.read(reinterpret_cast<char *> (&header_format), 1*sizeof(int));
    fstream.read(reinterpret_cast<char *> (&type), 1*sizeof(int));
    fstream.read(reinterpret_cast<char *> (&Nheader), 1*sizeof(size_t));
    fstream.read(reinterpret_cast<char *> (&Ndims), 1*sizeof(size_t));
    size_t i;
    size_t N;
    for(i=0; i<MAXDIMS; i++){
        fstream.read(reinterpret_cast<char *> (&N), 1*sizeof(size_t));
        geometry->setN(i+1,N);
    }
    double val;
    for(i=0; i<MAXDIMS; i++){
        fstream.read(reinterpret_cast<char *> (&val), 1*sizeof(double));
        geometry->setD(i+1, val);
    }
    for(i=0; i<MAXDIMS; i++){
        fstream.read(reinterpret_cast<char *> (&val), 1*sizeof(double));
        geometry->setO(i+1,val);
    }
    // Setting type
    this->setType(static_cast<rockseis::rs_datatype>(type));
}

bool File::checksum()
{
	size_t current_pos = fstream.tellg();
    size_t size_from_file, size_from_geom;
    size_t Ntot, N1;
    fstream.seekg(0, fstream.end);
    size_from_file = fstream.tellg();
    fstream.seekg(current_pos);
    Ntot= geometry->getNtot();
    N1= geometry->getN(1);
    size_from_geom =  startofdata + (Nheader*header_format*Ntot/N1) + (data_format*Ntot);
    if(size_from_file != size_from_geom){
        return FILE_ERR;
    }else {
        return FILE_OK;
    }
}

// Write functions
void File::write(char *buffer, size_t n)
{
	fstream.write(buffer, n);
}

void File::write(char *buffer, size_t n, size_t pos)
{
	fstream.seekp(pos + startofdata);
	fstream.write(buffer, n);
}

void File::write(float *buffer, size_t n)
{
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(float));
}

void File::write(float *buffer, size_t n, size_t pos)
{
	fstream.seekp(pos + startofdata);
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(float));
}

void File::write(double *buffer, size_t n)
{
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(double));
}

void File::write(double *buffer, size_t n, size_t pos)
{
	fstream.seekp(pos + startofdata);
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(double));
}

void File::write(int *buffer, size_t n)
{
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(int));
}

void File::write(int *buffer, size_t n, size_t pos)
{
	fstream.seekp(pos + startofdata);
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(int));
}

void File::write(size_t *buffer, size_t n)
{
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(size_t));
}

void File::write(size_t *buffer, size_t n, size_t pos)
{
	fstream.seekp(pos + startofdata);
	fstream.write(reinterpret_cast<char *> (buffer), n*sizeof(size_t));
}

//Read functions
void File::read(char *buffer, size_t n)
{
	fstream.read(buffer, n);
}

void File::read(char *buffer, size_t n, size_t pos)
{
	fstream.seekg(pos + startofdata);
	fstream.read(buffer, n);
}

void File::read(float *buffer, size_t n)
{
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(float));
}

void File::read(float *buffer, size_t n, size_t pos)
{
	fstream.seekg(pos + startofdata);
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(float));
}

void File::read(double *buffer, size_t n)
{
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(double));
}

void File::read(double *buffer, size_t n, size_t pos)
{
	fstream.seekg(pos + startofdata);
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(double));
}

void File::read(int *buffer, size_t n)
{
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(int));
}

void File::read(int *buffer, size_t n, size_t pos)
{
	fstream.seekg(pos + startofdata);
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(int));
}

void File::read(size_t *buffer, size_t n)
{
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(size_t));
}

void File::read(size_t *buffer, size_t n, size_t pos)
{
	fstream.seekg(pos + startofdata);
	fstream.read(reinterpret_cast<char *> (buffer), n*sizeof(size_t));
}


// destructor
File::~File(){
    fstream.close();
}
}


