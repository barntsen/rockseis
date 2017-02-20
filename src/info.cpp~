#include "file.h"

#define MAXDIM 8

int main(int argc, char* argv[])
{
	std::shared_ptr<rockseis::File> in (new rockseis::File());

    bool status;
    status = in->input();
	if(status == FILE_ERR){
		std::cerr << "Error reading from input file. \n";
		exit(1);
	}
    in->printGeometry();
    std::cerr << "esize: " << in->getData_format() << "\n"; 
    std::cerr << "type: " << in->getType() << "\n"; 
    std::cerr << "Nheader: " << in->getNheader() << "\n"; 

    in->close();
    exit (0);
}

