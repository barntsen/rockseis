#include "file.h"
#include "utils.h"

#define MAXDIM 8

int main(int argc, char* argv[])
{
	std::shared_ptr<rockseis::File> in (new rockseis::File());

    bool status;
    status = in->input();
	if(status == FILE_ERR){
        rockseis::rs_error("Error reading from input file.");
	}
    in->printGeometry();

    std::cerr << "esize: " << in->getData_format() << "\n"; 
    rockseis::rs_datatype type = static_cast<rockseis::rs_datatype>(in->getType());
    std::cerr << "type: ";
    switch(type){
        case rockseis::GENERIC:
            std::cerr << "GENERIC.";
            break;
        case rockseis::REGULAR:
            std::cerr << "REGULAR.";
            break;
        case rockseis::DATA2D:
            std::cerr << "DATA2D.";
            break;
        case rockseis::DATA3D:
            std::cerr << "DATA3D.";
            break;
        case rockseis::SNAPSHOT:
            std::cerr << "SNAPSHOT.";
            break;
        case rockseis::KEYMAP:
            std::cerr << "KEYMAP.";
            break;
        case rockseis::SORTMAP:
            std::cerr << "SORTMAP.";
            break;


    }
    std::cerr <<  std::endl;
    std::cerr << "Nheader: " << in->getNheader() << "\n"; 

    in->close();
    exit (0);
}
