#include <valarray>
#include <rsf.hh>
#include "file.h"

#define MAXDIM 8

int main(int argc, char* argv[])
{
    sf_init(argc,argv);
	std::shared_ptr<rockseis::File> out (new rockseis::File());
    out->output();

    madagascar::iRSF par(0);
    madagascar::iRSF in;

    int n[MAXDIM];
    float d[MAXDIM];
    float o[MAXDIM];
    for(int i=0; i< MAXDIM; i++){
        n[i]=0;
        d[i]=0.0;
        o[i]=0.0;
    }
    
    // Get parameters
    in.get("n1", n[0]);
    in.get("n2", n[1]);
    in.get("n3", n[2]);
    in.get("n4", n[3]);
    in.get("n5", n[4]);
    in.get("n6", n[5]);
    in.get("n7", n[6]);

    in.get("d1", d[0]);
    in.get("d2", d[1]);
    in.get("d3", d[2]);
    in.get("d4", d[3]);
    in.get("d5", d[4]);
    in.get("d6", d[5]);
    in.get("d7", d[6]);

    in.get("o1", o[0]);
    in.get("o2", o[1]);
    in.get("o3", o[2]);
    in.get("o4", o[3]);
    in.get("o5", o[4]);
    in.get("o6", o[5]);
    in.get("o7", o[6]);

    int filetype = in.type();
    
    if( n[0] < 1)
    {
        std::cout << "Error in input data.\n";
        exit(1);
    }

    for(int i=0; i < MAXDIM; i++)
    {
        if( n[i] )
        {
            out->setN(i+1,n[i]);
            out->setD(i+1,d[i]);
            out->setO(i+1,o[i]);
        }
    }
    out->printGeometry();

    // Find total number of traces in data
    size_t ntot = 1;
    for(int i=1; i < MAXDIM; i++)
    {
        if( n[i] ) ntot *= n[i];
    }

    // Allocate trace for reading RSF data
    size_t esize;
    if(filetype == 2) esize=sizeof(int);
    if(filetype == 3) esize=sizeof(float);

    // Make traces for both float and int types
    std::valarray<int> inttrace(n[0]);
    std::valarray<float> floattrace(n[0]);

    out->setData_format(esize);
    out->writeHeader();
    out->seekp(out->getStartofdata());
    for (int i=0; i < ntot; i++) {
        switch(filetype)
        {
            case 2:
                in >> inttrace;
                out->write(&(inttrace[0]), inttrace.size());
                break;
            case 3:
                in >> floattrace;
                out->write(&(floattrace[0]), floattrace.size());
                break;
            default:
                in >> floattrace;
                out->write(&(floattrace[0]), floattrace.size());
                break;
        }

    }
    out->close();


    exit (0);
}

