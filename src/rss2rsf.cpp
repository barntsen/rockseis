#include <valarray>
#include <rsf.hh>
#include "file.h"

#define MAXDIM 8

int main(int argc, char* argv[])
{
    sf_init(argc,argv);
    std::shared_ptr<rockseis::File> in (new rockseis::File());
    bool status;
    status = in->input();
    if(status == FILE_ERR){
	    std::cout << "Error reading from input file. \n";
	    exit(1);
    }

    madagascar::oRSF out;

    int n[MAXDIM];
    float d[MAXDIM];
    float o[MAXDIM];
    for(int i=0; i< MAXDIM; i++){
        n[i]=0;
        d[i]=0.0;
        o[i]=0.0;
    }
    
    for(int i=0; i < MAXDIM; i++)
    {
	n[i] = in->getN(i+1);
	d[i] = (float) in->getD(i+1);
	o[i] = (float) in->getO(i+1);
    }

    // Get parameters
    if(n[0] > 0) out.put("n1", n[0]);
    if(n[1] > 0) out.put("n2", n[1]);
    if(n[2] > 0) out.put("n3", n[2]);
    if(n[3] > 0) out.put("n4", n[3]);
    if(n[4] > 0) out.put("n5", n[4]);
    if(n[5] > 0) out.put("n6", n[5]);
    if(n[6] > 0) out.put("n7", n[6]);

    if(n[0] > 0) out.put("d1", d[0]);
    if(n[1] > 0) out.put("d2", d[1]);
    if(n[2] > 0) out.put("d3", d[2]);
    if(n[3] > 0) out.put("d4", d[3]);
    if(n[4] > 0) out.put("d5", d[4]);
    if(n[5] > 0) out.put("d6", d[5]);
    if(n[6] > 0) out.put("d7", d[6]);

    if(n[0] > 0) out.put("o1", o[0]);
    if(n[1] > 0) out.put("o2", o[1]);
    if(n[2] > 0) out.put("o3", o[2]);
    if(n[3] > 0) out.put("o4", o[3]);
    if(n[4] > 0) out.put("o5", o[4]);
    if(n[5] > 0) out.put("o6", o[5]);
    if(n[6] > 0) out.put("o7", o[6]);

    // Find total number of traces in data
    size_t ntot = 1;
    for(int i=1; i < MAXDIM; i++)
    {
        if( n[i] ) ntot *= n[i];
    }

    // Allocate trace for reading RSF data
    int dsize, hsize;
    int filetype;
    sf_datatype type;
    dsize = in->getData_format();
    hsize = in->getHeader_format();
    if(dsize == sizeof(int))
    {
	    filetype = 2;
	    type = SF_INT;
    }
    if(dsize == sizeof(float))
    {
	    filetype = 3;
	    type = SF_FLOAT;
    }
    if(dsize == sizeof(double))
    {
	    filetype = 4;
	    type = SF_FLOAT; // We will convert to float on output
    }

    int Nheader = in->getNheader();

    // Make traces for both float and int types
    std::valarray<int> inttrace(n[0]);
    std::valarray<float> floattrace(n[0]);
    std::valarray<double> doubletrace(n[0]);

    out.type(type);
    for (int i=0; i < ntot; i++) {
        in->seekg(in->getStartofdata() + i*(n[0]*dsize+Nheader*hsize) + Nheader*hsize);
        switch(filetype)
        {
            case 2:
                in->read(&(inttrace[0]), inttrace.size());
                out << inttrace;
                break;
            case 3:
                in->read(&(floattrace[0]), floattrace.size());
                out << floattrace;
                break;
            case 4:
                in->read(&(doubletrace[0]), doubletrace.size());
                for (int j=0; j < n[0]; j++) floattrace[j] = (float) doubletrace[j];  // Converting from double to float
                out << floattrace;
                break;
            default:
                rockseis::rs_error("Unknown numerical precision in file.");
                break;
        }

    }

    in->close();

    exit (0);
}

