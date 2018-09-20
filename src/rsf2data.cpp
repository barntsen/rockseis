#include <valarray>
#include <rsf.hh>
#include "file.h"
#include "data.h"
#include "utils.h"
#include <args.hxx>

#define NTFILE 91
#define MAXDIM 8

int main(int argc, char* argv[])
{

    std::string hdrfile, Datatype;
    args::ArgumentParser parser("Program to convert RSF file to Data RSS file.", "");
    parser.LongPrefix("");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> parhdrfile(parser, "FILENAME", "Input tfile", {"tfile"});
    args::ValueFlag<std::string> partype(parser, "DATA2D or DATA3D", "Dimension", {"type"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cerr << parser;
        return 1;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (parhdrfile){
         hdrfile = args::get(parhdrfile);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing tfile.");
    }
    if (partype){
         Datatype = args::get(partype);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing type.");
    }


    sf_init(argc,argv);
    madagascar::iRSF par(0);
    madagascar::iRSF in;
    madagascar::iRSF *tfile;
    tfile = new madagascar::iRSF(hdrfile.c_str());

    off_t n[MAXDIM];
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

    if(in.type() != SF_FLOAT) rockseis::rs_error("File type not supported.");
    if((*tfile).type() != SF_INT) rockseis::rs_error("tfile file type must be SF_INT.");
    int nh;
    (*tfile).get("n1", nh);
    if(nh != NTFILE) rockseis::rs_error("Wrong shape in tfile, n1 should be 91");

    if( n[0] < 1)
    {
        std::cerr << "Error in input data.\n";
        exit(1);
    }

    /* Variables */
    int type = -1; 
    if(!Datatype.compare("DATA2D"))
    {
        type = 2;
    }
    if(!Datatype.compare("DATA3D"))
    {
        type = 3;
    }
    if(type == -1)
    {
        rockseis::rs_error("type must be either DATA2D or DATA3D");
    }

    rockseis::rs_datatype dtype = static_cast<rockseis::rs_datatype>(type);

    // Getting filesize
    if (n[1] == 0) n[1] = 1;
    long ntr = n[1];
    for (int i = 2; i < MAXDIM; i++){
        if(n[i] > 0) ntr *= n[i];
    }

    off_t nhx;
    (*tfile).get("n2", nhx);
    if(nhx != ntr) rockseis::rs_error("Mismatch in number of traces between data file and tfile.");

    // Rockseis data 
    std::shared_ptr<rockseis::Data2D<float>> Bdata2D = NULL;
    std::shared_ptr<rockseis::Data3D<float>> Bdata3D = NULL;
    rockseis::Point2D<float> *scoords2D = NULL;
    rockseis::Point2D<float> *gcoords2D = NULL;
    rockseis::Point3D<float> *scoords3D = NULL;
    rockseis::Point3D<float> *gcoords3D = NULL;
    float *ftrace;

    // Madagascar data
    std::valarray<float> floattrace(n[0]);
    std::valarray<int> headtrace(NTFILE);
    double fscalco, fscalel;
    double sx, sy, sz, gx, gy, gz;
    int status;
    switch(dtype){
        case rockseis::DATA2D:
            Bdata2D = std::make_shared<rockseis::Data2D<float>>(1, n[0], d[0], o[0]);
            scoords2D = (Bdata2D->getGeom())->getScoords();
            gcoords2D = (Bdata2D->getGeom())->getGcoords();
            //Using stdout as output 
            Bdata2D->setFile("stdout");
            // Open data for output
            status = Bdata2D->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

            ftrace = Bdata2D->getData();
            for (size_t i=0; i < ntr; i++)
            {
                in >> floattrace;
                *tfile >> headtrace;
                if(headtrace[19] < 0){
                    fscalel = -1.0/headtrace[19];
                }else{
                    fscalel = headtrace[19];
                }
                if(headtrace[20] < 0){
                    fscalco = -1.0/headtrace[20];
                }else{
                    fscalco = headtrace[20];
                }
                sx = (double) headtrace[21];
                sz = (double) headtrace[13];

                gx = (double) headtrace[23];
                gz = (double) headtrace[12];

                scoords2D[0].x = fscalco * sx;
                scoords2D[0].y = fscalel * sz;
                gcoords2D[0].x = fscalco * gx;
                gcoords2D[0].y = fscalel * gz;
                for(long j=0; j < n[0]; j++) {
                    ftrace[j] = floattrace[j];
                }
                Bdata2D->writeTraces();
            }	
            Bdata2D->close();
            break;
        case rockseis::DATA3D:
            Bdata3D = std::make_shared<rockseis::Data3D<float>>(1, n[0], d[0], o[0]);
            scoords3D = (Bdata3D->getGeom())->getScoords();
            gcoords3D = (Bdata3D->getGeom())->getGcoords();
            //Using stdout as output 
            Bdata3D->setFile("stdout");
            // Open data for output
            status = Bdata3D->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

            ftrace = Bdata3D->getData();
            for (size_t i=0; i < ntr; i++)
            {
                in >> floattrace;
                *tfile >> headtrace;
                if(headtrace[19] < 0){
                    fscalel = -1.0/headtrace[19];
                }else{
                    fscalel = headtrace[19];
                }
                if(headtrace[20] < 0){
                    fscalco = -1.0/headtrace[20];
                }else{
                    fscalco = headtrace[20];
                }
                sx = (double) headtrace[21];
                sy = (double) headtrace[22];
                sz = (double) headtrace[13];

                gx = (double) headtrace[23];
                gy = (double) headtrace[24];
                gz = (double) headtrace[12];

                scoords3D[0].x = fscalco * sx;
                scoords3D[0].y = fscalco * sy;
                scoords3D[0].z = fscalel * sz;
                gcoords3D[0].x = fscalco * gx;
                gcoords3D[0].y = fscalco * gy;
                gcoords3D[0].z = fscalel * gz;
                for(long j=0; j < n[0]; j++) {
                    ftrace[j] = floattrace[j];
                }
                Bdata3D->writeTraces();
            }	
            Bdata3D->close();
            break;
        default:
            rockseis::rs_error("Invalid data type.");
            break;
    }

    exit (0);
}

