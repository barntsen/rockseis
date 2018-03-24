#include "file.h"
#include "data.h"
#include "utils.h"
#include "inparse.h"

#define MAXDIM 8
#define ATTR 2
using namespace rockseis;

int main(int argc, char* argv[])
{
	std::shared_ptr<rockseis::File> in (new rockseis::File());

    bool status;

    /* Parameters */
    Point3D<float> s[ATTR];
    Point3D<float> g[ATTR];
    Point3D<float> offset[ATTR];
    Point3D<float> aoffset[ATTR];

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
        case rockseis::EDGESNAP:
            std::cerr << "EDGESNAP.";
            break;
        case rockseis::CHECKPOINT:
            std::cerr << "CHECKPOINT.";
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
    std::cerr << "Header format: " << in->getHeader_format() << "\n"; 

    if((type == rockseis::DATA2D || type == rockseis::DATA3D))
    {
        std::shared_ptr<rockseis::Data2D<float>> Bdata2d;
        std::shared_ptr<rockseis::Data3D<float>> Bdata3d;
        rockseis::Point2D<float> *scoords2d;
        rockseis::Point2D<float> *gcoords2d;
        rockseis::Point3D<float> *scoords3d;
        rockseis::Point3D<float> *gcoords3d;
        size_t ntr;
        switch(type)
        {
            case DATA2D:
                Bdata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Bdata2d->setFdata(in);
                scoords2d = (Bdata2d->getGeom())->getScoords();
                gcoords2d = (Bdata2d->getGeom())->getGcoords();

                if(Bdata2d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                s[0].x = scoords2d[0].x;
                s[1].x = scoords2d[0].x;
                g[0].x = gcoords2d[0].x;
                g[1].x = gcoords2d[0].x;
                s[0].y = scoords2d[0].y;
                s[1].y = scoords2d[0].y;
                g[0].y = gcoords2d[0].y;
                g[1].y = gcoords2d[0].y;               
                offset[0].x = gcoords2d[0].x - scoords2d[0].x;
                offset[1].x = gcoords2d[0].x - scoords2d[0].x;
                offset[0].y = gcoords2d[0].y - scoords2d[0].y;
                offset[1].y = gcoords2d[0].y - scoords2d[0].y;

                aoffset[0].x = fabsf(gcoords2d[0].x - scoords2d[0].x);
                aoffset[1].x = fabsf(gcoords2d[0].x - scoords2d[0].x);
                aoffset[0].y = fabsf(gcoords2d[0].y - scoords2d[0].y);
                aoffset[1].y = fabsf(gcoords2d[0].y - scoords2d[0].y);

                ntr = in->getN(2);
                for(size_t i=1; i< ntr; i++)
                {
                    if(Bdata2d->readTraces() == FILE_ERR) rs_error("Error reading from input file");

                    if (scoords2d[0].x < s[0].x ) s[0].x = scoords2d[0].x;
                    if (scoords2d[0].x > s[1].x ) s[1].x = scoords2d[0].x;
                    if (gcoords2d[0].x < g[0].x ) g[0].x = gcoords2d[0].x;
                    if (gcoords2d[0].x > g[1].x ) g[1].x = gcoords2d[0].x;
                    if (scoords2d[0].y < s[0].y ) s[0].y = scoords2d[0].y;
                    if (scoords2d[0].y > s[1].y ) s[1].y = scoords2d[0].y;
                    if (gcoords2d[0].y < g[0].y ) g[0].y = gcoords2d[0].y;
                    if (gcoords2d[0].y > g[1].y ) g[1].y = gcoords2d[0].y;
                    if((gcoords2d[0].x - scoords2d[0].x) < offset[0].x) offset[0].x = gcoords2d[0].x - scoords2d[0].x;
                    if((gcoords2d[0].x - scoords2d[0].x) > offset[1].x) offset[1].x = gcoords2d[0].x - scoords2d[0].x;
                    if((gcoords2d[0].y - scoords2d[0].y) < offset[0].y) offset[0].y = gcoords2d[0].y - scoords2d[0].y;
                    if((gcoords2d[0].y - scoords2d[0].y) > offset[1].y) offset[1].y = gcoords2d[0].y - scoords2d[0].y;
                    if(fabsf(gcoords2d[0].x - scoords2d[0].x) < aoffset[0].x) aoffset[0].x = fabsf(gcoords2d[0].x - scoords2d[0].x);
                    if(fabsf(gcoords2d[0].x - scoords2d[0].x) > aoffset[1].x) aoffset[1].x = fabsf(gcoords2d[0].x - scoords2d[0].x);
                    if(fabsf(gcoords2d[0].y - scoords2d[0].y) < aoffset[0].y) aoffset[0].y = fabsf(gcoords2d[0].y - scoords2d[0].y);
                    if(fabsf(gcoords2d[0].y - scoords2d[0].y) > aoffset[1].y) aoffset[1].y = fabsf(gcoords2d[0].y - scoords2d[0].y);

                }
                std::cerr << std::endl;
                std::cerr <<"Atributes" << std::endl;
                std::cerr <<"*******************************************************************************" << std::endl;
                std::cerr <<"key      \t\t\t\t" << "min\t\t\t\t" << "max\t\t\t\t\t" << std::endl;
                std::cerr <<"-------------------------------------------------------------------------------" << std::endl;
                std::cerr <<"sx:      \t\t\t\t" << s[0].x << "\t\t\t\t" << s[1].x << std::endl;
                std::cerr <<"sz:      \t\t\t\t" << s[0].y << "\t\t\t\t" << s[1].y << std::endl;

                std::cerr <<"gx:      \t\t\t\t" << g[0].x << "\t\t\t\t" << g[1].x << std::endl;
                std::cerr <<"gz:      \t\t\t\t" << g[0].y << "\t\t\t\t" << g[1].y << std::endl;

                std::cerr <<"offsetx: \t\t\t\t" << offset[0].x << "\t\t\t\t" << offset[1].x << std::endl;
                std::cerr <<"offsetz: \t\t\t\t" << offset[0].y << "\t\t\t\t" << offset[1].y << std::endl;
                std::cerr <<"Aoffsetx:\t\t\t\t" << aoffset[0].x << "\t\t\t\t" << aoffset[1].x << std::endl;
                std::cerr <<"Aoffsetz:\t\t\t\t" << aoffset[0].y << "\t\t\t\t" << aoffset[1].y << std::endl;
                std::cerr <<"*******************************************************************************" << std::endl;
                break;
            case DATA3D:
                Bdata3d = std::make_shared<rockseis::Data3D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Bdata3d->setFdata(in);
                scoords3d = (Bdata3d->getGeom())->getScoords();
                gcoords3d = (Bdata3d->getGeom())->getGcoords();

                if(Bdata3d->readTraces() == FILE_ERR)  rs_error("Error reading from input file ");
                s[0].x = scoords3d[0].x;
                s[1].x = scoords3d[0].x;
                g[0].x = gcoords3d[0].x;
                g[1].x = gcoords3d[0].x;
                s[0].y = scoords3d[0].y;
                s[1].y = scoords3d[0].y;
                g[0].y = gcoords3d[0].y;
                g[1].y = gcoords3d[0].y;
                s[0].z = scoords3d[0].z;
                s[1].z = scoords3d[0].z;
                g[0].z = gcoords3d[0].z;
                g[1].z = gcoords3d[0].z;
                offset[0].x = gcoords3d[0].x - scoords3d[0].x;
                offset[1].x = gcoords3d[0].x - scoords3d[0].x;
                offset[0].y = gcoords3d[0].y - scoords3d[0].y;
                offset[1].y = gcoords3d[0].y - scoords3d[0].y;
                offset[0].z = gcoords3d[0].z - scoords3d[0].z;
                offset[1].z = gcoords3d[0].z - scoords3d[0].z;

                aoffset[0].x = fabsf(gcoords3d[0].x - scoords3d[0].x);
                aoffset[1].x = fabsf(gcoords3d[0].x - scoords3d[0].x);
                aoffset[0].y = fabsf(gcoords3d[0].y - scoords3d[0].y);
                aoffset[1].y = fabsf(gcoords3d[0].y - scoords3d[0].y);
                aoffset[0].z = fabsf(gcoords3d[0].z - scoords3d[0].z);
                aoffset[1].z = fabsf(gcoords3d[0].z - scoords3d[0].z);

                ntr = in->getN(2);
                for(size_t i=1; i< ntr; i++)
                {
                    if(Bdata3d->readTraces() == FILE_ERR) rs_error("Error reading from input file");

                    if (scoords3d[0].x < s[0].x ) s[0].x = scoords3d[0].x;
                    if (scoords3d[0].x > s[1].x ) s[1].x = scoords3d[0].x;
                    if (gcoords3d[0].x < g[0].x ) g[0].x = gcoords3d[0].x;
                    if (gcoords3d[0].x > g[1].x ) g[1].x = gcoords3d[0].x;
                    if (scoords3d[0].y < s[0].y ) s[0].y = scoords3d[0].y;
                    if (scoords3d[0].y > s[1].y ) s[1].y = scoords3d[0].y;
                    if (gcoords3d[0].y < g[0].y ) g[0].y = gcoords3d[0].y;
                    if (gcoords3d[0].y > g[1].y ) g[1].y = gcoords3d[0].y;
                    if (scoords3d[0].z < s[0].z ) s[0].z = scoords3d[0].z;
                    if (scoords3d[0].z > s[1].z ) s[1].z = scoords3d[0].z;
                    if (gcoords3d[0].z < g[0].z ) g[0].z = gcoords3d[0].z;
                    if (gcoords3d[0].z > g[1].z ) g[1].z = gcoords3d[0].z;
                    if((gcoords3d[0].x - scoords3d[0].x) < offset[0].x) offset[0].x = gcoords3d[0].x - scoords3d[0].x;
                    if((gcoords3d[0].x - scoords3d[0].x) > offset[1].x) offset[1].x = gcoords3d[0].x - scoords3d[0].x;
                    if((gcoords3d[0].y - scoords3d[0].y) < offset[0].y) offset[0].y = gcoords3d[0].y - scoords3d[0].y;
                    if((gcoords3d[0].y - scoords3d[0].y) > offset[1].y) offset[1].y = gcoords3d[0].y - scoords3d[0].y;
                    if((gcoords3d[0].z - scoords3d[0].z) < offset[0].z) offset[0].z = gcoords3d[0].z - scoords3d[0].z;
                    if((gcoords3d[0].z - scoords3d[0].z) > offset[1].z) offset[1].z = gcoords3d[0].z - scoords3d[0].z;
                    if(fabsf(gcoords3d[0].x - scoords3d[0].x) < aoffset[0].x) aoffset[0].x = fabsf(gcoords3d[0].x - scoords3d[0].x);
                    if(fabsf(gcoords3d[0].x - scoords3d[0].x) > aoffset[1].x) aoffset[1].x = fabsf(gcoords3d[0].x - scoords3d[0].x);
                    if(fabsf(gcoords3d[0].y - scoords3d[0].y) < aoffset[0].y) aoffset[0].y = fabsf(gcoords3d[0].y - scoords3d[0].y);
                    if(fabsf(gcoords3d[0].y - scoords3d[0].y) > aoffset[1].y) aoffset[1].y = fabsf(gcoords3d[0].y - scoords3d[0].y);
                    if(fabsf(gcoords3d[0].z - scoords3d[0].z) < aoffset[0].z) aoffset[0].z = fabsf(gcoords3d[0].z - scoords3d[0].z);
                    if(fabsf(gcoords3d[0].z - scoords3d[0].z) > aoffset[1].z) aoffset[1].z = fabsf(gcoords3d[0].z - scoords3d[0].z);

                }
                std::cerr << std::endl;
                std::cerr <<"Atributes" << std::endl;
                std::cerr <<"*******************************************************************************" << std::endl;
                std::cerr <<"key      \t\t\t\t" << "min\t\t\t\t" << "max\t\t\t\t\t" << std::endl;
                std::cerr <<"-------------------------------------------------------------------------------" << std::endl;
                std::cerr <<"sx:      \t\t\t\t" << s[0].x << "\t\t\t\t" << s[1].x << std::endl;
                std::cerr <<"sy:      \t\t\t\t" << s[0].y << "\t\t\t\t" << s[1].y << std::endl;
                std::cerr <<"sz:      \t\t\t\t" << s[0].z << "\t\t\t\t" << s[1].z << std::endl;

                std::cerr <<"gx:      \t\t\t\t" << g[0].x << "\t\t\t\t" << g[1].x << std::endl;
                std::cerr <<"gy:      \t\t\t\t" << g[0].y << "\t\t\t\t" << g[1].y << std::endl;
                std::cerr <<"gz:      \t\t\t\t" << g[0].z << "\t\t\t\t" << g[1].z << std::endl;

                std::cerr <<"offsetx: \t\t\t\t" << offset[0].x << "\t\t\t\t" << offset[1].x << std::endl;
                std::cerr <<"offsety: \t\t\t\t" << offset[0].y << "\t\t\t\t" << offset[1].y << std::endl;
                std::cerr <<"offsetz: \t\t\t\t" << offset[0].z << "\t\t\t\t" << offset[1].z << std::endl;

                std::cerr <<"Aoffsetx:\t\t\t\t" << aoffset[0].x << "\t\t\t\t" << aoffset[1].x << std::endl;
                std::cerr <<"Aoffsety:\t\t\t\t" << aoffset[0].y << "\t\t\t\t" << aoffset[1].y << std::endl;
                std::cerr <<"Aoffsetz:\t\t\t\t" << aoffset[0].z << "\t\t\t\t" << aoffset[1].z << std::endl;
                std::cerr <<"*******************************************************************************" << std::endl;
                break;
            default:
                break;
        }
    }

    in->close();
    exit (0);
}
