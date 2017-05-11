#include "file.h"
#include "data.h"
#include "inparse.h"
using namespace rockseis;

int main(int argc, char* argv[])
{
    /* Variables */
    unsigned long int fldr; 
    unsigned long int n2; 
	float sxpos;
	float sypos;
	float szpos;
	float gxpos;
	float gypos;
	float gzpos;
	float dx,dy;

	float sx0, sy0, sz0, dsx, dsy;
	float gx0, gy0, gz0, dgx, dgy;
	float ds1x,ds1y,ds2x,ds2y;
	int nsx, nsy, ngx, ngy;
	int isx, isy, igx, igy;
    int dim;
	bool OBC;
	bool verb;
    bool status;

    std::shared_ptr<rockseis::Data2D<float>> Bdata2d;
    std::shared_ptr<rockseis::Data3D<float>> Bdata3d;
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());

    if(argc < 2) rs_error("No configuration file given");


    /* Get parameters from input */
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("sx0", &sx0) == INPARSE_ERR) status = true;
    if(Inpar->getPar("sy0", &sy0) == INPARSE_ERR) status = true;
    if(Inpar->getPar("sz0", &sz0) == INPARSE_ERR) status = true;

    if(Inpar->getPar("dsx", &dsx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("dsy", &dsy) == INPARSE_ERR) status = true;

    if(Inpar->getPar("ds1x", &ds1x) == INPARSE_ERR) status = true;
    if(Inpar->getPar("ds1y", &ds1y) == INPARSE_ERR) status = true;

    if(Inpar->getPar("ds2x", &ds2x) == INPARSE_ERR) status = true;
    if(Inpar->getPar("ds2y", &ds2y) == INPARSE_ERR) status = true;

    if(Inpar->getPar("nsx", &nsx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nsy", &nsy) == INPARSE_ERR) status = true;

    if(Inpar->getPar("gx0", &gx0) == INPARSE_ERR) status = true;
    if(Inpar->getPar("gy0", &gy0) == INPARSE_ERR) status = true;
    if(Inpar->getPar("gz0", &gz0) == INPARSE_ERR) status = true;

    if(Inpar->getPar("dgx", &dgx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("dgy", &dgy) == INPARSE_ERR) status = true;

    if(Inpar->getPar("ngx", &ngx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("ngy", &ngy) == INPARSE_ERR) status = true;

    if(Inpar->getPar("dim", &dim) == INPARSE_ERR) status = true;
    if(Inpar->getPar("OBC", &OBC) == INPARSE_ERR) status = true;

    if(Inpar->getPar("verb", &verb) == INPARSE_ERR) status = true;

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    //Compute number of traces
    n2=nsx*nsy*ngx*ngy;
    rockseis::Point2D<float> *scoords2d;
    rockseis::Point2D<float> *gcoords2d;
    rockseis::Point3D<float> *scoords3d;
    rockseis::Point3D<float> *gcoords3d;

    if(dim == 3){
        Bdata3d = std::make_shared<rockseis::Data3D<float>>(1, 1, 1., 0.);
        Bdata3d->setFile("stdout");
        status = Bdata3d->open("o");
        if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
        scoords3d = (Bdata3d->getGeom())->getScoords();
        gcoords3d = (Bdata3d->getGeom())->getGcoords();
    }else {
        Bdata2d = std::make_shared<rockseis::Data2D<float>>(1, 1, 1., 0.);
        Bdata2d->setFile("stdout");
        status = Bdata2d->open("o");
        if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
        scoords2d = (Bdata2d->getGeom())->getScoords();
        gcoords2d = (Bdata2d->getGeom())->getGcoords();
    }

    fldr=0;
    for (isy=0; isy < nsy; isy++) {
        for (isx=0; isx < nsx; isx++) {

            // Increment fldr
            fldr++;
            if(fldr%2 == 0){
                dx=ds1x;
                dy=ds1y;
            }else{
                dx=ds2x;
                dy=ds2y;
            }

            // SX
            sxpos=sx0 + dx + isx*dsx;

            // SY
            sypos=sy0 + dy + isy*dsy;

            for (igy=0; igy < ngy; igy++) {
                // GY
                if( !OBC ){
                    gypos=sypos - dy + gy0 + igy*dgy;
                }else{
                    gypos=gy0 + igy*dgy;
                }
                for (igx=0; igx < ngx; igx++) {
                    // GX
                    if( !OBC ){
                        gxpos=sxpos - dx + gx0 + igx*dgx;
                    }else{
                        gxpos=gx0 + igx*dgx;
                    }
                    szpos= sz0;
                    gzpos= gz0;
                    if(dim == 3){
                        scoords3d[0].x = sxpos;
                        scoords3d[0].y = sypos;
                        scoords3d[0].z = szpos;
                        gcoords3d[0].x = gxpos;
                        gcoords3d[0].y = gypos;
                        gcoords3d[0].z = gzpos;
                        Bdata3d->writeTraces();
                    }else{
                        scoords2d[0].x = sxpos;
                        scoords2d[0].y = szpos;
                        gcoords2d[0].x = gxpos;
                        gcoords2d[0].y = gzpos;
                        Bdata2d->writeTraces();
                    }
                }
            }
        }
    }

    if(dim ==3){
        Bdata3d->close();
    }else{
        Bdata2d->close();
    }

    exit (0);
}

