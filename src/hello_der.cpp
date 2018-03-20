#include <iostream>
#include <memory>
#include <fstream>
#include <math.h>
#include "data.h"
#include "der.h"
#include "model.h"
#include "file.h"
#include "sort.h"
#include "utils.h"
#include "geometry.h"

int main()
{
    std::string modelname, shotname;
    modelname = "Model.rss";
    shotname = "survey.rss";
    int lpml = 18;
    bool fs = 0;
    float apertx = 450.;
    float aperty = 450.;

    std::shared_ptr<rockseis::ModelAcoustic3D<float>> model (new rockseis::ModelAcoustic3D<float>(modelname, modelname, lpml, fs));
    std::shared_ptr<rockseis::ModelAcoustic3D<float>> lmodel (new rockseis::ModelAcoustic3D<float>(modelname, modelname, lpml, fs));

    std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->setDatafile(shotname);
    Sort->createShotmap(shotname); 

    std::shared_ptr<rockseis::Data3D<float>> Shotgeom;
    Shotgeom = Sort->get3DGather(0);
    lmodel = model->getLocal(Shotgeom, apertx, aperty, SMAP);

    std::shared_ptr<rockseis::Der<float>> der (new rockseis::Der<float>(lmodel->getNx_pml(), lmodel->getNy_pml(), lmodel->getNz_pml(), lmodel->getDx(), lmodel->getDy(), lmodel->getDz(), 8));
    lmodel->staggerModels();


    std::shared_ptr<rockseis::File> Fout (new rockseis::File());
    size_t nx = lmodel->getNx();
    size_t ny = lmodel->getNy();
    size_t nz = lmodel->getNz();

    float *array = lmodel->getVp();

    Fout->output("Vp.rss");
    Fout->setN(1,lmodel->getNx());
    Fout->setN(2,lmodel->getNy());
    Fout->setN(3,lmodel->getNz());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(array, nx*ny*nz, 0);
    Fout->close();

    array = lmodel->getL();
    float *df = der->getDf();

    nx = lmodel->getNx_pml();
    ny = lmodel->getNy_pml();
    nz = lmodel->getNz_pml();

    Fout->output("L.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(array, nx*ny*nz, 0);
    Fout->close();

    // Test all derivatives
    // ddx_fw
    der->ddx_fw(array);
    Fout->output("ddx_fw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();

    // ddx_fw
    der->ddx_bw(array);
    Fout->output("ddx_bw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();

// ddy_fw
    der->ddy_fw(array);
    Fout->output("ddy_fw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();

// ddy_bw
    der->ddy_bw(array);
    Fout->output("ddy_bw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();

// ddz_fw
    der->ddz_fw(array);
    Fout->output("ddz_fw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();

// ddz_bw
    der->ddz_bw(array);
    Fout->output("ddz_bw.rss");
    Fout->setN(1,lmodel->getNx_pml());
    Fout->setN(2,lmodel->getNy_pml());
    Fout->setN(3,lmodel->getNz_pml());
    Fout->setD(1,lmodel->getDx());
    Fout->setD(2,lmodel->getDy());
    Fout->setD(3,lmodel->getDz());
    Fout->setO(1,lmodel->getOx());
    Fout->setO(2,lmodel->getOy());
    Fout->setO(3,lmodel->getOz());
    Fout->setType(rockseis::REGULAR);
    Fout->setData_format(sizeof(float));
    Fout->writeHeader();
    Fout->write(df, nx*ny*nz, 0);
    Fout->close();
}
