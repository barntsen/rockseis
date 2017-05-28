#include <iostream>
#include "model.h"
#include "pml.h"
#include "waves.h"
#include "utils.h"
#include "der.h"
#include "data.h"
#include "file.h"
#include <memory>
#include <fstream>
#include <math.h>



int main()
{
	// Parameters
	int lpml = 10;
	bool fs = 1; 

	std::shared_ptr<rockseis::File> Fmod (new rockseis::File());
	bool status;
	// Parameters from file
	int nx;
	int ny;
	int nz;
	float dx;
	float dy;
	float dz;
	float ox;
	float oy;
	float oz;

	// Set time variables for wavelet
	int nt;
	float dt;
	float ot;

	nx = 45;
	ny = 11;
	nz = 51;
	dx = 10.;
	dy = 10.;
	dz = 10.;
	ox = 0.;
	oy = 0.;
	oz = 0.;

	nt = 1001;
	ot = 0.;
	dt = 1e-3;

	// Create the classes 
	std::shared_ptr<rockseis::ModelElastic2D<float>> model2d (new rockseis::ModelElastic2D<float>(nx, nz, lpml, dx, dz, ox, oz, fs));
	std::shared_ptr<rockseis::ModelElastic3D<float>> model3d (new rockseis::ModelElastic3D<float>(nx, ny, nz, lpml, dx, dy, dz, ox, oy, oz, fs));


	/* Make a model */

    // Create empty models
    model2d->createModel();
    model3d->createModel();

    // Get pointers to the models
	float *R2d, *Vp2d, *Vs2d;
	float *R3d, *Vp3d, *Vs3d;
	R3d = model3d->getR();
	Vp3d = model3d->getVp();
	Vs3d = model3d->getVs();
	R2d = model2d->getR();
	Vp2d = model2d->getVp();
	Vs2d = model2d->getVs();
	rockseis::Index k2d(nx,nz);
	rockseis::Index k3d(nx,ny,nz);
	int ix,iy, iz;
    for(iz=0; iz<nz; iz++){
        for(iy=0; iy<ny; iy++){
            for(ix=0; ix<nx; ix++){
                R3d[k3d(ix,iy,iz)]= 1;
                Vp3d[k3d(ix,iy,iz)]= 2000;
                Vs3d[k3d(ix,iy,iz)]= 1100;
                R2d[k2d(ix,iz)]= 1;
                Vp2d[k2d(ix,iz)]= 2000;
                Vs2d[k2d(ix,iz)]= 1100;
            }
        }
    }

    // Adding a reflector
    for(iz=22; iz<nz; iz++){
        for(iy=0; iy<ny; iy++){
            for(ix=0; ix<nx; ix++){
                Vp2d[k2d(ix,iz)]= 2500;
                Vp3d[k3d(ix,iy,iz)]= 2500;
            }
        }
    }

    model2d->setVpfile("Vp2d.rss");
    model2d->setVsfile("Vs2d.rss");
    model2d->setRfile("Rho2d.rss");
    model2d->writeModel();

    model3d->setVpfile("Vp3d.rss");
    model3d->setVsfile("Vs3d.rss");
    model3d->setRfile("Rho3d.rss");
    model3d->writeModel();


	// Setup a 2d Wavelet 
	std::shared_ptr<rockseis::Data2D<float>> source2d (new rockseis::Data2D<float>(1, nt, dt, 0.0));
	float *wav = source2d->getData();
	float f0 = 20;
	float t0 = 1e-1;
	float t;
	float s=-3.1415*3.1415*f0*f0;
	rockseis::Index Idata(nt,1);
	for(int it=0; it < nt; it ++){
		// Computing wavelet 
		t = it * dt - t0;
		wav[Idata(it,0)] = (1.0 + 2.0*s*t*t)*exp(s*t*t); // Ricker function
		//wav[it] = sin(2*3.1415*f0*t);  // Harmonic function 
	}

 	// Set source2d position
	rockseis::Point2D<float> *scoords2d = (source2d->getGeom())->getScoords();
	scoords2d[0].x = 220;
	scoords2d[0].y = 10; 

	//Output the wavelet
	source2d->setFile("Wav2d.rss");
	status = source2d->write();
	if(status == FILE_ERR){
		std::cout << "Failed to write wavelet file. \n";
	}

   

	// Setup a 3d Wavelet 
	std::shared_ptr<rockseis::Data3D<float>> source3d (new rockseis::Data3D<float>(1, nt, dt, 0.0));
	wav = source3d->getData();
	for(int it=0; it < nt; it ++){
		// Computing wavelet 
		t = it * dt - t0;
		wav[it] = (1.0 + 2.0*s*t*t)*exp(s*t*t); // Ricker function
		//wav[it] = sin(2*3.1415*f0*t);  // Harmonic function 
	}

 	// Set source3d position
	rockseis::Point3D<float> *scoords3d = (source3d->getGeom())->getScoords();
	scoords3d[0].x = 220;
	scoords3d[0].y = 50; 
	scoords3d[0].z = 10; 

	//Output the wavelet
	source3d->setFile("Wav3d.rss");
	status = source3d->write();
	if(status == FILE_ERR){
		std::cout << "Failed to write wavelet file. \n";
	}

	// Setup a 2d record
	std::shared_ptr<rockseis::Data2D<float>> record2d (new rockseis::Data2D<float>(nx, 1, dt, 0.0));
	// Set receiver positions
	scoords2d = (record2d->getGeom())->getScoords();
	rockseis::Point2D<float> *gcoords2d = (record2d->getGeom())->getGcoords();
    for(int i=0; i<nx; i++)
    {
        gcoords2d[i].x = i*dx;
        gcoords2d[i].y = 10; 
        scoords2d[i].x = 220;
        scoords2d[i].y = 200;
    }
	record2d->setFile("Shot2d.rss");
	status = record2d->write();
	if(status == FILE_ERR){
		std::cout << "Failed to write record2d file. \n";
	}

	// Setup a 3d record
	std::shared_ptr<rockseis::Data3D<float>> record3d (new rockseis::Data3D<float>(nx*ny, 1, dt, 0.0));
	// Set receiver positions
	scoords3d = (record3d->getGeom())->getScoords();
	rockseis::Point3D<float> *gcoords3d = (record3d->getGeom())->getGcoords();
    rockseis::Index I(nx,ny);
    for(int j=0; j<ny; j++)
    {
        for(int i=0; i<nx; i++)
        {
            gcoords3d[I(i,j)].x = i*dx;
            gcoords3d[I(i,j)].y = j*dy; 
            gcoords3d[I(i,j)].z = 10; 
            scoords3d[I(i,j)].x = 220; 
            scoords3d[I(i,j)].y = 50;
            scoords3d[I(i,j)].z = 200;
        }
    }
	record3d->setFile("Shot3d.rss");
	status = record3d->write();
	if(status == FILE_ERR){
		std::cout << "Failed to write record3d file. \n";
	}


    /* Test Get local model function */
    std::shared_ptr<rockseis::ModelElastic2D<float>> lmodel2d;
    lmodel2d = model2d->getLocal(source2d, 900, SMAP);
    lmodel2d->setVpfile("LVp2d.rss");
    lmodel2d->setVsfile("LVs2d.rss");
    lmodel2d->setRfile("LRho2d.rss");
    lmodel2d->writeModel();

    std::shared_ptr<rockseis::ModelElastic3D<float>> lmodel3d;
    lmodel3d = model3d->getLocal(source3d, 900, 900, SMAP);
    lmodel3d->setVpfile("LVp3d.rss");
    lmodel3d->setVsfile("LVs3d.rss");
    lmodel3d->setRfile("LRho3d.rss");
    lmodel3d->writeModel();


    return 0;
}
