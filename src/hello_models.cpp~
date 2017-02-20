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
	nz = 41;
	dx = 10.;
	dy = 10.;
	dz = 10.;
	ox = 0.;
	oy = 0.;
	oz = 0.;

	nt = 501;
	ot = 0.;
	dt = 1e-3;

	// Create the classes 
	std::shared_ptr<rockseis::ModelElastic2D<float>> model2d (new rockseis::ModelElastic2D<float>(nx, nz, lpml, dx, dz, ox, oz, fs));
	std::shared_ptr<rockseis::ModelElastic3D<float>> model3d (new rockseis::ModelElastic3D<float>(nx, ny, nz, lpml, dx, dy, dz, ox, oy, oz, fs));

	// Make a model
	int ix,iy, iz;
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

	// Output 2d models
	Fmod->output("Vp2d.rss");
	Fmod->clearGeometry();
	Fmod->setN(1,nx);
	Fmod->setN(3,nz);
	Fmod->setD(1,dx);
	Fmod->setD(3,dz);
	Fmod->setO(1,ox);
	Fmod->setO(3,oz);
	Fmod->writeHeader();
	Fmod->floatwrite(Vp2d, nx*nz,0);
	Fmod->close();

	Fmod->output("Vs2d.rss");
	Fmod->writeHeader();
	Fmod->floatwrite(Vs2d, nx*nz,0);
	Fmod->close();

	Fmod->output("Rho2d.rss");
	Fmod->writeHeader();
	Fmod->floatwrite(R2d, nx*nz,0);
	Fmod->close();

	// Output 3d models
	Fmod->output("Vp3d.rss");
	Fmod->clearGeometry();
	Fmod->setN(1,nx);
	Fmod->setN(2,ny);
	Fmod->setN(3,nz);
	Fmod->setD(1,dx);
	Fmod->setD(2,dy);
	Fmod->setD(3,dz);
	Fmod->setO(1,ox);
	Fmod->setO(2,oy);
	Fmod->setO(3,oz);
	Fmod->writeHeader();
	Fmod->floatwrite(Vp3d, nx*ny*nz,0);
	Fmod->close();

	Fmod->output("Vs3d.rss");
	Fmod->writeHeader();
	Fmod->floatwrite(Vs3d, nx*ny*nz,0);
	Fmod->close();

	Fmod->output("Rho3d.rss");
	Fmod->writeHeader();
	Fmod->floatwrite(R3d, nx*ny*nz,0);
	Fmod->close();
	
	// Setup a 2d Wavelet 
	std::shared_ptr<rockseis::Data2D<float>> source2d (new rockseis::Data2D<float>(1, nt, dt));
	float *wav = source2d->getData();
	float f0 = 15;
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
	scoords2d[0].y = 200; 

	//Output the wavelet
	std::shared_ptr<rockseis::File> Fwav (new rockseis::File());
	Fwav->output("Wav2d.rss");
	Fwav->clearGeometry();
	status = source2d->writefloatData(Fwav);
	if(status == FILE_ERR){
		std::cout << "Failed to write wavelet file. \n";
	}
	Fwav->close();

	// Setup a 3d Wavelet 
	std::shared_ptr<rockseis::Data3D<float>> source3d (new rockseis::Data3D<float>(1, nt, dt));
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
	scoords3d[0].z = 200; 

	//Output the wavelet
	Fwav->output("Wav3d.rss");
	Fwav->clearGeometry();
	status = source3d->writefloatData(Fwav);
	if(status == FILE_ERR){
		std::cout << "Failed to write wavelet file. \n";
	}
	Fwav->close();



	// Setup a 2d record
	std::shared_ptr<rockseis::Data2D<float>> record2d (new rockseis::Data2D<float>(nx, 1, dt));
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
	Fwav->output("Geom2d.rss");
	Fwav->clearGeometry();
	status = record2d->writefloatData(Fwav);
	if(status == FILE_ERR){
		std::cout << "Failed to write record2d file. \n";
	}
	Fwav->close();

	// Setup a 3d record
	std::shared_ptr<rockseis::Data3D<float>> record3d (new rockseis::Data3D<float>(nx*ny, 1, dt));
	// Set receiver positions
	scoords3d = (record3d->getGeom())->getGcoords();
	rockseis::Point3D<float> *gcoords3d = (record3d->getGeom())->getGcoords();
    for(int j=0; j<ny; j++)
    {
        for(int i=0; i<nx; i++)
        {
            gcoords3d[j*ny + i].x = i*dx;
            gcoords3d[j*ny + i].y = j*dy; 
            gcoords3d[j*ny + i].z = 10; 
            scoords3d[j*ny + i].x = 220; 
            scoords3d[j*ny + i].y = 50;
            scoords3d[j*ny + i].z = 200;
        }
    }
    Fwav->output("Geom3d.rss");
	Fwav->clearGeometry();
	status = record3d->writefloatData(Fwav);
	if(status == FILE_ERR){
		std::cout << "Failed to write record3d file. \n";
	}
	Fwav->close();

	return 0;
}
