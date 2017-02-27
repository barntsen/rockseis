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
	bool fs = 0; 
	int order=2;

	std::shared_ptr<rockseis::File> file (new rockseis::File());
	bool status;
	status = file->input("Vp3d.rss");
	if(status == FILE_ERR){
		std::cout << "Error reading from Vp file. \n";
	}
	file->close();
	// Parameters from file
	int nx = file->getN(1);
	int ny = file->getN(2);
	int nz = file->getN(3);
	float dx = (float) file->getD(1);
	float dy = (float) file->getD(2);
	float dz = (float) file->getD(3);
	float ox = (float) file->getO(1);
	float oy = (float) file->getO(2);
	float oz = (float) file->getO(3);

	// Get time variables from wav file
	int nt;
	float dt;
	float ot;
	status = file->input("Wav3d.rss");
	if(status == FILE_ERR){
		std::cout << "Error reading from wavelet file. \n";
	}
	file->close();
	nt = file->getN(1);
	dt = (float) file->getD(1);
	ot = (float) file->getO(1);

	/*
	nx = 45;
	ny = 11;
	nz = 41;
	dx = 10.;
	dy = 10.;
	dz = 10.;
	ox = 0.;
	oy = 0.;
	oz = 0.;
	*/

/*	nt = 501;
	ot = 0.;
	dt = 1e-3;
	*/

	std::cout << "Nt = " << nt << "\n";
	std::cout << "Dt = " << dt << "\n";
	std::cout << "Ot = " << ot << "\n";

	std::cout << "Nx = " << nx << "\n";
	std::cout << "Ny = " << ny << "\n";
	std::cout << "Nz = " << nz << "\n";
	std::cout << "Dx = " << dx << "\n";
	std::cout << "Dy = " << dy << "\n";
	std::cout << "Dz = " << dz << "\n";
	std::cout << "Ox = " << ox << "\n";
	std::cout << "Oy = " << oy << "\n";
	std::cout << "Oz = " << oz << "\n";


	// Create the classes 
	std::shared_ptr<rockseis::WavesElastic3D<float>> waves (new rockseis::WavesElastic3D<float>(nx, ny, nz, nt, lpml, dx, dy, dz, dt, ox, oy, oz, ot));
	std::shared_ptr<rockseis::ModelElastic3D<float>> model (new rockseis::ModelElastic3D<float>(nx, ny, nz, lpml, dx, dy, dz, ox, oy, oz, fs));
	std::shared_ptr<rockseis::Der<float>> der (new rockseis::Der<float>(nx+2*lpml, ny+2*lpml, nz+2*lpml, dx, dy, dz, order));
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(1, nt, dt));


	/*
	// Make an acoustic model
	int ix,iy, iz;
	float *R, *Vp, *Vs;
	R = model->getR();
	Vp = model->getVp();
	Vs = model->getVs();
	rockseis::Index k1(nx,ny,nz);
	for(iz=0; iz<nz; iz++){
		for(iy=0; iy<ny; iy++){
		for(ix=0; ix<nx; ix++){
			R[k1(ix,iy,iz)]= 1;
			Vp[k1(ix,iy,iz)]= 2000;
			Vs[k1(ix,iy,iz)]= 1100;
		}
	}
	}

	file->output("Vp3d.rss");
	file->setN(1,nx);
	file->setN(2,ny);
	file->setN(3,nz);
	file->setD(1,dx);
	file->setD(2,dy);
	file->setD(3,dz);
	file->setO(1,ox);
	file->setO(2,oy);
	file->setO(3,oz);
	file->writeHeader();
	file->write(Vp, nx*ny*nz,0);
	file->close();

	file->output("Vs3d.rss");
	file->writeHeader();
	file->write(Vs, nx*ny*nz,0);
	file->close();

	file->output("Rho3d.rss");
	file->writeHeader();
	file->write(R, nx*ny*nz,0);
	file->close();

	*/
	
	// Get models from files
	float *R, *Vp, *Vs;
	R = model->getR();
	Vp = model->getVp();
	Vs = model->getVs();
	// Read vp model
	status = file->input("Vp3d.rss");
	file->read(Vp, nx*ny*nz);
	file->close();

	// Read vs model
	status = file->input("Vs3d.rss");
	file->read(Vs, nx*ny*nz);
	file->close();

	// Read rho model
	status = file->input("Rho3d.rss");
	file->read(R, nx*ny*nz);
	file->close();

	// Stagger model
	model->staggerModels();

	/*
	// Setup a Wavelet
	float *wav = source->getData();
	float f0 = 15;
	float t0 = 1e-1;
	float t;
	float s=-3.1415*3.1415*f0*f0;
	for(int it=0; it < nt; it ++){
		// Computing wavelet 
		t = it * dt - t0;
		wav[it] = (1.0 + 2.0*s*t*t)*exp(s*t*t); // Ricker function
		//wav[it] = sin(2*3.1415*f0*t);  // Harmonic function 
	}

 	// Set source position
	rockseis::Point3D<float> *coords = (source->getGeom())->getScoords();
	coords[0].x = 220;
	coords[0].y = 50; 
	coords[0].z = 200; 
	source->makeMap(model->getGeom());
	rockseis::Point3D<int> *map = (source->getGeom())->getSmap();
	//Output the wavelet
	file->output("Wav3d.rss");
	status = source->writefloatData(file);
	if(status == FILE_ERR){
		std::cout << "Failed to write wavelet file. \n";
	}
	file->close();
	*/

	// Load wavelet from file 
	status = file->input("Wav3d.rss");
	source->readData();
	file->close();
	source->makeMap(model->getGeom());

	// Output snapshots to a binary file
	std::ofstream myFile;
        myFile.open ("data.bin", std::ios::out | std::ios::binary);

	// Pointer to one of the Fields
	float *Sxx;
	Sxx = waves->getSxx();

	// Loop over time
	for(int it=0; it < nt; it++)
	{
		// Time stepping
		waves->forwardstepVelocity(model, der);
		waves->forwardstepStress(model, der);

		// Inserting source
		waves->insertSource(model, source, 3, 0, it);

		//Writting out results to binary file
		myFile.write ((char *) Sxx, (nx+2*lpml) * (ny+2*lpml) * (nz+2*lpml) * sizeof(float));
		
	}	

	myFile.close();

	return 0;
}
