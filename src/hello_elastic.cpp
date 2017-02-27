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
#include <config4cpp/Configuration.h>

int main()
{
	bool status;
	std::shared_ptr<rockseis::File> file (new rockseis::File());
	// Get geometry variables from Vp file
	status = file->input("Vp2d.rss");
	if(status == FILE_ERR){
		std::cout << "Error reading from Vp file. \n";
	}
	file->close();
	
	// Parameters from file
	int nx = file->getN(1);
	int nz = file->getN(3);
	float dx = (float) file->getD(1);
	float dz = (float) file->getD(3);
	float ox = (float) file->getO(1);
	float oz = (float) file->getO(3);

	// Other parameters
	int lpml = 16;
	bool fs = 0; 
	int order=6;

	// Get time variables from wav file
	int nt;
	float dt;
	float ot;
	status = file->input("Wav2d.rss");
	if(status == FILE_ERR){
		std::cout << "Error reading from wavelet file. \n";
	}
	file->close();
	nt = file->getN(1);
	dt = file->getD(1);
	ot = file->getO(1);

	// Create the classes 
	std::shared_ptr<rockseis::ModelElastic2D<float>> model (new rockseis::ModelElastic2D<float>(nx, nz, lpml, dx, dz, ox, oz, fs));
	std::shared_ptr<rockseis::WavesElastic2D<float>> waves (new rockseis::WavesElastic2D<float>(model, nt, dt, ot));
	std::shared_ptr<rockseis::Der<float>> der (new rockseis::Der<float>(nx+2*lpml, 1, nz+2*lpml, dx, 1., dz, order));
	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>(1, nt, dt));

	// Make an acoustic model
	float *R, *Vp, *Vs;
	R = model->getR();
	Vp = model->getVp();
	Vs = model->getVs();

	// Read vp model
	status = file->input("Vp2d.rss");
	file->read(Vp, nx*nz);
	file->close();

	// Read vs model
	status = file->input("Vs2d.rss");
	file->read(Vs, nx*nz);
	file->close();

	// Read rho model
	status = file->input("Rho2d.rss");
	file->read(R, nx*nz);
	file->close();

	// Stagger model
	model->staggerModels();

	// Load wavelet from file
	status = file->input("Wav2d.rss");
	source->readData();
	file->close();
	source->makeMap(model->getGeom());

	// Output snapshots to a binary file
	std::fstream myFile;
        myFile.open ("data.bin", std::ios::out | std::ios::binary);

	float *Szz;
	Szz = waves->getSzz();
	// Loop over time
	for(int it=0; it < nt; it++)
	{
		// Time stepping
		waves->forwardstepVelocity(model, der);
		waves->forwardstepStress(model, der);

		// Inserting source (Pressure)
		waves->insertSource(model, source, 2, 0, it);

		//Writting out results to binary file
		myFile.write (reinterpret_cast<char *> (Szz), (nx+2*lpml) * (nz+2*lpml) * sizeof(float));
	}	

	myFile.close();
	return 0;
}
