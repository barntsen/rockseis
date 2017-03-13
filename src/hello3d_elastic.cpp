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
	// Parameters
	bool status;
	int lpml=0;
	bool fs=0;
	int order=0;

	// Parse parameters from file
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "mod2d.cfg";

    status = 0;
    try {
        cfg->parse(configFile);
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        return 1;
    }

    try {
        lpml = cfg->lookupInt(scope, "lpml");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        order = cfg->lookupInt(scope, "order");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        fs = cfg->lookupBoolean(scope, "freesurface");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    // Destroy cfg
    cfg->destroy();

    if(status == 1){
        std::cerr << "Program terminated due to input errors." << std::endl;
        return 1;
    }

    std::shared_ptr<rockseis::File> file (new rockseis::File());

	// Create the classes 
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>("Wav3d.rss"));
	std::shared_ptr<rockseis::ModelElastic3D<float>> model (new rockseis::ModelElastic3D<float>("Vp3d.rss", "Vs3d.rss", "Rho3d.rss", lpml, fs));
	std::shared_ptr<rockseis::WavesElastic3D<float>> waves (new rockseis::WavesElastic3D<float>(model, source->getNt(), source->getDt(), source->getOt(),1));
	std::shared_ptr<rockseis::Der<float>> der (new rockseis::Der<float>(model->getNx_pml(), model->getNy_pml(), model->getNz_pml(), model->getDx(), model->getDy(), model->getDz(), order));
	
	// Get models from files
	// Read model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Load wavelet from file 
	source->readData();
	source->makeMap(model->getGeom());

	// Output snapshots to a binary file
	std::shared_ptr<rockseis::File> Fsnap (new rockseis::File());
	Fsnap->output("snaps.rss");
	Fsnap->setN(1,model->getNx_pml());
	Fsnap->setD(1,model->getDx());
	Fsnap->setO(1,model->getOx());
	Fsnap->setN(2,model->getNy_pml());
	Fsnap->setD(2,model->getDy());
	Fsnap->setO(2,model->getOy());
	Fsnap->setN(3,model->getNz_pml());
	Fsnap->setD(3,model->getDz());
	Fsnap->setO(3,model->getOz());
	Fsnap->setN(4,waves->getNt());
	Fsnap->setD(4,waves->getDt());
	Fsnap->setO(4,waves->getOt());
	Fsnap->setData_format(sizeof(float));
	Fsnap->writeHeader();
	Fsnap->seekp(Fsnap->getStartofdata());

	// Pointer to one of the Fields
	float *Sxx;
	Sxx = waves->getSxx();

	// Loop over time
	for(int it=0; it < waves->getNt(); it++)
	{
		// Time stepping
		waves->forwardstepVelocity(model, der);
		waves->forwardstepStress(model, der);

		// Inserting source
		waves->insertSource(model, source, 3, 0, it);

		//Writting out results to binary file
		Fsnap->write(Sxx, model->getNx_pml() * model->getNy_pml() * model->getNz_pml());
	}	

	Fsnap->close();

	return 0;
}
