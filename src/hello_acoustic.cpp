#include <iostream>
#include "model.h"
#include "modelling.h"
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

	// Parameters
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

	// Create the classes 
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> model (new rockseis::ModelAcoustic2D<float>("Vp2d.rss", "Rho2d.rss", lpml ,fs));
	std::shared_ptr<rockseis::Modelling<float>> modelling (new rockseis::Modelling<float>(order));
	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>("Wav2d.rss"));
	std::shared_ptr<rockseis::Data2D<float>> record (new rockseis::Data2D<float>("Geom2d.rss", source->getNt(), source->getDt()));

	// Read acoustic model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Read wavelet data and coordinates and make a map
	source->readData();
	source->makeMap(model->getGeom());

	// Load recording geometry from file
	record->readCoords();
	record->makeMap(model->getGeom());

	// Run modelling 
	modelling->Acoustic2D(model, source, record);

	// Output record
	record->setFile("Shot.rss");
	record->write();

	return 0;
}
