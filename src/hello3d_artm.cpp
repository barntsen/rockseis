#include <iostream>
#include "model.h"
#include "rtm.h"
#include "pml.h"
#include "waves.h"
#include "utils.h"
#include "der.h"
#include "data.h"
#include "file.h"
#include "image.h"
#include <memory>
#include <fstream>
#include <math.h>
#include <config4cpp/Configuration.h>

int main()
{

	bool status;

	// Parameters
	int lpml=0;
    int nhx, nhy, nhz;
	bool fs=0;
	int order=0;
	int snapinc=0;
    std::string Sourcefile;
    std::string Precordfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Psnapfile;
    std::string Pimagefile;
    std::shared_ptr<rockseis::Data3D<float>> Pdata;

	// Parse parameters from file
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "acurtm3d.cfg";

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
        snapinc = cfg->lookupInt(scope, "snapinc");
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
    try {
        Vpfile = cfg->lookupString(scope, "Vp");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        Rhofile = cfg->lookupString(scope, "Rho");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Sourcefile = cfg->lookupString(scope, "source");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Psnapfile = cfg->lookupString(scope, "Psnapfile");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Pimagefile = cfg->lookupString(scope, "Pimagefile");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    try {
        Precordfile = cfg->lookupString(scope, "Precordfile");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        nhx = cfg->lookupInt(scope, "nhx");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        nhy = cfg->lookupInt(scope, "nhy");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        nhz = cfg->lookupInt(scope, "nhz");
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
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> model (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));
    std::shared_ptr<rockseis::Image3D<float>> pimage (new rockseis::Image3D<float>(Pimagefile, model, nhx, nhy, nhz));
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Sourcefile));
	std::shared_ptr<rockseis::Data3D<float>> data (new rockseis::Data3D<float>(Precordfile));
	std::shared_ptr<rockseis::RtmAcoustic3D<float>> rtm (new rockseis::RtmAcoustic3D<float>(model, pimage, source, data, order, snapinc));

    // Setting Snapshot file 
    rtm->setSnapfile(Psnapfile);

    rtm->setNcheck(11);
    rtm->setIncore(true);

   	// Read acoustic model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Read wavelet data and coordinates and make a map
	source->read();
	source->makeMap(model->getGeom());

    //Read input data 
	data->read();
	data->makeMap(model->getGeom());
    //rockseis::rs_snapmethod snapmethod = rockseis::OPTIMAL;
    rockseis::rs_snapmethod snapmethod = rockseis::FULL;
    //rockseis::rs_snapmethod snapmethod = rockseis::EDGES;

	// Run modelling 
    switch(snapmethod){
        case rockseis::FULL:
            rtm->run();
            break;
        case rockseis::OPTIMAL:
            rtm->run_optimal();
            break;
        default:
            rockseis::rs_error("Invalid option of snapshot saving."); 
    }

    // Write out image file
    pimage->write();

	return 0;
}
