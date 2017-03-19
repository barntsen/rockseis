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
	int snapinc=0;
    std::string Sourcefile;
    std::string Vpfile;
    std::string Rhofile;
    bool Psnap=0, Precord=0;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data3D<float>> Pdata;

    bool Axsnap=0, Axrecord=0;
    std::string Axsnapfile;
    std::string Axrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Axdata;

    bool Aysnap=0, Ayrecord=0;
    std::string Aysnapfile;
    std::string Ayrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Aydata;

    bool Azsnap=0, Azrecord=0;
    std::string Azsnapfile;
    std::string Azrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Azdata;

	// Parse parameters from file
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "mod3d.cfg";

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
        Psnap = cfg->lookupBoolean(scope, "Psnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Psnap){
        try {
            Psnapfile = cfg->lookupString(scope, "Psnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Precord = cfg->lookupBoolean(scope, "Precord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Precord){
        try {
            Precordfile = cfg->lookupString(scope, "Precordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
   try {
        Axsnap = cfg->lookupBoolean(scope, "Axsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Axsnap){
        try {
            Axsnapfile = cfg->lookupString(scope, "Axsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Axrecord = cfg->lookupBoolean(scope, "Axrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Axrecord){
        try {
            Axrecordfile = cfg->lookupString(scope, "Axrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
   }
   try {
        Aysnap = cfg->lookupBoolean(scope, "Aysnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Aysnap){
        try {
            Aysnapfile = cfg->lookupString(scope, "Aysnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Ayrecord = cfg->lookupBoolean(scope, "Ayrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Ayrecord){
        try {
            Ayrecordfile = cfg->lookupString(scope, "Ayrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }

   try {
        Azsnap = cfg->lookupBoolean(scope, "Azsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Azsnap){
        try {
            Azsnapfile = cfg->lookupString(scope, "Azsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Azrecord = cfg->lookupBoolean(scope, "Azrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Azrecord){
        try {
            Azrecordfile = cfg->lookupString(scope, "Azrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }

	// Destroy cfg
	cfg->destroy();

	if(status == 1){
		std::cerr << "Program terminated due to input errors." << std::endl;
		return 1;
	}

	// Create the classes 
	std::shared_ptr<rockseis::ModelAcoustic3D<float>> model (new rockseis::ModelAcoustic3D<float>(Vpfile, Rhofile, lpml ,fs));
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Sourcefile));
	std::shared_ptr<rockseis::ModellingAcoustic3D<float>> modelling (new rockseis::ModellingAcoustic3D<float>(model, source, order, snapinc));

    // Setting Snapshot file 
    if(Psnap){
        modelling->setSnapP(Psnapfile);
    }
    if(Axsnap){
        modelling->setSnapAx(Axsnapfile);
    }
    if(Aysnap){
        modelling->setSnapAy(Aysnapfile);
    }
    if(Azsnap){
        modelling->setSnapAz(Azsnapfile);
    }

    // Setting Record
    if(Precord){
        Pdata = std::make_shared<rockseis::Data3D<float>>(Precordfile, source->getNt(), source->getDt());
        // Load data geometry from file
        Pdata->readCoords();
        Pdata->makeMap(model->getGeom());
        modelling->setRecP(Pdata);
    }
    // Setting Record
    if(Axrecord){
        Axdata = std::make_shared<rockseis::Data3D<float>>(Axrecordfile, source->getNt(), source->getDt());
        // Load data geometry from file
        Axdata->readCoords();
        Axdata->makeMap(model->getGeom());
        modelling->setRecAx(Axdata);
    }
    // Setting Record
    if(Ayrecord){
        Aydata = std::make_shared<rockseis::Data3D<float>>(Ayrecordfile, source->getNt(), source->getDt());
        // Load data geometry from file
        Aydata->readCoords();
        Aydata->makeMap(model->getGeom());
        modelling->setRecAy(Aydata);
    }
    // Setting Record
    if(Azrecord){
        Azdata = std::make_shared<rockseis::Data3D<float>>(Azrecordfile, source->getNt(), source->getDt());
        // Load data geometry from file
        Azdata->readCoords();
        Azdata->makeMap(model->getGeom());
        modelling->setRecAz(Azdata);
    }

	// Read acoustic model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Read wavelet data and coordinates and make a map
	source->readData();
	source->makeMap(model->getGeom());

	// Run modelling 
	modelling->run();

	// Output record
    if(Precord){
        Pdata->write();
    }

    if(Axrecord){
        Axdata->write();
    }
    if(Ayrecord){
        Aydata->write();
    }
    if(Azrecord){
        Azdata->write();
    }


	return 0;
}
