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
    std::string Vsfile;
    std::string Rhofile;
    bool Psnap=0, Precord=0;
    std::string Psnapfile;
    std::string Precordfile;
    std::shared_ptr<rockseis::Data3D<float>> Pdata;

    bool Vxsnap=0, Vxrecord=0;
    std::string Vxsnapfile;
    std::string Vxrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vxdata;

    bool Vysnap=0, Vyrecord=0;
    std::string Vysnapfile;
    std::string Vyrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vydata;

    bool Vzsnap=0, Vzrecord=0;
    std::string Vzsnapfile;
    std::string Vzrecordfile;
    std::shared_ptr<rockseis::Data3D<float>> Vzdata;

	// Parse parameters from file
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "elamod3d.cfg";

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
        Vsfile = cfg->lookupString(scope, "Vs");
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
        Vxsnap = cfg->lookupBoolean(scope, "Vxsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vxsnap){
        try {
            Vxsnapfile = cfg->lookupString(scope, "Vxsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Vxrecord = cfg->lookupBoolean(scope, "Vxrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vxrecord){
        try {
            Vxrecordfile = cfg->lookupString(scope, "Vxrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
   try {
        Vysnap = cfg->lookupBoolean(scope, "Vysnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vysnap){
        try {
            Vysnapfile = cfg->lookupString(scope, "Vysnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Vyrecord = cfg->lookupBoolean(scope, "Vyrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vyrecord){
        try {
            Vyrecordfile = cfg->lookupString(scope, "Vyrecordfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
   try {
        Vzsnap = cfg->lookupBoolean(scope, "Vzsnap");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vzsnap){
        try {
            Vzsnapfile = cfg->lookupString(scope, "Vzsnapfile");
        } catch(const config4cpp::ConfigurationException & ex) {
            std::cerr << ex.c_str() << std::endl;
            status = 1;
        }
    }
    try {
        Vzrecord = cfg->lookupBoolean(scope, "Vzrecord");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
    if(Vzrecord){
        try {
            Vzrecordfile = cfg->lookupString(scope, "Vzrecordfile");
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
	std::shared_ptr<rockseis::ModelElastic3D<float>> model (new rockseis::ModelElastic3D<float>(Vpfile, Vsfile, Rhofile, lpml ,fs));
	std::shared_ptr<rockseis::Data3D<float>> source (new rockseis::Data3D<float>(Sourcefile));
	std::shared_ptr<rockseis::ModellingElastic3D<float>> modelling (new rockseis::ModellingElastic3D<float>(model, source, order, snapinc));

    // Setting Snapshot file 
    if(Psnap){
        modelling->setSnapP(Psnapfile);
    }
    if(Vxsnap){
        modelling->setSnapVx(Vxsnapfile);
    }
    if(Vysnap){
        modelling->setSnapVy(Vysnapfile);
    }
    if(Vzsnap){
        modelling->setSnapVz(Vzsnapfile);
    }

    //Setting sourcetype to VX (x-Force) type
    source->setField(rockseis::VY);

    // Setting Record
    if(Precord){
        Pdata = std::make_shared<rockseis::Data3D<float>>(Precordfile, source->getNt(), source->getDt(), 0.0);
        Pdata->setField(rockseis::PRESSURE);
        // Load data geometry from file
        Pdata->readCoords();
        Pdata->makeMap(model->getGeom());
        modelling->setRecP(Pdata);
    }
    // Setting Record
    if(Vxrecord){
        Vxdata = std::make_shared<rockseis::Data3D<float>>(Vxrecordfile, source->getNt(), source->getDt(), 0.0);
        Vxdata->setField(rockseis::VX);
        // Load data geometry from file
        Vxdata->readCoords();
        Vxdata->makeMap(model->getGeom());
        modelling->setRecVx(Vxdata);
    }
    if(Vyrecord){
        Vydata = std::make_shared<rockseis::Data3D<float>>(Vyrecordfile, source->getNt(), source->getDt(), 0.0);
        Vydata->setField(rockseis::VY);
        // Load data geometry from file
        Vydata->readCoords();
        Vydata->makeMap(model->getGeom());
        modelling->setRecVy(Vydata);
    }
    // Setting Record
    if(Vzrecord){
        Vzdata = std::make_shared<rockseis::Data3D<float>>(Vzrecordfile, source->getNt(), source->getDt(), 0.0);
        Vzdata->setField(rockseis::VZ);
        // Load data geometry from file
        Vzdata->readCoords();
        Vzdata->makeMap(model->getGeom());
        modelling->setRecVz(Vzdata);
    }

	// Read acoustic model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Read wavelet data and coordinates and make a map
	source->read();
	source->makeMap(model->getGeom());

	// Run modelling 
	modelling->run();

	// Output records
    if(Precord){
        Pdata->write();
    }

    if(Vxrecord){
        Vxdata->write();
    }

    if(Vyrecord){
        Vydata->write();
    }

    if(Vzrecord){
        Vzdata->write();
    }

	return 0;
}
