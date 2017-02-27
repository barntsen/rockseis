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

    std::shared_ptr<rockseis::File> file (new rockseis::File());

//    // Get time variables from wav file
//    int nt;
//    float dt;
//    float ot;
//    int nsou;
//    status = file->input("Wav2d.rss");
//    if(status == FILE_ERR){
//	    std::cout << "Error reading from wavelet file. \n";
//	    exit(1);
//    }
//    file->close();
//    nsou = file->getN(2);
//    nt = file->getN(1);
//    dt = file->getD(1);
//    ot = file->getO(1);
//

    //Get information from data geometry file
	status = file->input("Geom2d.rss");
    if(status == FILE_ERR){
		std::cout << "Error reading from geometry file. \n";
		exit(1);
	}
    int nrec = file->getN(2);
    file->close();

	// Create the classes 
	std::shared_ptr<rockseis::ModelAcoustic2D<float>> model (new rockseis::ModelAcoustic2D<float>("Vp2d.rss", "Rho2d.rss", lpml ,fs));
	std::shared_ptr<rockseis::Modelling<float>> modelling (new rockseis::Modelling<float>(order));
	std::shared_ptr<rockseis::Data2D<float>> source (new rockseis::Data2D<float>("Wav2d.rss"));
	std::shared_ptr<rockseis::Data2D<float>> record (new rockseis::Data2D<float>(nrec, source->getNt(), source->getDt()));

	// Read acoustic model
	model->readModel();

	// Stagger model
	model->staggerModels();

	// Load wavelet from file
//	status = file->input("Wav2d.rss");
//	if(status == FILE_ERR){
//		std::cout << "Error reading from Wav file. \n";
//		exit(1);
//	}
	source->readData();
//	file->close();
	source->makeMap(model->getGeom());

	// Load data geometry from file
	record->setFile("Geom2d.rss");
        record->readCoords();
	record->makeMap(model->getGeom());
	//rockseis::Point2D<float> *scoords2d = (record->getGeom())->getScoords();
	//rockseis::Point2D<float> *gcoords2d = (record->getGeom())->getGcoords();
    //DEBUG, print receiver coordinates
//    for (int i=0; i< nrec; i++){
  //      std::cout << "sx[" << i <<"]=" << scoords2d[i].x << "\n";  
    //    std::cout << "sy[" << i <<"]=" << scoords2d[i].y << "\n";  
      //  std::cout << "gx[" << i <<"]=" << gcoords2d[i].x << "\n";  
     //   std::cout << "gy[" << i <<"]=" << gcoords2d[i].y << "\n";  
   // }

	// Run modelling 
	modelling->Acoustic2D(model, source, record);

	return 0;
}
