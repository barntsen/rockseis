#include <iostream>
#include <memory>
#include <fstream>
#include <math.h>
#include "data.h"
#include "file.h"
#include "sort.h"
#include "utils.h"
#include "geometry.h"

int main()
{
    std::string filename;
    filename = "Shots.rss";
	std::shared_ptr<rockseis::Data3D<float>> Shots (new rockseis::Data3D<float>(filename));
	std::shared_ptr<rockseis::Data3D<float>> OneShot;
	std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->createShotmap(filename); 
    std::cerr << "Number of shot gathers: " << Sort->getNensemb() << std::endl;
    rockseis::Point3D<float> *scoords;
    while((OneShot = Sort->get3DGather()) != nullptr ){
        scoords = (OneShot->getGeom())->getScoords();
        std::cerr << "Got shot with coordinates x: " << scoords[0].x << " y: " << scoords[0].y << " z: " << scoords[0].z << std::endl;
        OneShot.reset();
    }

    Sort->createReceivermap(filename); 
    std::cerr << "Number of receiver gathers: " << Sort->getNensemb() << std::endl;
    rockseis::Point3D<float> *gcoords;
    while((OneShot = Sort->get3DGather()) != nullptr ){
        gcoords = (OneShot->getGeom())->getGcoords();
        std::cerr << "Got receiver gather with coordinates x: " << gcoords[0].x << " y: " << gcoords[0].y << " z: " << gcoords[0].z << std::endl;
        OneShot.reset();
    }

    Sort->createCMPmap(filename, 50.0, 50.0); 
    std::cerr << "Number of CMP gathers: " << Sort->getNensemb() << std::endl;
    float cmpx, cmpy;
    while((OneShot = Sort->get3DGather()) != nullptr ){
        scoords = (OneShot->getGeom())->getScoords();
        gcoords = (OneShot->getGeom())->getGcoords();
        cmpx = 0.5*(scoords[0].x + gcoords[0].x);
        cmpy = 0.5*(scoords[0].y + gcoords[0].y);
        std::cerr << "Got CMP with coordinates x: " << cmpx << " y: " << cmpy << std::endl;
        OneShot.reset();
    }

    Sort->writeKeymap();
    Sort->writeSortmap();

	std::shared_ptr<rockseis::Sort<float>> Sort2 (new rockseis::Sort<float>());
    Sort2->readKeymap();
    Sort2->readSortmap();

    Sort2->setKmapfile("kmap2.rss");
    Sort2->setSmapfile("smap2.rss");
    Sort2->writeKeymap();
    Sort2->writeSortmap();
	return 0;
}
