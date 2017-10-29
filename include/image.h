#ifndef IMAGE_H
#define IMAGE_H

// Include statements
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "model.h"
#include "file.h"
#include "utils.h"


namespace rockseis {
// =============== 2D IMAGE CLASS =============== //
/** The 2D image class
 *
 */
template<typename T>
class Image2D: public Model<T> {
public:
    Image2D(std::string imagefile); ///< Constructor
    Image2D(std::string imagefile, std::shared_ptr<ModelAcoustic2D<T>> model, int nhx, int nhz); 	///< Constructor
    Image2D(std::string imagefile, std::shared_ptr<ModelElastic2D<T>> model, int nhx, int nhz); 	///< Constructor
    //Image2D(std::string imagefile, ModelAcoustic2D<T> model, T apertx); 	///< Constructor
    ~Image2D();       	///< Destructor

    // Get functions
    int getNhx() { return (this->getGeom())->getN(4); }
    int getNhz() { return (this->getGeom())->getN(6); }
    T *getImagedata() { return imagedata; }  ///< Get image array
    bool getAllocated() { return allocated; } ///< Returns true if image is allocated

    // Set functions
    void setNhx(int nhx) { (this->getGeom())->setN(4, nhx); }
    void setNhz(int nhz) { (this->getGeom())->setN(6, nhz); }
    void setAllocated(bool flag) { allocated = flag; } ///< Sets allocated image flag 
    void setImagefile(std::string file) { imagefile = file; }  ///< Sets name of file

    // File input/output
    /** Read and write image files */
    bool read(); ///< Read image and coordinates
    bool write(); ///< Writes image
    void allocateImage(); /// Allocate memory for image
    void freeImage(); /// Free memory for imagedata
    bool createEmpty(); ///< Create empty image for stacking
    bool stackImage(std::string infile); ///< Stack image

private:
    T *imagedata; // Image array
    std::string imagefile;
    bool allocated;
};

// =============== 3D IMAGE CLASS =============== //
/** The 3D image class
 *
 */
template<typename T>
class Image3D: public Model<T> {
public:
    Image3D(std::string imagefile); ///< Constructor
    Image3D(std::string imagefile, std::shared_ptr<ModelAcoustic3D<T>> model, int nhx, int nhy, int nhz); 	///< Constructor
    Image3D(std::string imagefile, std::shared_ptr<ModelElastic3D<T>> model, int nhx, int nhy, int nhz); 	///< Constructor
    ~Image3D();       	///< Destructor

    // Get functions
    int getNhx() { return (this->getGeom())->getN(4); }
    int getNhy() { return (this->getGeom())->getN(5); }
    int getNhz() { return (this->getGeom())->getN(6); }
    T *getImagedata() { return imagedata; }  ///< Get image array
    bool getAllocated() { return allocated; } ///< Returns true if image is allocated

    // Set functions
    void setNhx(int nhx) { (this->getGeom())->setN(4, nhx); }
    void setNhy(int nhy) { (this->getGeom())->setN(5, nhy); }
    void setNhz(int nhz) { (this->getGeom())->setN(6, nhz); }
    void setAllocated(bool flag) { allocated = flag; } ///< Sets allocated image flag 
    void setImagefile(std::string file) { imagefile = file; }  ///< Sets name of file

    // File input/output
    /** Read and write image files */
    bool read(); ///< Read image and coordinates
    bool write(); ///< Writes image
    void allocateImage(); /// Allocate memory for image
    void freeImage(); /// Free memory for imagedata
    bool createEmpty(); ///< Create empty image for stacking
    bool stackImage(std::string infile); ///< Stack image

private:
    T *imagedata; // Image array
    std::string imagefile;
    bool allocated;
};

}

#endif //IMAGE_H
