#ifndef SORT_H
#define SORT_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "data.h"
#include "file.h"

#define SORT_OK 1
#define SORT_ERR 0

namespace rockseis {

    typedef struct {
        size_t i0;
        size_t n;
        rs_status status;
    } key;

typedef struct
{
	    double x;
	    double offx;
	    double y;
	    double offy;
	    double z;
	    double offz;
        double foff;
	    size_t ind;

} position_t;


// =============== ABSTRACT SORT CLASS =============== //
/** The abstract sort class
 *
 */
template<typename T>
class Sort {
public:
    Sort(); ///<Constructor
    virtual ~Sort();	///< Destructor
    bool createShotmap(std::string filename) { return createSort(filename, SOURCE, 0., 0., 0); } ///< Create a shot map using a data file
    bool createShotmap(std::string filename, int sort_order) { return createSort(filename, SOURCE, 0., 0., sort_order); } ///< Create a shot map using a data file
    bool createReceivermap(std::string filename) { return createSort(filename, RECEIVER, 0., 0., 0); } ///< Create a receiver map using a data file
    bool createReceivermap(std::string filename, int sort_order) { return createSort(filename, RECEIVER, 0., 0., sort_order); } ///< Create a receiver map using a data file
    bool createCMPmap(std::string filename, T dx, T dy) { return createSort(filename, CMP, dx, dy, 0); } ///< Create a CMP map using a data file
    bool createCMPmap(std::string filename, T dx, T dy, int sort_order) { return createSort(filename, CMP, dx, dy, sort_order); } ///< Create a CMP map using a data file
    std::shared_ptr<Data2D<T>> get2DGather(); ///< Get a gather using the keymap and sortmap
    std::shared_ptr<Data2D<T>> get2DGather(size_t num); ///< Get gather number using the keymap and sortmap
    void put2DGather(std::shared_ptr<Data2D<T>>, size_t num); ///< Put 2d gather into a data file number using the keymap and sortmap
    void put2DGather(std::shared_ptr<Data2D<T>>, size_t num, Point2D<int> *mask); ///< Put 2d gather into a data file number using the keymap and sortmap and a mask
    void put3DGather(std::shared_ptr<Data3D<T>>, size_t num); ///< Put 3d gather into a data file number using the keymap and sortmap
    std::shared_ptr<Data3D<T>> get3DGather(); ///< Get a gather using the keymap and sortmap
    std::shared_ptr<Data3D<T>> get3DGather(size_t num); ///< Get gather number using the keymap and sortmap

    // Functions to read and write maps
    void readKeymap();  ///< Read a key map
    void writeKeymap(); ///< Write a key map
    void readSortmap(); ///< Read a sort map
    void writeSortmap(); ///< Write a sort map

    // Create empty dataset 
    void createEmptydataset(std::string filename, size_t n1, T d1, T o1);
    
    //Get functions
    size_t getNensemb() { return ngathers; } ///< Get number of gathers
    size_t getNtraces() { return ntraces; }  ///< Get number of traces in file
    size_t getMaxtraces();   ///< Get maximum number of traces in ensembles
    size_t getMintraces();   ///< Get minimum number of traces in ensembles
    std::string getKmapfile() { return kmapfile; } ///< Get Kmapfile
    std::string getSmapfile() { return smapfile; }  ///< Get Smapfile
    rs_status getStatus(size_t key) { if(key < ngathers && key >= 0) return keymap[key].status; else rs_error("key out of bounds."); return FAILED;} ///< Set status for a key 
    bool getReciprocity() { return reciprocity; } ///< Get reciprocity flag
    rs_key getSortkey() {  return sortkey; } ///< return sortkey

    //Set functions
    void setKmapfile(std::string name) { kmapfile = name; }
    void setSmapfile(std::string name) { smapfile = name; }
    void setDatafile(std::string name) { datafile = name; }
    void setStatus(size_t key, rs_status status); ///< Set status for a key 
    void setReciprocity(bool val) { reciprocity = val; } ///< Set reciprocity flag

private:
    size_t ngathers;
    size_t ntraces;
    key *keymap; // size nensemble, position to first trace and number of traces of ensemble, and finally a status of that ensemble 
    size_t *sortmap; // size ntrace, positions of the ensemble traces in file
    rs_key sortkey;
    std::string kmapfile;
    std::string smapfile;
    std::string datafile;
    bool reciprocity;
    bool sortset;

    bool createSort(std::string filename, rs_key _sortkey, T dx, T dy, int sort_order); ///< Create a sort map using a data file
};
}
#endif //SORT_H
