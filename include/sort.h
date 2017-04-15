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
	    double y;
	    double z;
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
    void setStatus(size_t key, rs_status status); ///< Set status for a key 
    rs_status getStatus(size_t key) { if(key < nensembles && key >= 0) return keymap[key].status; else rs_error("key out of bounds."); return FAILED;} ///< Set status for a key 
    
    bool createShotmap(std::string filename) { return createSort(filename, SOURCE, 0., 0.); } ///< Create a shot map using a data file
    bool createReceivermap(std::string filename) { return createSort(filename, RECEIVER, 0., 0.); } ///< Create a receiver map using a data file
    bool createCMPmap(std::string filename, T dx, T dy) { return createSort(filename, CMP, dx, dy); } ///< Create a CMP map using a data file
    bool getGather(std::shared_ptr<Data2D<T>> gather);

    std::shared_ptr<Data3D<T>> getGather();

    // Read and write maps functions
    ///< Read a keymap to file
    ///< Write a keymap to file
    ///< Read a sortmap to file
    ///< Write a sortmap to file
    
    //Get functions
    size_t getNensemb() { return nensembles; } 

private:
    size_t nensembles;
    size_t ntraces;
    key *keymap; // size nensemble, position to first trace and number of traces of ensemble, and finally a status of that ensemble 
    size_t *sortmap; // size ntrace, positions of the ensemble traces in file
    rs_key sortkey;
    std::string kmapfile;
    std::string smapfile;
    std::string datafile;

    bool createSort(std::string filename, rs_key _sortkey, T dx, T dy); ///< Create a sort map using a data file
};
}
#endif //SORT_H
