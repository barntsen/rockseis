#ifndef SORT_H
#define SORT_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "file.h"


namespace rockseis {
// =============== ABSTRACT SORT CLASS =============== //
/** The abstract sort class
 *
 */
template<typename T>
class Sort {
public:
    Sort(); ///<Constructor
    virtual ~Sort();	///< Destructor

    // Create a sort map using a data file

private:
    int nensembles;
    int ntraces;
    int *keymap; // size 2*nensemble, position to first trace and number of traces of ensemble 
    int *sortmap; // size ntrace, positions of the ensemble traces in file
    rs_key sortkey;
}
#endif //SORT_H
