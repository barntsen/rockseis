#ifndef MODELLING_H
#define MODELLING_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "geometry.h"
#include "utils.h"
#include "file.h"
#include "model.h"
#include "data.h"
#include "waves.h"
#include "der.h"

#define MOD_OK 1
#define MOD_ERR 0

namespace rockseis {
// =============== ABSTRACT MODELLING CLASS =============== //
/** The abstract modelling class
 *
 */

template<typename T>
class Modelling {
public:
    Modelling();	///< Constructor
    Modelling(int order);  ///< Constructor
    
    // Modelling functions
    int Acoustic2D(std::shared_ptr<ModelAcoustic2D<T>> model, std::shared_ptr<Data2D<T>> source, std::shared_ptr<Data2D<T>> recP); ///< Runs modelling with parameters (model, data and parameters)
    int getOrder() { return order; } ///< Get order of FD stencil
    void setOrder(int _order) { if(_order > 1 && _order < 9)  order = _order;} ///< Set order of FD stencil

    ~Modelling();	///< Constructor


private:
    int order;  ///< Order of the FD stencil 
};


}
#endif //MODELLING_H
