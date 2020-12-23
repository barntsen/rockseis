#include "domain.h"

namespace rockseis {
// constructor
template<typename T>
Domain<T>::Domain()
{
    geometry = std::make_shared<Geometry<T>>(); 
    lpml = 0;
    low = -1;
    high = -1;
    padl = 0;
    padh = 0;
    nd = 1;
    d = 0;
}

template<typename T>
Domain<T>::~Domain() {
   // Do nothing
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Domain<float>;
template class Domain<double>;

}
