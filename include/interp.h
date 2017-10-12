#ifndef INTERP_H
#define INTERP_H

// Include statements
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <utils.h>
#include <data.h>

#define PI 3.14159265358979323846
#define LN_PAD 1
#define BS_PAD 2
#define SN_PAD 3

namespace rockseis {

// =============== INTERP CLASS =============== //
/** The Interp class.
 * Interp for seeking memory in 2D and 3D arrays.
 * */
template<typename T>
class Interp{
public:
    Interp(); ///< Constructor
    Interp(rs_interpmode _mode); ///< Constructor
    ~Interp(); ///< Destructor
    T interp_bs_int1d(T *mod, T d1);
    T interp_ln_int1d(T *mod, T d1);
    T interp_sinc_int1d(T *mod, T d1);
    int interp_seekindex(const T *vector, const int n, const T value);
    void interp(std::shared_ptr<Data2D<T>> from, std::shared_ptr<Data2D<T>> to);
    void interp(std::shared_ptr<Data3D<T>> from, std::shared_ptr<Data3D<T>> to);

private:
    T sinc(const T x);
    rs_interpmode mode;
    int pad;
    T wrk[6];
};


// ================ Mathematical functions ============//
}
#endif //INTERP_H
