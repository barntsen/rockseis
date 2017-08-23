#ifndef PML_H
#define PML_H

// Include statements
#include <vector>
#include <string>
#include <stdlib.h>
#include <math.h>

namespace rockseis {
// =============== ABSTRACT PML CLASS =============== //
/** The abstract pml class
 *
 */
template<typename T>
class Pml {
public:
    Pml();		///< Constructor
    Pml(const int _L, const T _dt);	///< Constructor with dimension
    virtual ~Pml();	///< Destructor
    
    // Get functions
    T getAmax() { return Amax; }	///< Get Amax
    T getSmax() { return Smax; }   ///< Get Smax
    T getKmax() { return Kmax; }	///< Get Kmax
    int getLpml() { return Lpml; }	///< Get Lpml
    T getDt() { return dt; }	///< Get dt
    
    // Set functions
    void setAmax(T _Amax) { Amax= _Amax;}	///< Set Amax
    void setKmax(T _Kmax) { Kmax= _Kmax;}	///< Set Kmax
    void setSmax(T _Smax) { Smax= _Smax;}	///< Set Smax
    void setDt(T _dt) { dt= _dt;}	///< Set dt
    void setLpml(int _L) { Lpml = _L;}	///< Set Lpml
    
    /** Compute A,B and C constants. 
     * Uses Amax, Kmax and Smax to compute the PML variables.
    * */
    void computeABC();
    
    // Left, top and front constants 
    T *A_ltf; // Non-staggered
    T *B_ltf; // Non-staggered
    T *C_ltf; // Non-staggered
    T *A_ltf_stag; // Staggered 
    T *B_ltf_stag; // Staggered
    T *C_ltf_stag; // Staggered
    // Right, bottom and back constants 
    T *A_rbb; // Non-staggered
    T *B_rbb; // Non-staggered
    T *C_rbb; // Non-staggered
    T *A_rbb_stag; // Staggered
    T *B_rbb_stag; // Staggered
    T *C_rbb_stag; // Staggered
private:
    int Lpml; // Length of PML boundary
    T dt;     // Time sampling interval
    T Amax; // PML constant (usually = pi*f0) where f0 is the dominant frequency
    T Smax; // Pml constant (usually = f0*3*vpmax/(10*lpml))
    T Kmax; // Pml constant (usually = 2)
};

// =============== 1D ACOUSTIC PML CLASS =============== //
/** The 1D Acoustic PML model class
 *
 */
template<typename T>
class PmlAcoustic1D: public Pml<T> {
public:
    PmlAcoustic1D();	///< Constructor
    PmlAcoustic1D(const int Lpml, const T dt);	///< Constructor
    void callcompABC() { this->computeABC(); }  ///< Interface to computeABC()
    ~PmlAcoustic1D();	///< Destructor
    
    T *P_top;
    T *P_bottom;
    T *Azz_top;
    T *Azz_bottom;
};

// =============== 2D ACOUSTIC PML CLASS =============== //
/** The 2D Acoustic PML model class
 *
 */
template<typename T>
class PmlAcoustic2D: public Pml<T> {
public:
    PmlAcoustic2D();	///< Constructor
    PmlAcoustic2D(const int nx, const int nz, const int Lpml, const T dt);	///< Constructor
    void callcompABC() { this->computeABC(); }  ///< Interface to computeABC()
    ~PmlAcoustic2D();	///< Destructor
    
    T *P_left;
    T *P_right;
    T *P_top;
    T *P_bottom;
    T *Axx_left;
    T *Axx_right;
    T *Azz_top;
    T *Azz_bottom;
};

// =============== 3D ACOUSTIC PML CLASS =============== //
/** The 3D Acoustic PML model class
 *
 */
template<typename T>
class PmlAcoustic3D: public Pml<T> {
public:
    PmlAcoustic3D();	///< Constructor
    PmlAcoustic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt);	///< Constructor
    void callcompABC() { this->computeABC(); }  ///< Interface to computeABC()
    ~PmlAcoustic3D();	///< Destructor
    
    T *P_left;
    T *P_right;
    T *P_top;
    T *P_bottom;
    T *P_front;
    T *P_back;
    T *Axx_left;
    T *Axx_right;
    T *Ayy_front;
    T *Ayy_back;
    T *Azz_top;
    T *Azz_bottom;
};

// =============== 2D ELASTIC PML CLASS =============== //
/** The 2D Elastic PML model class
 *
 */
template<typename T>
class PmlElastic2D: public Pml<T> {
public:
    PmlElastic2D();	///< Constructor
    PmlElastic2D(const int nx, const int nz, const int Lpml, const T dt);	///< Constructor
    ~PmlElastic2D();	///< Destructor
    
    T *Sxx_left;
    T *Sxx_right;
    T *Szz_top;
    T *Szz_bottom;
    
    T *Sxzz_top;
    T *Sxzz_bottom;
    T *Sxzx_left;
    T *Sxzx_right;
    
    T *Vxx_left;
    T *Vxx_right;
    T *Vzx_left;
    T *Vzx_right;
    
    T *Vzz_top;
    T *Vzz_bottom;
    T *Vxz_top;
    T *Vxz_bottom;
    
};

// =============== 3D ELASTIC PML CLASS =============== //
/** The 3D Elastic PML model class
 *
 */
template<typename T>
class PmlElastic3D: public Pml<T> {
public:
    PmlElastic3D();	///< Constructor
    PmlElastic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt); 	///< Constructor
    ~PmlElastic3D();	///< Destructor
    
    T *Sxx_left;
    T *Sxx_right;
    T *Syy_front;
    T *Syy_back;
    T *Szz_top;
    T *Szz_bottom;
    
    T *Sxzz_top;
    T *Sxzz_bottom;
    T *Syzz_top;
    T *Syzz_bottom;
    
    T *Sxzx_left;
    T *Sxzx_right;
    T *Sxyx_left;
    T *Sxyx_right;
    
    T *Sxyy_front;
    T *Sxyy_back;
    T *Syzy_front;
    T *Syzy_back;
    
    T *Vxx_left;
    T *Vxx_right;
    T *Vzx_left;
    T *Vzx_right;
    T *Vyx_left;
    T *Vyx_right;
    
    T *Vzz_top;
    T *Vzz_bottom;
    T *Vxz_top;
    T *Vxz_bottom;
    T *Vyz_top;
    T *Vyz_bottom;
    
    T *Vxy_front;
    T *Vxy_back;
    T *Vzy_front;
    T *Vzy_back;
    T *Vyy_front;
    T *Vyy_back;
};

}

#endif //PML_H
