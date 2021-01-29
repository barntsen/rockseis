// Include statements
#include "pml.h"

namespace rockseis {

// =============== ABSTRACT PML CLASS =============== //
template<typename T>
Pml<T>::Pml() {
    Amax = AMAX;
    Smax = 1200.;
    Kmax = KMAX;
    dt = 0.0;
    Lpml=10;
    A_ltf = (T *) malloc(Lpml*sizeof(T));
    B_ltf = (T *) malloc(Lpml*sizeof(T));
    C_ltf = (T *) malloc(Lpml*sizeof(T));
    A_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    B_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    C_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    A_rbb = (T *) malloc(Lpml*sizeof(T));
    B_rbb = (T *) malloc(Lpml*sizeof(T));
    C_rbb = (T *) malloc(Lpml*sizeof(T));
    A_rbb_stag = (T *) malloc(Lpml*sizeof(T));
    B_rbb_stag = (T *) malloc(Lpml*sizeof(T));
    C_rbb_stag = (T *) malloc(Lpml*sizeof(T));
 
    for(int i=0; i<6; i++) setApplypml(i,true);

    computeABC();

}

template<typename T>
Pml<T>::~Pml() {
    /* Free variables */
    free(A_ltf);
    free(B_ltf);
    free(C_ltf);
    free(A_ltf_stag);
    free(B_ltf_stag);
    free(C_ltf_stag);
    free(A_rbb);
    free(B_rbb);
    free(C_rbb);
    free(A_rbb_stag);
    free(B_rbb_stag);
    free(C_rbb_stag);

}

template<typename T>
Pml<T>::Pml(const int _Lpml, const T _dt) {
    Amax = AMAX;
    Smax = 1200.;
    Kmax = KMAX;
    dt = _dt;
    Lpml = _Lpml;
    A_ltf = (T *) malloc(Lpml*sizeof(T));
    B_ltf = (T *) malloc(Lpml*sizeof(T));
    C_ltf = (T *) malloc(Lpml*sizeof(T));
    A_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    B_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    C_ltf_stag = (T *) malloc(Lpml*sizeof(T));
    A_rbb = (T *) malloc(Lpml*sizeof(T));
    B_rbb = (T *) malloc(Lpml*sizeof(T));
    C_rbb = (T *) malloc(Lpml*sizeof(T));
    A_rbb_stag = (T *) malloc(Lpml*sizeof(T));
    B_rbb_stag = (T *) malloc(Lpml*sizeof(T));
    C_rbb_stag = (T *) malloc(Lpml*sizeof(T));

    for(int i=0; i<6; i++) setApplypml(i,true);
    
    computeABC();

}


template<typename T>
void Pml<T>::computeABC()
{
    T func, a, k, s;
    int i;
    for(i=0; i<Lpml; i++){
	/* Left and Top*/
        /* Non Staggered */
        func = (T) (Lpml - i)/Lpml;
        func = pow(func,3); 
        a = Amax*func;
        k = 1. + (Kmax-1)*func;
        s = Smax*func;
        C_ltf[i] = 1.0 - 1.0/k;
        B_ltf[i] = exp(-(a + s/k)*dt);
        A_ltf[i] = (s*(1.-B_ltf[i]))/(k*(k*a + s));

        /* Staggered */
        func = (T) (Lpml - i - 0.5)/Lpml;
        func = pow(func,3); 
        a = Amax*func;
        k = 1. + (Kmax-1)*func;
        s = Smax*func;
        C_ltf_stag[i] = 1.0 - 1.0/k;
        B_ltf_stag[i] = exp(-(a + s/k)*dt);
        A_ltf_stag[i] = (s*(1.-B_ltf_stag[i]))/(k*(k*a + s));

	/* Right and Bottom */
        /* Non Staggered */
        func = (T) (i + 1.0)/Lpml;
        func = pow(func,3); 
        a = Amax*func;
        k = 1. + (Kmax-1)*func;
        s = Smax*func;
        C_rbb[i] = 1.0 - 1.0/k;
        B_rbb[i] = exp(-(a + s/k)*dt);
        A_rbb[i] = (s*(1.-B_rbb[i]))/(k*(k*a + s));

        /* Staggered */
        func = (T) (i + 1.0 + 0.5)/Lpml;
        func = pow(func,3); 
        a = Amax*func;
        k = 1. + (Kmax-1)*func;
        s = Smax*func;
        C_rbb_stag[i] = 1.0 - 1.0/k;
        B_rbb_stag[i] = exp(-(a + s/k)*dt);
        A_rbb_stag[i] = (s*(1.-B_rbb_stag[i]))/(k*(k*a + s));
    }
}


// =============== 1D ACOUSTIC PML CLASS =============== //
template<typename T>
PmlAcoustic1D<T>::PmlAcoustic1D(): Pml<T>() {
    /* Default constructor, not to be used */ 
    /* Minimal allocation to avoid free from crashing */
    P_top=(T *) malloc(1);
    P_bottom=(T *) malloc(1);
    Azz_top=(T *) malloc(1);
    Azz_bottom=(T *) malloc(1);

}

template<typename T>
PmlAcoustic1D<T>::PmlAcoustic1D(const int Lpml, const T dt): Pml<T>(Lpml, dt) {

    /* Allocate variables */
    P_top=(T *) calloc(Lpml,sizeof(T));
    P_bottom=(T *) calloc(Lpml,sizeof(T));
    Azz_top=(T *) calloc(Lpml,sizeof(T));
    Azz_bottom=(T *) calloc(Lpml,sizeof(T));
}

template<typename T>
PmlAcoustic1D<T>::~PmlAcoustic1D() {
    /* Free variables */
    free(P_top);
    free(P_bottom);
    free(Azz_top);
    free(Azz_bottom);
}



// =============== 2D ACOUSTIC PML CLASS =============== //
template<typename T>
PmlAcoustic2D<T>::PmlAcoustic2D(): Pml<T>() {
    /* Default constructor, not to be used */ 
    /* Minimal allocation to avoid free from crashing */
    P_left=(T *) malloc(1);
    P_right=(T *) malloc(1);
    P_top=(T *) malloc(1);
    P_bottom=(T *) malloc(1);
    Axx_left=(T *) malloc(1);
    Axx_right=(T *) malloc(1);
    Azz_top=(T *) malloc(1);
    Azz_bottom=(T *) malloc(1);

}

template<typename T>
PmlAcoustic2D<T>::PmlAcoustic2D(const int nx, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    P_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    P_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    P_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    P_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Axx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Axx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Azz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Azz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
}

template<typename T>
PmlAcoustic2D<T>::PmlAcoustic2D(const int nx, const int nz, const int Lpml, const T dt, const int dim, const bool low, const bool high): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    switch (dim){
       case 0:
          if(low){
             P_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
             Axx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
             this->setApplypml(0,true);
          }
          if(high){
             P_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
             Axx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
             this->setApplypml(1,true);
          }
          P_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
          P_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
          Azz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
          Azz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
          this->setApplypml(4,true);
          this->setApplypml(5,true);
          break;
       case 2:
          if(low){
             P_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
             Azz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
             this->setApplypml(4,true);
          }
          if(high){
             P_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
             Azz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
             this->setApplypml(5,true);
          }
          P_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
          P_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
          Axx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
          Axx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
          this->setApplypml(0,true);
          this->setApplypml(1,true);
          break;
       default:
          rs_error("PmlAcoustic2D::Invalid value for dim");
          break;
    }
}

template<typename T>
PmlAcoustic2D<T>::~PmlAcoustic2D() {
   /* Free variables */
   if(this->getApplypml(0)){
      free(P_left);
      free(Axx_left);
   }
   if(this->getApplypml(1)){
      free(P_right);
      free(Axx_right);
   }
   if(this->getApplypml(4)){
      free(P_top);
      free(Azz_top);
   }
   if(this->getApplypml(5)){
      free(P_bottom);
      free(Azz_bottom);
   }
}

// =============== 3D ACOUSTIC PML CLASS =============== //
template<typename T>
PmlAcoustic3D<T>::PmlAcoustic3D(): Pml<T>() {
    /* Default constructor, not to be used */ 
    /* Minimal allocation to avoid free from crashing */
    P_left=(T *) malloc(1);
    P_right=(T *) malloc(1);
    P_top=(T *) malloc(1);
    P_bottom=(T *) malloc(1);
    P_front=(T *) malloc(1);
    P_back=(T *) malloc(1);
    Axx_left=(T *) malloc(1);
    Axx_right=(T *) malloc(1);
    Ayy_front=(T *) malloc(1);
    Ayy_back=(T *) malloc(1);
    Azz_top=(T *) malloc(1);
    Azz_bottom=(T *) malloc(1);
}

template<typename T>
PmlAcoustic3D<T>::PmlAcoustic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, ny_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    ny_pml= ny + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    P_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    P_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    P_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    P_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    P_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    P_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Axx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Axx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Ayy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Ayy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Azz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Azz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
}

template<typename T>
PmlAcoustic3D<T>::PmlAcoustic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt, const int dim, const bool low, const bool high): Pml<T>(Lpml, dt) {
    int nx_pml, ny_pml, nz_pml;
    nx_pml= nx;
    ny_pml= ny;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    switch (dim){
        case 0:
            if(low){
                P_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
                Axx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
                this->setApplypml(0,true);
            }
            if(high){
                P_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
                Axx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
                this->setApplypml(1,true);
            }

            P_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            P_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            P_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            P_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            Ayy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            Ayy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            Azz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            Azz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            this->setApplypml(2,true);
            this->setApplypml(3,true);
            this->setApplypml(4,true);
            this->setApplypml(5,true);
            break;
        case 1:
            if(low){
                P_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
                Ayy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
                this->setApplypml(2,true);
            }
            if(high){
                P_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
                Ayy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
                this->setApplypml(3,true);
            }

            P_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            Axx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            P_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            Axx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            P_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            P_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            Azz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            Azz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
            this->setApplypml(0,true);
            this->setApplypml(1,true);
            this->setApplypml(4,true);
            this->setApplypml(5,true);
            break;

        case 2:
            if(low){
                P_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
                Azz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
                this->setApplypml(4,true);
            }
            if(high){
                P_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
                Azz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
                this->setApplypml(5,true);
            }

            P_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            P_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            P_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            P_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            Axx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            Axx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
            Ayy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            Ayy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
            this->setApplypml(0,true);
            this->setApplypml(1,true);
            this->setApplypml(2,true);
            this->setApplypml(3,true);
            break;
        default:
            rs_error("PmlAcoustic2D::Invalid value for dim");
            break;
    }
}

template<typename T>
PmlAcoustic3D<T>::~PmlAcoustic3D() {
    /* Free variables */

    if(this->getApplypml(0)){
        free(P_left);
        free(Axx_left);

    }
    if(this->getApplypml(1)){
        free(P_right);
        free(Axx_right);
    }
    if(this->getApplypml(2)){
        free(P_front);
        free(Ayy_front);
    }
    if(this->getApplypml(3)){
        free(P_back);
        free(Ayy_back);
    }
    if(this->getApplypml(4)){
        free(P_top);
        free(Azz_top);
    }
    if(this->getApplypml(5)){
        free(P_bottom);
        free(Azz_bottom);
    }
}

// =============== 2D ELASTIC PML CLASS =============== //
template<typename T>
PmlElastic2D<T>::PmlElastic2D(const int nx, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    Sxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Sxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Sxzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Sxzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Szz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Szz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Sxzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Sxzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Vxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Vxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Vzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Vzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Vzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Vzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Vxz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Vxz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
}

template<typename T>
PmlElastic2D<T>::~PmlElastic2D() {
    /* Free variables */
    free(Sxx_left);
    free(Sxx_right);
    free(Szz_top);
    free(Szz_bottom);
    free(Sxzz_top);
    free(Sxzz_bottom);
    free(Sxzx_left);
    free(Sxzx_right);
    free(Vxx_left);
    free(Vxx_right);
    free(Vzx_left);
    free(Vzx_right);
    free(Vzz_top);
    free(Vzz_bottom);
    free(Vxz_top);
    free(Vxz_bottom);
}


// =============== 3D ELASTIC PML CLASS =============== //
template<typename T>
PmlElastic3D<T>::PmlElastic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, ny_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    ny_pml= ny + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    Sxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Sxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Sxzx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Sxzx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Sxyx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Sxyx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vzx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vzx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vyx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vyx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));

    Szz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Szz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Sxzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Sxzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Syzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Syzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vxz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vxz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vyz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vyz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));

    Syy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Syy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Sxyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Sxyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Syzy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Syzy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vxy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vxy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vzy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vzy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
}

template<typename T>
PmlElastic3D<T>::~PmlElastic3D() {
    /* Free variables */
    free(Sxx_left);
    free(Sxx_right);
    free(Syy_front);
    free(Syy_back);
    free(Szz_top);
    free(Szz_bottom);
    free(Sxzz_top);
    free(Sxzz_bottom);
    free(Syzz_top);
    free(Syzz_bottom);
    free(Sxzx_left);
    free(Sxzx_right);
    free(Sxyx_left);
    free(Sxyx_right);
    free(Sxyy_front);
    free(Sxyy_back);
    free(Syzy_front);
    free(Syzy_back);
    free(Vxx_left);
    free(Vxx_right);
    free(Vzx_left);
    free(Vzx_right);
    free(Vyx_left);
    free(Vyx_right);
    free(Vzz_top);
    free(Vzz_bottom);
    free(Vxz_top);
    free(Vxz_bottom);
    free(Vyz_top);
    free(Vyz_bottom);
    free(Vxy_front);
    free(Vxy_back);
    free(Vzy_front);
    free(Vzy_back);
    free(Vyy_front);
    free(Vyy_back);
}


// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Pml<float>;
template class PmlAcoustic2D<float>;
template class PmlAcoustic3D<float>;
template class PmlElastic2D<float>;
template class PmlElastic3D<float>;

template class Pml<double>;
template class PmlAcoustic2D<double>;
template class PmlAcoustic3D<double>;
template class PmlElastic2D<double>;
template class PmlElastic3D<double>;


}
