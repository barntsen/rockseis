// Include statements
#include "pml.h"
#include "balloc.h"

namespace rockseis {

// =============== ABSTRACT PML CLASS =============== //
template<typename T>
Pml<T>::Pml() {
    Amax = AMAX;
    Smax = 1200.;
    Kmax = KMAX;
    dt = 0.0;
    Lpml=10;
    A_ltf = (T *) BallocNew(Lpml,sizeof(T));
    B_ltf = (T *) BallocNew(Lpml,sizeof(T));
    C_ltf = (T *) BallocNew(Lpml,sizeof(T));
    A_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    B_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    C_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    A_rbb = (T *) BallocNew(Lpml,sizeof(T));
    B_rbb = (T *) BallocNew(Lpml,sizeof(T));
    C_rbb = (T *) BallocNew(Lpml,sizeof(T));
    A_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));
    B_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));
    C_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));
 
    for(int i=0; i<6; i++) setApplypml(i,true);

    computeABC();

}

template<typename T>
Pml<T>::~Pml() {
    /* Free variables */
    BallocDelete(A_ltf);
    BallocDelete(B_ltf);
    BallocDelete(C_ltf);
    BallocDelete(A_ltf_stag);
    BallocDelete(B_ltf_stag);
    BallocDelete(C_ltf_stag);
    BallocDelete(A_rbb);
    BallocDelete(B_rbb);
    BallocDelete(C_rbb);
    BallocDelete(A_rbb_stag);
    BallocDelete(B_rbb_stag);
    BallocDelete(C_rbb_stag);
}

template<typename T>
Pml<T>::Pml(const int _Lpml, const T _dt) {
    Amax = AMAX;
    Smax = 1200.;
    Kmax = KMAX;
    dt = _dt;
    Lpml = _Lpml;
    A_ltf = (T *) BallocNew(Lpml,sizeof(T));
    B_ltf = (T *) BallocNew(Lpml,sizeof(T));
    C_ltf = (T *) BallocNew(Lpml,sizeof(T));
    A_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    B_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    C_ltf_stag = (T *) BallocNew(Lpml,sizeof(T));
    A_rbb = (T *) BallocNew(Lpml,sizeof(T));
    B_rbb = (T *) BallocNew(Lpml,sizeof(T));
    C_rbb = (T *) BallocNew(Lpml,sizeof(T));
    A_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));
    B_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));
    C_rbb_stag = (T *) BallocNew(Lpml,sizeof(T));

    for(int i=0; i<6; i++) setApplypml(i,true);
    
    computeABC();

}


template<typename T>
void Pml<T>::computeABC(int sign)
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
        func = (T) (Lpml - i - sign*0.5)/Lpml;
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
        func = (T) (i + 1.0 + sign*0.5)/Lpml;
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
    Vzz_top=(T *) malloc(1);
    Vzz_bottom=(T *) malloc(1);

}

template<typename T>
PmlAcoustic1D<T>::PmlAcoustic1D(const int Lpml, const T dt): Pml<T>(Lpml, dt) {

    /* Allocate variables */
    P_top=(T *) calloc(Lpml,sizeof(T));
    P_bottom=(T *) calloc(Lpml,sizeof(T));
    Vzz_top=(T *) calloc(Lpml,sizeof(T));
    Vzz_bottom=(T *) calloc(Lpml,sizeof(T));
}

template<typename T>
PmlAcoustic1D<T>::~PmlAcoustic1D() {
    /* Free variables */
    free(P_top);
    free(P_bottom);
    free(Vzz_top);
    free(Vzz_bottom);
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
    Vxx_left=(T *) malloc(1);
    Vxx_right=(T *) malloc(1);
    Vzz_top=(T *) malloc(1);
    Vzz_bottom=(T *) malloc(1);

}

template<typename T>
PmlAcoustic2D<T>::PmlAcoustic2D(const int nx, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    P_left=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
    P_right=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
    P_top=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
    P_bottom=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
    Vxx_left=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
    Vxx_right=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
    Vzz_top=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
    Vzz_bottom=(T *) BallocNew(nx_pml*Lpml,sizeof(T));

   //DEBUG
   printf("nx_pml,nz_pml,Lpml: %d %d %d \n", nx_pml,nz_pml,Lpml);
}

template<typename T>
PmlAcoustic2D<T>::PmlAcoustic2D(const int nx, const int nz, const int Lpml, const T dt, const bool *low, const bool *high): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    if(low[0]){
       P_left=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
       Vxx_left=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
       this->setApplypml(0,true);
    }
    if(high[0]){
       P_right=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
       Vxx_right=(T *) BallocNew(nz_pml*Lpml,sizeof(T));
       this->setApplypml(1,true);
    }
    if(low[1]){
       P_top=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
       Vzz_top=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
       this->setApplypml(4,true);
    }
    if(high[1]){
       P_bottom=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
       Vzz_bottom=(T *) BallocNew(nx_pml*Lpml,sizeof(T));
       this->setApplypml(5,true);
    }
}

template<typename T>
PmlAcoustic2D<T>::~PmlAcoustic2D() {
   /* Free variables */
   if(this->getApplypml(0)){
      BallocDelete(P_left);
      BallocDelete(Vxx_left);
   }
   if(this->getApplypml(1)){
      BallocDelete(P_right);
      BallocDelete(Vxx_right);
   }
   if(this->getApplypml(4)){
      BallocDelete(P_top);
      BallocDelete(Vzz_top);
   }
   if(this->getApplypml(5)){
      BallocDelete(P_bottom);
      BallocDelete(Vzz_bottom);
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
    Vxx_left=(T *) malloc(1);
    Vxx_right=(T *) malloc(1);
    Vyy_front=(T *) malloc(1);
    Vyy_back=(T *) malloc(1);
    Vzz_top=(T *) malloc(1);
    Vzz_bottom=(T *) malloc(1);
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
    Vxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
    Vyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
    Vzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
    Vzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
}

template<typename T>
PmlAcoustic3D<T>::PmlAcoustic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt, const bool *low, const bool *high): Pml<T>(Lpml, dt) {
    int nx_pml, ny_pml, nz_pml;
    nx_pml= nx;
    ny_pml= ny;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    if(low[0]){
       P_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(0,true);
    }
    if(high[0]){
       P_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(1,true);
    }

    if(low[1]){
       P_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(2,true);
    }
    if(high[1]){
       P_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(3,true);
    }

    if(low[2]){
       P_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       this->setApplypml(4,true);
    }
    if(high[2]){
       P_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       this->setApplypml(5,true);
    }
}

template<typename T>
PmlAcoustic3D<T>::~PmlAcoustic3D() {
    /* Free variables */

    if(this->getApplypml(0)){
        free(P_left);
        free(Vxx_left);

    }
    if(this->getApplypml(1)){
        free(P_right);
        free(Vxx_right);
    }
    if(this->getApplypml(2)){
        free(P_front);
        free(Vyy_front);
    }
    if(this->getApplypml(3)){
        free(P_back);
        free(Vyy_back);
    }
    if(this->getApplypml(4)){
        free(P_top);
        free(Vzz_top);
    }
    if(this->getApplypml(5)){
        free(P_bottom);
        free(Vzz_bottom);
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
PmlElastic2D<T>::PmlElastic2D(const int nx, const int nz, const int Lpml, const T dt, const bool *low, const bool *high): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    if(low[0]){
       Sxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       this->setApplypml(0,true);
    }
    if(high[0]){
       Sxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       this->setApplypml(1,true);
    }
    if(low[1]){
       Szz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Sxzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vxz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       this->setApplypml(4,true);
    }
    if(high[1]){
       Szz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Sxzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vxz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       this->setApplypml(5,true);
    }
}

template<typename T>
PmlElastic2D<T>::~PmlElastic2D() {
    /* Free variables */
   if(this->getApplypml(0)){
      free(Sxx_left);
      free(Vxx_left);
      free(Sxzx_left);
      free(Vzx_left);
   }
   if(this->getApplypml(1)){
      free(Sxx_right);
      free(Sxzx_right);
      free(Vxx_right);
      free(Vzx_right);
   }
   if(this->getApplypml(4)){
      free(Szz_top);
      free(Sxzz_top);
      free(Vzz_top);
      free(Vxz_top);
   }
   if(this->getApplypml(5)){
      free(Szz_bottom);
      free(Sxzz_bottom);
      free(Vzz_bottom);
      free(Vxz_bottom);
   }
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
PmlElastic3D<T>::PmlElastic3D(const int nx, const int ny, const int nz, const int Lpml, const T dt, const bool *low, const bool *high): Pml<T>(Lpml, dt) {
    int nx_pml, ny_pml, nz_pml;
    nx_pml= nx;
    ny_pml= ny;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    if(low[0]){
       Sxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Sxzx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Sxyx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vxx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vzx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vyx_left=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(0,true);
    }
    if(high[0]){
       Sxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Sxzx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Sxyx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vxx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vzx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       Vyx_right=(T *) calloc(ny_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(1,true);
    }

    if(low[1]){
       Syy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Sxyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Syzy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vxy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vzy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vyy_front=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(2,true);
    }
    if(high[1]){
       Syy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Sxyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Syzy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vxy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vzy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       Vyy_back=(T *) calloc(nx_pml*nz_pml*Lpml,sizeof(T));
       this->setApplypml(3,true);
    }

    if(low[2]){
       Szz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Sxzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Syzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vzz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vxz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vyz_top=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       this->setApplypml(4,true);
    }
    if(high[2]){
       Szz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Sxzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Syzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vzz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vxz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       Vyz_bottom=(T *) calloc(nx_pml*ny_pml*Lpml,sizeof(T));
       this->setApplypml(5,true);
    }

}


template<typename T>
PmlElastic3D<T>::~PmlElastic3D() {
    /* Free variables */
   if(this->getApplypml(0)){
      free(Sxx_left);
      free(Sxzx_left);
      free(Sxyx_left);
      free(Vxx_left);
      free(Vzx_left);
      free(Vyx_left);
   }
   if(this->getApplypml(1)){
      free(Sxx_right);
      free(Sxzx_right);
      free(Sxyx_right);
      free(Vxx_right);
      free(Vzx_right);
      free(Vyx_right);
   }
   if(this->getApplypml(2)){
      free(Syy_front);
      free(Sxyy_front);
      free(Syzy_front);
      free(Vxy_front);
      free(Vzy_front);
      free(Vyy_front);
   }
   if(this->getApplypml(3)){
      free(Sxyy_back);
      free(Syzy_back);
      free(Vxy_back);
      free(Vzy_back);
      free(Vyy_back);
      free(Syy_back);
   }
   if(this->getApplypml(4)){
      free(Szz_top);
      free(Sxzz_top);
      free(Syzz_top);
      free(Vzz_top);
      free(Vxz_top);
      free(Vyz_top);
   }
   if(this->getApplypml(5)){
      free(Szz_bottom);
      free(Sxzz_bottom);
      free(Syzz_bottom);
      free(Vzz_bottom);
      free(Vxz_bottom);
      free(Vyz_bottom);
   }
}

// =============== 2D POROELASTIC PML CLASS =============== //
template<typename T>
PmlPoroelastic2D<T>::PmlPoroelastic2D(const int nx, const int nz, const int Lpml, const T dt): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx + 2*Lpml;
    nz_pml= nz + 2*Lpml;

    /* Allocate variables */
    P_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    P_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    P_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    P_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Qxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Qxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
    Qzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
    Qzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
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
PmlPoroelastic2D<T>::PmlPoroelastic2D(const int nx, const int nz, const int Lpml, const T dt, const bool *low, const bool *high): Pml<T>(Lpml, dt) {
    int nx_pml, nz_pml;
    nx_pml= nx;
    nz_pml= nz;

    int i;
    for (i=0; i<6; i++) this->setApplypml(i, false);

    /* Allocate variables */
    if(low[0]){
       P_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Qxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vxx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vzx_left=(T *) calloc(nz_pml*Lpml,sizeof(T));
       this->setApplypml(0,true);
    }
    if(high[0]){
       P_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Sxzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Qxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vxx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       Vzx_right=(T *) calloc(nz_pml*Lpml,sizeof(T));
       this->setApplypml(1,true);
    }
    if(low[1]){
       P_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Szz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Sxzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Qzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vzz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vxz_top=(T *) calloc(nx_pml*Lpml,sizeof(T));
       this->setApplypml(4,true);
    }
    if(high[1]){
       P_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Szz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Sxzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Qzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vzz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       Vxz_bottom=(T *) calloc(nx_pml*Lpml,sizeof(T));
       this->setApplypml(5,true);
    }
}

template<typename T>
PmlPoroelastic2D<T>::~PmlPoroelastic2D() {
    /* Free variables */
   if(this->getApplypml(0)){
      free(P_left);
      free(Sxx_left);
      free(Vxx_left);
      free(Qxx_left);
      free(Sxzx_left);
      free(Vzx_left);
   }
   if(this->getApplypml(1)){
      free(P_right);
      free(Sxx_right);
      free(Sxzx_right);
      free(Qxx_right);
      free(Vxx_right);
      free(Vzx_right);
   }
   if(this->getApplypml(4)){
      free(P_top);
      free(Szz_top);
      free(Sxzz_top);
      free(Qzz_top);
      free(Vzz_top);
      free(Vxz_top);
   }
   if(this->getApplypml(5)){
      free(P_bottom);
      free(Szz_bottom);
      free(Sxzz_bottom);
      free(Qzz_bottom);
      free(Vzz_bottom);
      free(Vxz_bottom);
   }
}

// =============== INITIALIZING TEMPLATE CLASSES =============== //
template class Pml<float>;
template class PmlAcoustic2D<float>;
template class PmlAcoustic3D<float>;
template class PmlElastic2D<float>;
template class PmlElastic3D<float>;
template class PmlPoroelastic2D<float>;

template class Pml<double>;
template class PmlAcoustic2D<double>;
template class PmlAcoustic3D<double>;
template class PmlElastic2D<double>;
template class PmlElastic3D<double>;
template class PmlPoroelastic2D<double>;


}
