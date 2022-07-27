//==============================================================================
// Differentiator Interface}
//==============================================================================

struct diff {
  int  l;        // Differentiator length
  int lmax;
  float [*,*] coeffs; 
  float [*]   w;// Differentiator weights 
}
struct diff DiffNew(int l){}
int DiffDxminus(float [*,*] A, float [*,*] dA, float [*] w,float dx,int l){}
int DiffDyminus(float [*,*] A, float [*,*] dA, float [*] w,float dx,int l){}
int DiffDxplus(float [*,*] A, float [*,*] dA,  float [*] w,float dx,int l){}
int DiffDyplus(float [*,*] A, float [*,*] dA, float [*]  w,float dx, int l){}
