
// Interface to the diff gpu kernels

typedef struct nctempfloat1 { int d[1]; float *a;} nctempfloat1; 
typedef struct nctempfloat2 { int d[2]; float *a;} nctempfloat2; 
typedef struct nctempint1   { int d[1]; int *a;} nctempint1; 

int DiffDxminus(nctempfloat2 *df, nctempfloat2 *ddf, nctempfloat1 *dcoeffs, 
                float dx, int order);

int DiffDxplus(nctempfloat2 *df, nctempfloat2 *ddf, nctempfloat1 *dcoeffs, 
                float dx, int order);

int DiffDyminus(nctempfloat2 *df, nctempfloat2 *ddf, nctempfloat1 *dcoeffs, 
                float dx, int order);

int DiffDyplus(nctempfloat2 *df, nctempfloat2 *ddf, nctempfloat1 *dcoeffs, 
                float dx, int order);
