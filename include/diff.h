
typedef struct nctempfloat1 { int d[1]; float *a;} nctempfloat1; 
typedef struct nctempfloat2 { int d[2]; float *a;} nctempfloat2; 

int DiffDxminus (float *f, float *df, float *w, float dx, int l, int nx, int ny);
int DiffDyminus (float *f, float *df, float *w, float dx, int l, int nx, int ny);
int DiffDxplus (float *f, float *df, float *w, float dx, int l, int nx, int ny);
int DiffDyplus (float *f, float *df, float *w, float dx, int l, int nx, int ny);
