// Interface to the diffc functions


int DiffDxplusc(float *f, float *df, float *coeffs,
               float dx, int order, int nx, int nz);
 
int DiffDxminusc(float *f, float *df, float *coeffs,
               float dx, int order, int nx, int nz);

int DiffDyplusc(float *f, float *df, float *coeffs,
               float dx, int order, int nx, int nz);

int DiffDyminusc(float *f, float *df, float *coeffs,
               float dx, int order, int nx, int nz);
