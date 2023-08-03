              
 
int Ac2dFwstepvxc(float *Pleft, float *Pright,   
              float *Vx,float *Rx, float *df,    
              float *Alftstag,float *Blftstag,float *Clftstag, 
              float *Arbbstag,float *Brbbstag,float *Crbbstag, 
              int * Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng); 
              
int Ac2dFwstepvzc(float *Ptop, float *Pbot,   
              float *Vz,float *Rx, float *df,    
              float *Alftstag,float *Blftstag,float *Clftstag, 
              float *Arbbstag,float *Brbbstag,float *Crbbstag, 
              int * Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng); 

int Ac2dFwstepstressxc(float *Vxxleft, float *Vxxright,   
              float *P,float *L, float *df,    
              float *Alft,float *Blft,float *Clft, 
              float *Arbb,float *Brbb,float *Crbb, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng); 

int Ac2dFwstepstresszc(float *Vzztop, float *Vzzbot,   
              float *P,float *L, float *df,    
              float *Alft,float *Blft,float *Clft, 
              float *Arbb,float *Brbb,float *Crbb, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng); 


int Ac2dXcorrc(float *Vp,  int nx, int nz, int padr, int pads, 
              float *Rho, float *Rx, float *Rz,int nxs, int nys, 
	      int nxr, int nyr, float *wsp, float *wrx, float *wrz, 
              float *vpgraddata, float *rhograddata, 
	      float dx, float dz, int srcilumset, float *srcilumdata);  

int Ac2dMemcpyc( void * s, void *t, int n);
