
// Interface to the ac2d gpu kernels

typedef struct nctempfloat1 { int d[1]; float *a;} nctempfloat1; 
typedef struct nctempfloat2 { int d[2]; float *a;} nctempfloat2; 
typedef struct nctempint1   { int d[1]; int *a;}   nctempint1; 
typedef struct nctempchar1  { int d[1]; char *a;}  nctempchar1; 


int Ac2dFwstepvx(nctempfloat2 *Pleft, nctempfloat2 *Pright,   
         nctempfloat2 *Vx,      nctempfloat2 *Rx, nctempfloat2 *df,    
         nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag, 
         nctempfloat1 *Arbbstag,nctempfloat1*Brbbstag, nctempfloat1 *Crbbstag, 
         nctempint1   *Getapplypml,
         int ix0,int iz0,int nxo,int nzo,float dt); 

int Ac2dFwstepvx2(nctempfloat2 *Pleft, nctempfloat2 *Pright,   
         nctempfloat2 *Vx,      nctempfloat2 *Rx, nctempfloat2 *df,    
         nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag, 
         nctempfloat1 *Arbbstag,nctempfloat1*Brbbstag, nctempfloat1 *Crbbstag, 
         nctempint1   *Getapplypml,
         int ix0,int iz0,int nxo,int nzo,float dt); 
              
int Ac2dFwstepvz(nctempfloat2 *Ptop, nctempfloat2 *Pbot,   
         nctempfloat2 *Vz,      nctempfloat2 *Rz, nctempfloat2 *df,    
         nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag, 
         nctempfloat1 *Arbbstag,nctempfloat1*Brbbstag, nctempfloat1 *Crbbstag, 
         nctempint1   *Getapplypml,
         int ix0,int iz0,int nxo,int nzo,float dt); 

int Ac2dFwstepvz2(nctempfloat2 *Ptop, nctempfloat2 *Pbot,   
         nctempfloat2 *Vz,      nctempfloat2 *Rz, nctempfloat2 *df,    
         nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag, 
         nctempfloat1 *Arbbstag,nctempfloat1*Brbbstag, nctempfloat1 *Crbbstag, 
         nctempint1   *Getapplypml,
         int ix0,int iz0,int nxo,int nzo,float dt); 

int Ac2dFwstepstressx(nctempfloat2 *Vxxleft, nctempfloat2 *Vxxright,   
              nctempfloat2 *P,nctempfloat2 *L, nctempfloat2 *df,    
              nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft, 
              nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb, 
              nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt); 

int Ac2dFwstepstressx2(nctempfloat2 *Vxxleft, nctempfloat2 *Vxxright,   
              nctempfloat2 *P,nctempfloat2 *L, nctempfloat2 *df,    
              nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft, 
              nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb, 
              nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt); 
              
int Ac2dFwstepstressz(nctempfloat2 *Vzztop, nctempfloat2 *Vzzbot,   
              nctempfloat2 *P,nctempfloat2 *L, nctempfloat2 *df,    
              nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft, 
              nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1*Crbb, 
              nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt); 
int Ac2dFwstepstressz2(nctempfloat2 *Vzztop, nctempfloat2 *Vzzbot,   
              nctempfloat2 *P,nctempfloat2 *L, nctempfloat2 *df,    
              nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft, 
              nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1*Crbb, 
              nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt); 

int Ac2dXcorr(nctempfloat2 *Vp, int padr, int pads, 
              nctempfloat2 *Rho, nctempfloat2 *dRx,  nctempfloat2 *dRz, 
	      nctempfloat2 *dwsp, nctempfloat2 *dwrx, nctempfloat2 *dwrz, 
              nctempfloat2 *dvpgraddata, nctempfloat2 *drhograddata, 
	      float dx, float dz, int srcilumset, nctempfloat2 *dsrcilumdata);  

int Ac2dMemcpy(nctempchar1 *s, nctempchar1 *t);
