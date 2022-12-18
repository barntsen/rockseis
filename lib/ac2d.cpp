/*  Translated by epsc  version December 2021 */
typedef struct { float r; float i;} complex; 
typedef struct nctempfloat1 { int d[1]; float *a;} nctempfloat1; 
typedef struct nctempint1 { int d[1]; int *a;} nctempint1; 
typedef struct nctempchar1 { int d[1]; char *a;} nctempchar1; 
typedef struct nctempcomplex1 { int d[1]; complex *a;} nctempcomplex1; 
typedef struct nctempfloat2 { int d[2]; float *a;} nctempfloat2; 
typedef struct nctempint2 { int d[2]; int *a;} nctempint2; 
typedef struct nctempchar2 { int d[2]; char *a;} nctempchar2; 
typedef struct nctempcomplex2 { int d[2]; complex *a;} nctempcomplex2; 
typedef struct nctempfloat3 { int d[3]; float *a;} nctempfloat3; 
typedef struct nctempint3 { int d[3]; int *a;} nctempint3; 
typedef struct nctempchar3 { int d[3]; char *a;} nctempchar3; 
typedef struct nctempcomplex3 { int d[3]; complex *a;} nctempcomplex3; 
typedef struct nctempfloat4 { int d[4]; float *a;} nctempfloat4; 
typedef struct nctempint4 { int d[4]; int *a;} nctempint4; 
typedef struct nctempchar4 { int d[4]; char *a;} nctempchar4; 
typedef struct nctempcomplex4 { int d[4]; complex *a;} nctempcomplex4; 
#include <stdio.h>
extern "C" {
#include <stdlib.h>
#include <string.h>
}
#define NBLOCKS 1024
#define NTHREADS 1024
void *GpuNew(int n);
void *GpuDelete(void *f);
void *GpuError();
int Ac2dFwstepvx (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepvx2 (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepvz2 (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressx2 (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressz2 (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepvx (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepvx (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp5=Vx->d[0];nx =nctemp5;
int nctemp13=Vx->d[1];nz =nctemp13;
int nctemp21=Alftstag->d[0];lpml =nctemp21;
int nctemp27=0;
int nctemp29=nx;
int nctemp32=0;
int nctemp34=nz;
int nctemp25=(nctemp29-nctemp27)*(nctemp34-nctemp32);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp25;nctempno+=blockDim.x*gridDim.x){
iz=nctemp32+nctempno/(nctemp29-nctemp27);
ix=nctemp27+nctempno%(nctemp29-nctemp27);
{
int nctemp39=ix;
nctemp39=iz*Vx->d[0]+nctemp39;
int nctemp46=ix;
nctemp46=iz*Vx->d[0]+nctemp46;
int nctemp57=ix;
nctemp57=iz*Rx->d[0]+nctemp57;
float nctemp60 = dt * Rx->a[nctemp57];
int nctemp62=ix;
nctemp62=iz*df->d[0]+nctemp62;
float nctemp65 = nctemp60 * df->a[nctemp62];
float nctemp66 = Vx->a[nctemp46] + nctemp65;
Vx->a[nctemp39] =nctemp66;
}
}
}
int Ac2dFwstepvx (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvx<<<NBLOCKS,NTHREADS>>>(Pleft,Pright,Vx,Rx,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepvx2 (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepvx2 (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp71=Vx->d[0];nx =nctemp71;
int nctemp79=Vx->d[1];nz =nctemp79;
int nctemp87=Alftstag->d[0];lpml =nctemp87;
int nctemp94=0;
int nctemp97=1;
int nctemp91 = (Getapplypml->a[nctemp94] || Getapplypml->a[nctemp97]);
if(nctemp91)
{
int nctemp101=0;
int nctemp103=lpml;
int nctemp106=0;
int nctemp108=nz;
int nctemp99=(nctemp103-nctemp101)*(nctemp108-nctemp106);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp99;nctempno+=blockDim.x*gridDim.x){
iz=nctemp106+nctempno/(nctemp103-nctemp101);
ix=nctemp101+nctempno%(nctemp103-nctemp101);
{
int nctemp111=0;
if(Getapplypml->a[nctemp111])
{
int nctemp116 = (ix >= ix0);
int nctemp129 = ix0 + nx;
int nctemp121 = (ix < nctemp129);
int nctemp113 = (nctemp116 && nctemp121);
if(nctemp113)
{
int nctemp133=ix;
nctemp133=iz*Pleft->d[0]+nctemp133;
int nctemp143=ix;
int nctemp146=ix;
nctemp146=iz*Pleft->d[0]+nctemp146;
float nctemp149 = Blftstag->a[nctemp143] * Pleft->a[nctemp146];
int nctemp154=ix;
int nctemp162 = ix - ix0;
int nctemp157=nctemp162;
nctemp157=iz*df->d[0]+nctemp157;
float nctemp164 = Alftstag->a[nctemp154] * df->a[nctemp157];
float nctemp165 = nctemp149 + nctemp164;
Pleft->a[nctemp133] =nctemp165;
int nctemp174 = ix - ix0;
int nctemp169=nctemp174;
nctemp169=iz*Vx->d[0]+nctemp169;
int nctemp185 = ix - ix0;
int nctemp180=nctemp185;
nctemp180=iz*Vx->d[0]+nctemp180;
int nctemp200 = ix - ix0;
int nctemp195=nctemp200;
nctemp195=iz*Rx->d[0]+nctemp195;
float nctemp202 = dt * Rx->a[nctemp195];
int nctemp207=ix;
nctemp207=iz*Pleft->d[0]+nctemp207;
int nctemp214=ix;
int nctemp222 = ix - ix0;
int nctemp217=nctemp222;
nctemp217=iz*df->d[0]+nctemp217;
float nctemp224 = Clftstag->a[nctemp214] * df->a[nctemp217];
float nctemp225 = Pleft->a[nctemp207] + nctemp224;
float nctemp226 = nctemp202 * nctemp225;
float nctemp227 = Vx->a[nctemp180] - nctemp226;
Vx->a[nctemp169] =nctemp227;
}
}
int nctemp229=1;
if(Getapplypml->a[nctemp229])
{
int nctemp242 = ix + nxo;
int nctemp244 = nctemp242 - lpml;
i =nctemp244;
int nctemp248 = (i >= ix0);
int nctemp261 = ix0 + nx;
int nctemp253 = (i < nctemp261);
int nctemp245 = (nctemp248 && nctemp253);
if(nctemp245)
{
int nctemp265=ix;
nctemp265=iz*Pright->d[0]+nctemp265;
int nctemp275=ix;
int nctemp278=ix;
nctemp278=iz*Pright->d[0]+nctemp278;
float nctemp281 = Brbbstag->a[nctemp275] * Pright->a[nctemp278];
int nctemp286=ix;
int nctemp294 = i - ix0;
int nctemp289=nctemp294;
nctemp289=iz*df->d[0]+nctemp289;
float nctemp296 = Arbbstag->a[nctemp286] * df->a[nctemp289];
float nctemp297 = nctemp281 + nctemp296;
Pright->a[nctemp265] =nctemp297;
int nctemp306 = i - ix0;
int nctemp301=nctemp306;
nctemp301=iz*Vx->d[0]+nctemp301;
int nctemp317 = i - ix0;
int nctemp312=nctemp317;
nctemp312=iz*Vx->d[0]+nctemp312;
int nctemp332 = i - ix0;
int nctemp327=nctemp332;
nctemp327=iz*Rx->d[0]+nctemp327;
float nctemp334 = dt * Rx->a[nctemp327];
int nctemp339=ix;
nctemp339=iz*Pright->d[0]+nctemp339;
int nctemp346=ix;
int nctemp354 = i - ix0;
int nctemp349=nctemp354;
nctemp349=iz*df->d[0]+nctemp349;
float nctemp356 = Crbbstag->a[nctemp346] * df->a[nctemp349];
float nctemp357 = Pright->a[nctemp339] + nctemp356;
float nctemp358 = nctemp334 * nctemp357;
float nctemp359 = Vx->a[nctemp312] - nctemp358;
Vx->a[nctemp301] =nctemp359;
}
}
}
}
}
}
int Ac2dFwstepvx2 (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvx2<<<NBLOCKS,NTHREADS>>>(Pleft,Pright,Vx,Rx,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp364=Vz->d[0];nx =nctemp364;
int nctemp372=Vz->d[1];nz =nctemp372;
int nctemp380=Alftstag->d[0];lpml =nctemp380;
int nctemp386=0;
int nctemp388=nx;
int nctemp391=0;
int nctemp393=nz;
int nctemp384=(nctemp388-nctemp386)*(nctemp393-nctemp391);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp384;nctempno+=blockDim.x*gridDim.x){
iz=nctemp391+nctempno/(nctemp388-nctemp386);
ix=nctemp386+nctempno%(nctemp388-nctemp386);
{
int nctemp398=ix;
nctemp398=iz*Vz->d[0]+nctemp398;
int nctemp405=ix;
nctemp405=iz*Vz->d[0]+nctemp405;
int nctemp416=ix;
nctemp416=iz*Rz->d[0]+nctemp416;
float nctemp419 = dt * Rz->a[nctemp416];
int nctemp421=ix;
nctemp421=iz*df->d[0]+nctemp421;
float nctemp424 = nctemp419 * df->a[nctemp421];
float nctemp425 = Vz->a[nctemp405] + nctemp424;
Vz->a[nctemp398] =nctemp425;
}
}
}
int Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvz<<<NBLOCKS,NTHREADS>>>(Ptop,Pbot,Vz,Rz,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepvz2 (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepvz2 (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp430=Vz->d[0];nx =nctemp430;
int nctemp438=Vz->d[1];nz =nctemp438;
int nctemp446=Alftstag->d[0];lpml =nctemp446;
int nctemp453=4;
int nctemp456=5;
int nctemp450 = (Getapplypml->a[nctemp453] || Getapplypml->a[nctemp456]);
if(nctemp450)
{
int nctemp460=0;
int nctemp462=nx;
int nctemp465=0;
int nctemp467=lpml;
int nctemp458=(nctemp462-nctemp460)*(nctemp467-nctemp465);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp458;nctempno+=blockDim.x*gridDim.x){
iz=nctemp465+nctempno/(nctemp462-nctemp460);
ix=nctemp460+nctempno%(nctemp462-nctemp460);
{
int nctemp470=4;
if(Getapplypml->a[nctemp470])
{
int nctemp475 = (iz >= iz0);
int nctemp488 = iz0 + nz;
int nctemp480 = (iz < nctemp488);
int nctemp472 = (nctemp475 && nctemp480);
if(nctemp472)
{
int nctemp492=ix;
nctemp492=iz*Ptop->d[0]+nctemp492;
int nctemp502=iz;
int nctemp505=ix;
nctemp505=iz*Ptop->d[0]+nctemp505;
float nctemp508 = Blftstag->a[nctemp502] * Ptop->a[nctemp505];
int nctemp513=iz;
int nctemp516=ix;
int nctemp522 = iz - iz0;
nctemp516=nctemp522*df->d[0]+nctemp516;
float nctemp523 = Alftstag->a[nctemp513] * df->a[nctemp516];
float nctemp524 = nctemp508 + nctemp523;
Ptop->a[nctemp492] =nctemp524;
int nctemp528=ix;
int nctemp534 = iz - iz0;
nctemp528=nctemp534*Vz->d[0]+nctemp528;
int nctemp539=ix;
int nctemp545 = iz - iz0;
nctemp539=nctemp545*Vz->d[0]+nctemp539;
int nctemp554=ix;
int nctemp560 = iz - iz0;
nctemp554=nctemp560*Rz->d[0]+nctemp554;
float nctemp561 = dt * Rz->a[nctemp554];
int nctemp566=ix;
nctemp566=iz*Ptop->d[0]+nctemp566;
int nctemp573=iz;
int nctemp576=ix;
int nctemp582 = iz - iz0;
nctemp576=nctemp582*df->d[0]+nctemp576;
float nctemp583 = Clftstag->a[nctemp573] * df->a[nctemp576];
float nctemp584 = Ptop->a[nctemp566] + nctemp583;
float nctemp585 = nctemp561 * nctemp584;
float nctemp586 = Vz->a[nctemp539] - nctemp585;
Vz->a[nctemp528] =nctemp586;
}
}
int nctemp588=5;
if(Getapplypml->a[nctemp588])
{
int nctemp601 = iz + nzo;
int nctemp603 = nctemp601 - lpml;
i =nctemp603;
int nctemp610 = (i >= iz0);
int nctemp607 = (nctemp610 && i);
int nctemp620 = iz0 + nz;
int nctemp604 = (nctemp607 < nctemp620);
if(nctemp604)
{
int nctemp624=ix;
nctemp624=iz*Pbot->d[0]+nctemp624;
int nctemp634=iz;
int nctemp637=ix;
nctemp637=iz*Pbot->d[0]+nctemp637;
float nctemp640 = Brbbstag->a[nctemp634] * Pbot->a[nctemp637];
int nctemp645=iz;
int nctemp648=ix;
int nctemp654 = i - iz0;
nctemp648=nctemp654*df->d[0]+nctemp648;
float nctemp655 = Arbbstag->a[nctemp645] * df->a[nctemp648];
float nctemp656 = nctemp640 + nctemp655;
Pbot->a[nctemp624] =nctemp656;
int nctemp660=ix;
int nctemp666 = i - iz0;
nctemp660=nctemp666*Vz->d[0]+nctemp660;
int nctemp671=ix;
int nctemp677 = i - iz0;
nctemp671=nctemp677*Vz->d[0]+nctemp671;
int nctemp686=ix;
int nctemp692 = i - iz0;
nctemp686=nctemp692*Rz->d[0]+nctemp686;
float nctemp693 = dt * Rz->a[nctemp686];
int nctemp698=ix;
nctemp698=iz*Pbot->d[0]+nctemp698;
int nctemp705=iz;
int nctemp708=ix;
int nctemp714 = i - iz0;
nctemp708=nctemp714*df->d[0]+nctemp708;
float nctemp715 = Crbbstag->a[nctemp705] * df->a[nctemp708];
float nctemp716 = Pbot->a[nctemp698] + nctemp715;
float nctemp717 = nctemp693 * nctemp716;
float nctemp718 = Vz->a[nctemp671] - nctemp717;
Vz->a[nctemp660] =nctemp718;
}
}
}
}
}
}
int Ac2dFwstepvz2 (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvz2<<<NBLOCKS,NTHREADS>>>(Ptop,Pbot,Vz,Rz,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp723=P->d[0];nx =nctemp723;
int nctemp731=P->d[1];nz =nctemp731;
int nctemp739=Alft->d[0];lpml =nctemp739;
int nctemp745=0;
int nctemp747=nx;
int nctemp750=0;
int nctemp752=nz;
int nctemp743=(nctemp747-nctemp745)*(nctemp752-nctemp750);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp743;nctempno+=blockDim.x*gridDim.x){
iz=nctemp750+nctempno/(nctemp747-nctemp745);
ix=nctemp745+nctempno%(nctemp747-nctemp745);
{
int nctemp757=ix;
nctemp757=iz*P->d[0]+nctemp757;
int nctemp764=ix;
nctemp764=iz*P->d[0]+nctemp764;
int nctemp775=ix;
nctemp775=iz*L->d[0]+nctemp775;
float nctemp778 = dt * L->a[nctemp775];
int nctemp780=ix;
nctemp780=iz*df->d[0]+nctemp780;
float nctemp783 = nctemp778 * df->a[nctemp780];
float nctemp784 = P->a[nctemp764] + nctemp783;
P->a[nctemp757] =nctemp784;
}
}
}
int Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressx<<<NBLOCKS,NTHREADS>>>(Vxxleft,Vxxright,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepstressx2 (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepstressx2 (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp789=P->d[0];nx =nctemp789;
int nctemp797=P->d[1];nz =nctemp797;
int nctemp805=Alft->d[0];lpml =nctemp805;
int nctemp812=0;
int nctemp815=1;
int nctemp809 = (Getapplypml->a[nctemp812] || Getapplypml->a[nctemp815]);
if(nctemp809)
{
int nctemp819=0;
int nctemp821=lpml;
int nctemp824=0;
int nctemp826=nz;
int nctemp817=(nctemp821-nctemp819)*(nctemp826-nctemp824);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp817;nctempno+=blockDim.x*gridDim.x){
iz=nctemp824+nctempno/(nctemp821-nctemp819);
ix=nctemp819+nctempno%(nctemp821-nctemp819);
{
int nctemp829=0;
if(Getapplypml->a[nctemp829])
{
int nctemp837 = (ix >= ix0);
int nctemp834 = (nctemp837 && ix);
int nctemp847 = ix0 + nx;
int nctemp831 = (nctemp834 < nctemp847);
if(nctemp831)
{
int nctemp851=ix;
nctemp851=iz*Vxxleft->d[0]+nctemp851;
int nctemp861=ix;
int nctemp864=ix;
nctemp864=iz*Vxxleft->d[0]+nctemp864;
float nctemp867 = Blft->a[nctemp861] * Vxxleft->a[nctemp864];
int nctemp872=ix;
int nctemp880 = ix - ix0;
int nctemp875=nctemp880;
nctemp875=iz*df->d[0]+nctemp875;
float nctemp882 = Alft->a[nctemp872] * df->a[nctemp875];
float nctemp883 = nctemp867 + nctemp882;
Vxxleft->a[nctemp851] =nctemp883;
int nctemp892 = ix - ix0;
int nctemp887=nctemp892;
nctemp887=iz*P->d[0]+nctemp887;
int nctemp903 = ix - ix0;
int nctemp898=nctemp903;
nctemp898=iz*P->d[0]+nctemp898;
int nctemp918 = ix - ix0;
int nctemp913=nctemp918;
nctemp913=iz*L->d[0]+nctemp913;
float nctemp920 = dt * L->a[nctemp913];
int nctemp925=ix;
nctemp925=iz*Vxxleft->d[0]+nctemp925;
int nctemp932=ix;
int nctemp940 = ix - ix0;
int nctemp935=nctemp940;
nctemp935=iz*df->d[0]+nctemp935;
float nctemp942 = Clft->a[nctemp932] * df->a[nctemp935];
float nctemp943 = Vxxleft->a[nctemp925] + nctemp942;
float nctemp944 = nctemp920 * nctemp943;
float nctemp945 = P->a[nctemp898] - nctemp944;
P->a[nctemp887] =nctemp945;
}
}
int nctemp947=1;
if(Getapplypml->a[nctemp947])
{
int nctemp960 = ix + nxo;
int nctemp962 = nctemp960 - lpml;
i =nctemp962;
int nctemp969 = (i >= ix0);
int nctemp966 = (nctemp969 && i);
int nctemp979 = ix0 + nx;
int nctemp963 = (nctemp966 < nctemp979);
if(nctemp963)
{
int nctemp983=ix;
nctemp983=iz*Vxxright->d[0]+nctemp983;
int nctemp993=ix;
int nctemp996=ix;
nctemp996=iz*Vxxright->d[0]+nctemp996;
float nctemp999 = Brbb->a[nctemp993] * Vxxright->a[nctemp996];
int nctemp1004=ix;
int nctemp1012 = i - ix0;
int nctemp1007=nctemp1012;
nctemp1007=iz*df->d[0]+nctemp1007;
float nctemp1014 = Arbb->a[nctemp1004] * df->a[nctemp1007];
float nctemp1015 = nctemp999 + nctemp1014;
Vxxright->a[nctemp983] =nctemp1015;
int nctemp1024 = i - ix0;
int nctemp1019=nctemp1024;
nctemp1019=iz*P->d[0]+nctemp1019;
int nctemp1035 = i - ix0;
int nctemp1030=nctemp1035;
nctemp1030=iz*P->d[0]+nctemp1030;
int nctemp1050 = i - ix0;
int nctemp1045=nctemp1050;
nctemp1045=iz*L->d[0]+nctemp1045;
float nctemp1052 = dt * L->a[nctemp1045];
int nctemp1057=ix;
nctemp1057=iz*Vxxright->d[0]+nctemp1057;
int nctemp1064=ix;
int nctemp1072 = i - ix0;
int nctemp1067=nctemp1072;
nctemp1067=iz*df->d[0]+nctemp1067;
float nctemp1074 = Crbb->a[nctemp1064] * df->a[nctemp1067];
float nctemp1075 = Vxxright->a[nctemp1057] + nctemp1074;
float nctemp1076 = nctemp1052 * nctemp1075;
float nctemp1077 = P->a[nctemp1030] - nctemp1076;
P->a[nctemp1019] =nctemp1077;
}
}
}
}
}
}
int Ac2dFwstepstressx2 (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressx2<<<NBLOCKS,NTHREADS>>>(Vxxleft,Vxxright,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp1082=P->d[0];nx =nctemp1082;
int nctemp1090=P->d[1];nz =nctemp1090;
int nctemp1098=Alft->d[0];lpml =nctemp1098;
int nctemp1104=0;
int nctemp1106=nx;
int nctemp1109=0;
int nctemp1111=nz;
int nctemp1102=(nctemp1106-nctemp1104)*(nctemp1111-nctemp1109);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1102;nctempno+=blockDim.x*gridDim.x){
iz=nctemp1109+nctempno/(nctemp1106-nctemp1104);
ix=nctemp1104+nctempno%(nctemp1106-nctemp1104);
{
int nctemp1116=ix;
nctemp1116=iz*P->d[0]+nctemp1116;
int nctemp1123=ix;
nctemp1123=iz*P->d[0]+nctemp1123;
int nctemp1134=ix;
nctemp1134=iz*L->d[0]+nctemp1134;
float nctemp1137 = dt * L->a[nctemp1134];
int nctemp1139=ix;
nctemp1139=iz*df->d[0]+nctemp1139;
float nctemp1142 = nctemp1137 * df->a[nctemp1139];
float nctemp1143 = P->a[nctemp1123] + nctemp1142;
P->a[nctemp1116] =nctemp1143;
}
}
}
int Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressz<<<NBLOCKS,NTHREADS>>>(Vzztop,Vzzbot,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
__global__ void kernel_Ac2dFwstepstressz2 (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
__global__ void kernel_Ac2dFwstepstressz2 (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
int nx;
int nz;
int ix;
int iz;
int i;
int lpml;
int nctemp1148=P->d[0];nx =nctemp1148;
int nctemp1156=P->d[1];nz =nctemp1156;
int nctemp1164=Alft->d[0];lpml =nctemp1164;
int nctemp1171=4;
int nctemp1174=5;
int nctemp1168 = (Getapplypml->a[nctemp1171] || Getapplypml->a[nctemp1174]);
if(nctemp1168)
{
int nctemp1178=0;
int nctemp1180=nx;
int nctemp1183=0;
int nctemp1185=lpml;
int nctemp1176=(nctemp1180-nctemp1178)*(nctemp1185-nctemp1183);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1176;nctempno+=blockDim.x*gridDim.x){
iz=nctemp1183+nctempno/(nctemp1180-nctemp1178);
ix=nctemp1178+nctempno%(nctemp1180-nctemp1178);
{
int nctemp1188=4;
if(Getapplypml->a[nctemp1188])
{
int nctemp1196 = (iz >= iz0);
int nctemp1193 = (nctemp1196 && iz);
int nctemp1206 = iz0 + nz;
int nctemp1190 = (nctemp1193 < nctemp1206);
if(nctemp1190)
{
int nctemp1210=ix;
nctemp1210=iz*Vzztop->d[0]+nctemp1210;
int nctemp1220=iz;
int nctemp1223=ix;
nctemp1223=iz*Vzztop->d[0]+nctemp1223;
float nctemp1226 = Blft->a[nctemp1220] * Vzztop->a[nctemp1223];
int nctemp1231=iz;
int nctemp1234=ix;
int nctemp1240 = iz - iz0;
nctemp1234=nctemp1240*df->d[0]+nctemp1234;
float nctemp1241 = Alft->a[nctemp1231] * df->a[nctemp1234];
float nctemp1242 = nctemp1226 + nctemp1241;
Vzztop->a[nctemp1210] =nctemp1242;
int nctemp1246=ix;
int nctemp1252 = iz - iz0;
nctemp1246=nctemp1252*P->d[0]+nctemp1246;
int nctemp1257=ix;
int nctemp1263 = iz - iz0;
nctemp1257=nctemp1263*P->d[0]+nctemp1257;
int nctemp1272=ix;
int nctemp1278 = iz - iz0;
nctemp1272=nctemp1278*L->d[0]+nctemp1272;
float nctemp1279 = dt * L->a[nctemp1272];
int nctemp1284=ix;
nctemp1284=iz*Vzztop->d[0]+nctemp1284;
int nctemp1291=iz;
int nctemp1294=ix;
int nctemp1300 = iz - iz0;
nctemp1294=nctemp1300*df->d[0]+nctemp1294;
float nctemp1301 = Clft->a[nctemp1291] * df->a[nctemp1294];
float nctemp1302 = Vzztop->a[nctemp1284] + nctemp1301;
float nctemp1303 = nctemp1279 * nctemp1302;
float nctemp1304 = P->a[nctemp1257] - nctemp1303;
P->a[nctemp1246] =nctemp1304;
}
}
int nctemp1306=5;
if(Getapplypml->a[nctemp1306])
{
int nctemp1319 = iz + nzo;
int nctemp1321 = nctemp1319 - lpml;
i =nctemp1321;
int nctemp1328 = (i >= iz0);
int nctemp1325 = (nctemp1328 && i);
int nctemp1338 = iz0 + nz;
int nctemp1322 = (nctemp1325 < nctemp1338);
if(nctemp1322)
{
int nctemp1342=ix;
nctemp1342=iz*Vzzbot->d[0]+nctemp1342;
int nctemp1352=iz;
int nctemp1355=ix;
nctemp1355=iz*Vzzbot->d[0]+nctemp1355;
float nctemp1358 = Brbb->a[nctemp1352] * Vzzbot->a[nctemp1355];
int nctemp1363=iz;
int nctemp1366=ix;
int nctemp1372 = i - iz0;
nctemp1366=nctemp1372*df->d[0]+nctemp1366;
float nctemp1373 = Arbb->a[nctemp1363] * df->a[nctemp1366];
float nctemp1374 = nctemp1358 + nctemp1373;
Vzzbot->a[nctemp1342] =nctemp1374;
int nctemp1378=ix;
int nctemp1384 = i - iz0;
nctemp1378=nctemp1384*P->d[0]+nctemp1378;
int nctemp1389=ix;
int nctemp1395 = i - iz0;
nctemp1389=nctemp1395*P->d[0]+nctemp1389;
int nctemp1404=ix;
int nctemp1410 = i - iz0;
nctemp1404=nctemp1410*L->d[0]+nctemp1404;
float nctemp1411 = dt * L->a[nctemp1404];
int nctemp1416=ix;
nctemp1416=iz*Vzzbot->d[0]+nctemp1416;
int nctemp1423=iz;
int nctemp1426=ix;
int nctemp1432 = i - iz0;
nctemp1426=nctemp1432*df->d[0]+nctemp1426;
float nctemp1433 = Crbb->a[nctemp1423] * df->a[nctemp1426];
float nctemp1434 = Vzzbot->a[nctemp1416] + nctemp1433;
float nctemp1435 = nctemp1411 * nctemp1434;
float nctemp1436 = P->a[nctemp1389] - nctemp1435;
P->a[nctemp1378] =nctemp1436;
}
}
}
}
}
}
int Ac2dFwstepstressz2 (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressz2<<<NBLOCKS,NTHREADS>>>(Vzztop,Vzzbot,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
