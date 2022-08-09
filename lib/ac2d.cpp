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
int Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
int Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt);
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
int nctemp70=0;
int nctemp73=1;
int nctemp67 = (Getapplypml->a[nctemp70] || Getapplypml->a[nctemp73]);
if(nctemp67)
{
int nctemp77=0;
int nctemp79=lpml;
int nctemp82=0;
int nctemp84=nz;
int nctemp75=(nctemp79-nctemp77)*(nctemp84-nctemp82);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp75;nctempno+=blockDim.x*gridDim.x){
iz=nctemp82+nctempno/(nctemp79-nctemp77);
ix=nctemp77+nctempno%(nctemp79-nctemp77);
{
int nctemp87=0;
if(Getapplypml->a[nctemp87])
{
int nctemp95 = (ix >= ix0);
int nctemp92 = (nctemp95 && ix);
int nctemp105 = ix0 + nx;
int nctemp89 = (nctemp92 < nctemp105);
if(nctemp89)
{
int nctemp109=ix;
nctemp109=iz*Pleft->d[0]+nctemp109;
int nctemp119=ix;
int nctemp122=ix;
nctemp122=iz*Pleft->d[0]+nctemp122;
float nctemp125 = Blftstag->a[nctemp119] * Pleft->a[nctemp122];
int nctemp130=ix;
int nctemp138 = ix - ix0;
int nctemp133=nctemp138;
nctemp133=iz*df->d[0]+nctemp133;
float nctemp140 = Alftstag->a[nctemp130] * df->a[nctemp133];
float nctemp141 = nctemp125 + nctemp140;
Pleft->a[nctemp109] =nctemp141;
int nctemp150 = ix - ix0;
int nctemp145=nctemp150;
nctemp145=iz*Vx->d[0]+nctemp145;
int nctemp161 = ix - iz0;
int nctemp156=nctemp161;
nctemp156=iz*Vx->d[0]+nctemp156;
int nctemp176 = ix - ix0;
int nctemp171=nctemp176;
nctemp171=iz*Rx->d[0]+nctemp171;
float nctemp178 = dt * Rx->a[nctemp171];
int nctemp183=ix;
nctemp183=iz*Pleft->d[0]+nctemp183;
int nctemp190=ix;
int nctemp198 = ix - ix0;
int nctemp193=nctemp198;
nctemp193=iz*df->d[0]+nctemp193;
float nctemp200 = Clftstag->a[nctemp190] * df->a[nctemp193];
float nctemp201 = Pleft->a[nctemp183] + nctemp200;
float nctemp202 = nctemp178 * nctemp201;
float nctemp203 = Vx->a[nctemp156] - nctemp202;
Vx->a[nctemp145] =nctemp203;
}
}
int nctemp205=1;
if(Getapplypml->a[nctemp205])
{
int nctemp218 = ix + nxo;
int nctemp220 = nctemp218 - lpml;
i =nctemp220;
int nctemp227 = (i >= ix0);
int nctemp224 = (nctemp227 && i);
int nctemp237 = ix0 + nx;
int nctemp221 = (nctemp224 < nctemp237);
if(nctemp221)
{
int nctemp241=ix;
nctemp241=iz*Pright->d[0]+nctemp241;
int nctemp251=ix;
int nctemp254=ix;
nctemp254=iz*Pright->d[0]+nctemp254;
float nctemp257 = Brbbstag->a[nctemp251] * Pright->a[nctemp254];
int nctemp262=ix;
int nctemp270 = i - ix0;
int nctemp265=nctemp270;
nctemp265=iz*df->d[0]+nctemp265;
float nctemp272 = Arbbstag->a[nctemp262] * df->a[nctemp265];
float nctemp273 = nctemp257 + nctemp272;
Pright->a[nctemp241] =nctemp273;
int nctemp282 = i - ix0;
int nctemp277=nctemp282;
nctemp277=iz*Vx->d[0]+nctemp277;
int nctemp293 = ix - ix0;
int nctemp288=nctemp293;
nctemp288=iz*Vx->d[0]+nctemp288;
int nctemp308 = i - ix0;
int nctemp303=nctemp308;
nctemp303=iz*Rx->d[0]+nctemp303;
float nctemp310 = dt * Rx->a[nctemp303];
int nctemp315=ix;
nctemp315=iz*Pright->d[0]+nctemp315;
int nctemp322=ix;
int nctemp330 = i - ix0;
int nctemp325=nctemp330;
nctemp325=iz*df->d[0]+nctemp325;
float nctemp332 = Crbbstag->a[nctemp322] * df->a[nctemp325];
float nctemp333 = Pright->a[nctemp315] + nctemp332;
float nctemp334 = nctemp310 * nctemp333;
float nctemp335 = Vx->a[nctemp288] - nctemp334;
Vx->a[nctemp277] =nctemp335;
}
}
}
}
}
}
int Ac2dFwstepvx (nctempfloat2 *Pleft,nctempfloat2 *Pright,nctempfloat2 *Vx,nctempfloat2 *Rx,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvx<<<NBLOCKS,NTHREADS>>>(Pleft,Pright,Vx,Rx,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
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
int nctemp340=Vz->d[0];nx =nctemp340;
int nctemp348=Vz->d[1];nz =nctemp348;
int nctemp356=Alftstag->d[0];lpml =nctemp356;
int nctemp362=0;
int nctemp364=nx;
int nctemp367=0;
int nctemp369=nz;
int nctemp360=(nctemp364-nctemp362)*(nctemp369-nctemp367);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp360;nctempno+=blockDim.x*gridDim.x){
iz=nctemp367+nctempno/(nctemp364-nctemp362);
ix=nctemp362+nctempno%(nctemp364-nctemp362);
{
int nctemp374=ix;
nctemp374=iz*Vz->d[0]+nctemp374;
int nctemp381=ix;
nctemp381=iz*Vz->d[0]+nctemp381;
int nctemp392=ix;
nctemp392=iz*Rz->d[0]+nctemp392;
float nctemp395 = dt * Rz->a[nctemp392];
int nctemp397=ix;
nctemp397=iz*df->d[0]+nctemp397;
float nctemp400 = nctemp395 * df->a[nctemp397];
float nctemp401 = Vz->a[nctemp381] + nctemp400;
Vz->a[nctemp374] =nctemp401;
}
}
int nctemp405=4;
int nctemp408=5;
int nctemp402 = (Getapplypml->a[nctemp405] || Getapplypml->a[nctemp408]);
if(nctemp402)
{
int nctemp412=0;
int nctemp414=lpml;
int nctemp417=0;
int nctemp419=nx;
int nctemp410=(nctemp414-nctemp412)*(nctemp419-nctemp417);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp410;nctempno+=blockDim.x*gridDim.x){
ix=nctemp417+nctempno/(nctemp414-nctemp412);
iz=nctemp412+nctempno%(nctemp414-nctemp412);
{
int nctemp422=4;
if(Getapplypml->a[nctemp422])
{
int nctemp427 = (iz >= iz0);
int nctemp440 = iz0 + nz;
int nctemp432 = (iz < nctemp440);
int nctemp424 = (nctemp427 && nctemp432);
if(nctemp424)
{
int nctemp444=ix;
nctemp444=iz*Ptop->d[0]+nctemp444;
int nctemp454=iz;
int nctemp457=ix;
nctemp457=iz*Ptop->d[0]+nctemp457;
float nctemp460 = Blftstag->a[nctemp454] * Ptop->a[nctemp457];
int nctemp465=iz;
int nctemp468=ix;
int nctemp474 = iz - iz0;
nctemp468=nctemp474*df->d[0]+nctemp468;
float nctemp475 = Alftstag->a[nctemp465] * df->a[nctemp468];
float nctemp476 = nctemp460 + nctemp475;
Ptop->a[nctemp444] =nctemp476;
int nctemp480=ix;
int nctemp486 = iz - iz0;
nctemp480=nctemp486*Vz->d[0]+nctemp480;
int nctemp491=ix;
int nctemp497 = iz - iz0;
nctemp491=nctemp497*Vz->d[0]+nctemp491;
int nctemp506=ix;
int nctemp512 = iz - iz0;
nctemp506=nctemp512*Rz->d[0]+nctemp506;
float nctemp513 = dt * Rz->a[nctemp506];
int nctemp518=ix;
nctemp518=iz*Ptop->d[0]+nctemp518;
int nctemp525=iz;
int nctemp528=ix;
int nctemp534 = iz - iz0;
nctemp528=nctemp534*df->d[0]+nctemp528;
float nctemp535 = Clftstag->a[nctemp525] * df->a[nctemp528];
float nctemp536 = Ptop->a[nctemp518] + nctemp535;
float nctemp537 = nctemp513 * nctemp536;
float nctemp538 = Vz->a[nctemp491] - nctemp537;
Vz->a[nctemp480] =nctemp538;
}
}
int nctemp540=5;
if(Getapplypml->a[nctemp540])
{
int nctemp553 = iz + nzo;
int nctemp555 = nctemp553 - lpml;
i =nctemp555;
int nctemp562 = (i >= iz0);
int nctemp559 = (nctemp562 && i);
int nctemp572 = iz0 + nz;
int nctemp556 = (nctemp559 < nctemp572);
if(nctemp556)
{
int nctemp576=ix;
nctemp576=iz*Pbot->d[0]+nctemp576;
int nctemp586=iz;
int nctemp589=ix;
nctemp589=iz*Pbot->d[0]+nctemp589;
float nctemp592 = Brbbstag->a[nctemp586] * Pbot->a[nctemp589];
int nctemp597=iz;
int nctemp600=ix;
int nctemp606 = i - iz0;
nctemp600=nctemp606*df->d[0]+nctemp600;
float nctemp607 = Arbbstag->a[nctemp597] * df->a[nctemp600];
float nctemp608 = nctemp592 + nctemp607;
Pbot->a[nctemp576] =nctemp608;
int nctemp612=ix;
int nctemp618 = i - iz0;
nctemp612=nctemp618*Vz->d[0]+nctemp612;
int nctemp623=ix;
int nctemp629 = i - iz0;
nctemp623=nctemp629*Vz->d[0]+nctemp623;
int nctemp638=ix;
int nctemp644 = i - iz0;
nctemp638=nctemp644*Rz->d[0]+nctemp638;
float nctemp645 = dt * Rz->a[nctemp638];
int nctemp650=ix;
nctemp650=iz*Pbot->d[0]+nctemp650;
int nctemp657=iz;
int nctemp660=ix;
int nctemp666 = i - iz0;
nctemp660=nctemp666*df->d[0]+nctemp660;
float nctemp667 = Crbbstag->a[nctemp657] * df->a[nctemp660];
float nctemp668 = Pbot->a[nctemp650] + nctemp667;
float nctemp669 = nctemp645 * nctemp668;
float nctemp670 = Vz->a[nctemp623] - nctemp669;
Vz->a[nctemp612] =nctemp670;
}
}
}
}
}
}
int Ac2dFwstepvz (nctempfloat2 *Ptop,nctempfloat2 *Pbot,nctempfloat2 *Vz,nctempfloat2 *Rz,nctempfloat2 *df,nctempfloat1 *Alftstag,nctempfloat1 *Blftstag,nctempfloat1 *Clftstag,nctempfloat1 *Arbbstag,nctempfloat1 *Brbbstag,nctempfloat1 *Crbbstag,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepvz<<<NBLOCKS,NTHREADS>>>(Ptop,Pbot,Vz,Rz,df,Alftstag,Blftstag,Clftstag,Arbbstag,Brbbstag,Crbbstag,Getapplypml,ix0,iz0,nxo,nzo,dt);
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
int nctemp675=P->d[0];nx =nctemp675;
int nctemp683=P->d[1];nz =nctemp683;
int nctemp691=Alft->d[0];lpml =nctemp691;
int nctemp697=0;
int nctemp699=nx;
int nctemp702=0;
int nctemp704=nz;
int nctemp695=(nctemp699-nctemp697)*(nctemp704-nctemp702);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp695;nctempno+=blockDim.x*gridDim.x){
iz=nctemp702+nctempno/(nctemp699-nctemp697);
ix=nctemp697+nctempno%(nctemp699-nctemp697);
{
int nctemp709=ix;
nctemp709=iz*P->d[0]+nctemp709;
int nctemp716=ix;
nctemp716=iz*P->d[0]+nctemp716;
int nctemp727=ix;
nctemp727=iz*L->d[0]+nctemp727;
float nctemp730 = dt * L->a[nctemp727];
int nctemp732=ix;
nctemp732=iz*df->d[0]+nctemp732;
float nctemp735 = nctemp730 * df->a[nctemp732];
float nctemp736 = P->a[nctemp716] + nctemp735;
P->a[nctemp709] =nctemp736;
}
}
int nctemp740=0;
int nctemp743=1;
int nctemp737 = (Getapplypml->a[nctemp740] || Getapplypml->a[nctemp743]);
if(nctemp737)
{
int nctemp747=0;
int nctemp749=lpml;
int nctemp752=0;
int nctemp754=nz;
int nctemp745=(nctemp749-nctemp747)*(nctemp754-nctemp752);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp745;nctempno+=blockDim.x*gridDim.x){
iz=nctemp752+nctempno/(nctemp749-nctemp747);
ix=nctemp747+nctempno%(nctemp749-nctemp747);
{
int nctemp757=0;
if(Getapplypml->a[nctemp757])
{
int nctemp765 = (ix >= ix0);
int nctemp762 = (nctemp765 && ix);
int nctemp775 = ix0 + nx;
int nctemp759 = (nctemp762 < nctemp775);
if(nctemp759)
{
int nctemp779=ix;
nctemp779=iz*Vxxleft->d[0]+nctemp779;
int nctemp789=ix;
int nctemp792=ix;
nctemp792=iz*Vxxleft->d[0]+nctemp792;
float nctemp795 = Blft->a[nctemp789] * Vxxleft->a[nctemp792];
int nctemp800=ix;
int nctemp808 = ix - ix0;
int nctemp803=nctemp808;
nctemp803=iz*df->d[0]+nctemp803;
float nctemp810 = Alft->a[nctemp800] * df->a[nctemp803];
float nctemp811 = nctemp795 + nctemp810;
Vxxleft->a[nctemp779] =nctemp811;
int nctemp820 = ix - ix0;
int nctemp815=nctemp820;
nctemp815=iz*P->d[0]+nctemp815;
int nctemp831 = ix - iz0;
int nctemp826=nctemp831;
nctemp826=iz*P->d[0]+nctemp826;
int nctemp846 = ix - ix0;
int nctemp841=nctemp846;
nctemp841=iz*L->d[0]+nctemp841;
float nctemp848 = dt * L->a[nctemp841];
int nctemp853=ix;
nctemp853=iz*Vxxleft->d[0]+nctemp853;
int nctemp860=ix;
int nctemp868 = ix - ix0;
int nctemp863=nctemp868;
nctemp863=iz*df->d[0]+nctemp863;
float nctemp870 = Clft->a[nctemp860] * df->a[nctemp863];
float nctemp871 = Vxxleft->a[nctemp853] + nctemp870;
float nctemp872 = nctemp848 * nctemp871;
float nctemp873 = P->a[nctemp826] - nctemp872;
P->a[nctemp815] =nctemp873;
}
}
int nctemp875=1;
if(Getapplypml->a[nctemp875])
{
int nctemp888 = ix + nxo;
int nctemp890 = nctemp888 - lpml;
i =nctemp890;
int nctemp897 = (i >= ix0);
int nctemp894 = (nctemp897 && i);
int nctemp907 = ix0 + nx;
int nctemp891 = (nctemp894 < nctemp907);
if(nctemp891)
{
int nctemp911=ix;
nctemp911=iz*Vxxright->d[0]+nctemp911;
int nctemp921=ix;
int nctemp924=ix;
nctemp924=iz*Vxxright->d[0]+nctemp924;
float nctemp927 = Brbb->a[nctemp921] * Vxxright->a[nctemp924];
int nctemp932=ix;
int nctemp940 = i - ix0;
int nctemp935=nctemp940;
nctemp935=iz*df->d[0]+nctemp935;
float nctemp942 = Arbb->a[nctemp932] * df->a[nctemp935];
float nctemp943 = nctemp927 + nctemp942;
Vxxright->a[nctemp911] =nctemp943;
int nctemp952 = i - ix0;
int nctemp947=nctemp952;
nctemp947=iz*P->d[0]+nctemp947;
int nctemp963 = ix - ix0;
int nctemp958=nctemp963;
nctemp958=iz*P->d[0]+nctemp958;
int nctemp978 = i - ix0;
int nctemp973=nctemp978;
nctemp973=iz*L->d[0]+nctemp973;
float nctemp980 = dt * L->a[nctemp973];
int nctemp985=ix;
nctemp985=iz*Vxxright->d[0]+nctemp985;
int nctemp992=ix;
int nctemp1000 = i - ix0;
int nctemp995=nctemp1000;
nctemp995=iz*df->d[0]+nctemp995;
float nctemp1002 = Crbb->a[nctemp992] * df->a[nctemp995];
float nctemp1003 = Vxxright->a[nctemp985] + nctemp1002;
float nctemp1004 = nctemp980 * nctemp1003;
float nctemp1005 = P->a[nctemp958] - nctemp1004;
P->a[nctemp947] =nctemp1005;
}
}
}
}
}
}
int Ac2dFwstepstressx (nctempfloat2 *Vxxleft,nctempfloat2 *Vxxright,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressx<<<NBLOCKS,NTHREADS>>>(Vxxleft,Vxxright,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
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
int nctemp1010=P->d[0];nx =nctemp1010;
int nctemp1018=P->d[1];nz =nctemp1018;
int nctemp1026=Alft->d[0];lpml =nctemp1026;
int nctemp1032=0;
int nctemp1034=nx;
int nctemp1037=0;
int nctemp1039=nz;
int nctemp1030=(nctemp1034-nctemp1032)*(nctemp1039-nctemp1037);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1030;nctempno+=blockDim.x*gridDim.x){
iz=nctemp1037+nctempno/(nctemp1034-nctemp1032);
ix=nctemp1032+nctempno%(nctemp1034-nctemp1032);
{
int nctemp1044=ix;
nctemp1044=iz*P->d[0]+nctemp1044;
int nctemp1051=ix;
nctemp1051=iz*P->d[0]+nctemp1051;
int nctemp1062=ix;
nctemp1062=iz*L->d[0]+nctemp1062;
float nctemp1065 = dt * L->a[nctemp1062];
int nctemp1067=ix;
nctemp1067=iz*df->d[0]+nctemp1067;
float nctemp1070 = nctemp1065 * df->a[nctemp1067];
float nctemp1071 = P->a[nctemp1051] + nctemp1070;
P->a[nctemp1044] =nctemp1071;
}
}
int nctemp1075=4;
int nctemp1078=5;
int nctemp1072 = (Getapplypml->a[nctemp1075] || Getapplypml->a[nctemp1078]);
if(nctemp1072)
{
int nctemp1082=0;
int nctemp1084=nx;
int nctemp1087=0;
int nctemp1089=lpml;
int nctemp1080=(nctemp1084-nctemp1082)*(nctemp1089-nctemp1087);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1080;nctempno+=blockDim.x*gridDim.x){
iz=nctemp1087+nctempno/(nctemp1084-nctemp1082);
ix=nctemp1082+nctempno%(nctemp1084-nctemp1082);
{
int nctemp1092=4;
if(Getapplypml->a[nctemp1092])
{
int nctemp1100 = (iz >= iz0);
int nctemp1097 = (nctemp1100 && iz);
int nctemp1110 = iz0 + nz;
int nctemp1094 = (nctemp1097 < nctemp1110);
if(nctemp1094)
{
int nctemp1114=ix;
nctemp1114=iz*Vzztop->d[0]+nctemp1114;
int nctemp1124=iz;
int nctemp1127=ix;
nctemp1127=iz*Vzztop->d[0]+nctemp1127;
float nctemp1130 = Blft->a[nctemp1124] * Vzztop->a[nctemp1127];
int nctemp1135=iz;
int nctemp1138=ix;
int nctemp1144 = iz - iz0;
nctemp1138=nctemp1144*df->d[0]+nctemp1138;
float nctemp1145 = Alft->a[nctemp1135] * df->a[nctemp1138];
float nctemp1146 = nctemp1130 + nctemp1145;
Vzztop->a[nctemp1114] =nctemp1146;
int nctemp1150=ix;
int nctemp1156 = iz - iz0;
nctemp1150=nctemp1156*P->d[0]+nctemp1150;
int nctemp1161=ix;
int nctemp1167 = iz - iz0;
nctemp1161=nctemp1167*P->d[0]+nctemp1161;
int nctemp1176=ix;
int nctemp1182 = iz - iz0;
nctemp1176=nctemp1182*L->d[0]+nctemp1176;
float nctemp1183 = dt * L->a[nctemp1176];
int nctemp1188=ix;
nctemp1188=iz*Vzztop->d[0]+nctemp1188;
int nctemp1195=iz;
int nctemp1198=ix;
int nctemp1204 = iz - iz0;
nctemp1198=nctemp1204*df->d[0]+nctemp1198;
float nctemp1205 = Clft->a[nctemp1195] * df->a[nctemp1198];
float nctemp1206 = Vzztop->a[nctemp1188] + nctemp1205;
float nctemp1207 = nctemp1183 * nctemp1206;
float nctemp1208 = P->a[nctemp1161] - nctemp1207;
P->a[nctemp1150] =nctemp1208;
}
}
int nctemp1210=5;
if(Getapplypml->a[nctemp1210])
{
int nctemp1223 = iz + nzo;
int nctemp1225 = nctemp1223 - lpml;
i =nctemp1225;
int nctemp1232 = (i >= iz0);
int nctemp1229 = (nctemp1232 && i);
int nctemp1242 = iz0 + nz;
int nctemp1226 = (nctemp1229 < nctemp1242);
if(nctemp1226)
{
int nctemp1246=ix;
nctemp1246=iz*Vzzbot->d[0]+nctemp1246;
int nctemp1256=iz;
int nctemp1259=ix;
nctemp1259=iz*Vzzbot->d[0]+nctemp1259;
float nctemp1262 = Brbb->a[nctemp1256] * Vzzbot->a[nctemp1259];
int nctemp1267=iz;
int nctemp1270=ix;
int nctemp1276 = i - iz0;
nctemp1270=nctemp1276*df->d[0]+nctemp1270;
float nctemp1277 = Arbb->a[nctemp1267] * df->a[nctemp1270];
float nctemp1278 = nctemp1262 + nctemp1277;
Vzzbot->a[nctemp1246] =nctemp1278;
int nctemp1282=ix;
int nctemp1288 = i - iz0;
nctemp1282=nctemp1288*P->d[0]+nctemp1282;
int nctemp1293=ix;
int nctemp1299 = i - iz0;
nctemp1293=nctemp1299*P->d[0]+nctemp1293;
int nctemp1308=ix;
int nctemp1314 = i - iz0;
nctemp1308=nctemp1314*L->d[0]+nctemp1308;
float nctemp1315 = dt * L->a[nctemp1308];
int nctemp1320=ix;
nctemp1320=iz*Vzzbot->d[0]+nctemp1320;
int nctemp1327=iz;
int nctemp1330=ix;
int nctemp1336 = i - iz0;
nctemp1330=nctemp1336*df->d[0]+nctemp1330;
float nctemp1337 = Crbb->a[nctemp1327] * df->a[nctemp1330];
float nctemp1338 = Vzzbot->a[nctemp1320] + nctemp1337;
float nctemp1339 = nctemp1315 * nctemp1338;
float nctemp1340 = P->a[nctemp1293] - nctemp1339;
P->a[nctemp1282] =nctemp1340;
}
}
}
}
}
}
int Ac2dFwstepstressz (nctempfloat2 *Vzztop,nctempfloat2 *Vzzbot,nctempfloat2 *P,nctempfloat2 *L,nctempfloat2 *df,nctempfloat1 *Alft,nctempfloat1 *Blft,nctempfloat1 *Clft,nctempfloat1 *Arbb,nctempfloat1 *Brbb,nctempfloat1 *Crbb,nctempint1 *Getapplypml,int ix0,int iz0,int nxo,int nzo,float dt)
{
  kernel_Ac2dFwstepstressz<<<NBLOCKS,NTHREADS>>>(Vzztop,Vzzbot,P,L,df,Alft,Blft,Clft,Arbb,Brbb,Crbb,Getapplypml,ix0,iz0,nxo,nzo,dt);
GpuError();
return(1);
}
