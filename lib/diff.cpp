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
struct diff {int l;
int lmax;
nctempfloat2 *coeffs;
nctempfloat1 *w;
};
struct nctempdiff1 {int d[1]; struct diff *a; } ;
struct nctempdiff2 {int d[2]; struct diff *a; } ;
struct nctempdiff3 {int d[3]; struct diff *a; } ;
struct nctempdiff4 {int d[4]; struct diff *a; } ;
struct diff* DiffNew (int l);
int DiffDxminus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
int DiffDyminus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
int DiffDxplus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
int DiffDyplus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
__global__ void kernel_DiffDxminus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
__global__ void kernel_DiffDxminus (float *_A, float *_dA, float *_w,float dx,int l, int _nx, int _ny)
{
int nx;
int ny;
int i;
int j;
int k;
float sum;

nctempfloat2 A_;
nctempfloat2 *A; 
A_.d[0] = _nx;
A_.d[1] = _ny;
A_.a = _A;
A= &A_;

nctempfloat2 dA_;
nctempfloat2 *dA; 
dA_.d[0] = _nx;
dA_.d[1] = _ny;
dA_.a = _dA;
dA= &dA_;

nctempfloat1 w_;
nctempfloat1 *w; 
w_.d[0] = l;
w_.a = _w;
w= &w_;

int nctemp5=A->d[0];nx =nctemp5;
int nctemp13=A->d[1];ny =nctemp13;
int nctemp19=0;
int nctemp21=l;
int nctemp24=0;
int nctemp26=ny;
int nctemp17=(nctemp21-nctemp19)*(nctemp26-nctemp24);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp17;nctempno+=blockDim.x*gridDim.x){
j=nctemp24+nctempno/(nctemp21-nctemp19);
i=nctemp19+nctempno%(nctemp21-nctemp19);
{
sum =0.0;
k =1;
int nctemp44 = i + 1;
int nctemp36 = (k < nctemp44);
while(nctemp36){
{
int nctemp60 = k - 1;
int nctemp55=nctemp60;
float nctemp54= -w->a[nctemp55];
int nctemp67 = i - k;
int nctemp62=nctemp67;
nctemp62=j*A->d[0]+nctemp62;
float nctemp69 = nctemp54 * A->a[nctemp62];
float nctemp71 = nctemp69 + sum;
sum =nctemp71;
}
int nctemp80 = k + 1;
k =nctemp80;
int nctemp89 = i + 1;
int nctemp81 = (k < nctemp89);
nctemp36=nctemp81;
}
k =1;
int nctemp102 = l + 1;
int nctemp94 = (k < nctemp102);
while(nctemp94){
{
int nctemp118 = k - 1;
int nctemp113=nctemp118;
int nctemp129 = k - 1;
int nctemp130 = i + nctemp129;
int nctemp120=nctemp130;
nctemp120=j*A->d[0]+nctemp120;
float nctemp132 = w->a[nctemp113] * A->a[nctemp120];
float nctemp134 = nctemp132 + sum;
sum =nctemp134;
}
int nctemp143 = k + 1;
k =nctemp143;
int nctemp152 = l + 1;
int nctemp144 = (k < nctemp152);
nctemp94=nctemp144;
}
int nctemp156=i;
nctemp156=j*dA->d[0]+nctemp156;
float nctemp164 = sum / dx;
dA->a[nctemp156] =nctemp164;
}
}
int nctemp167=l;
int nctemp174 = nx - l;
int nctemp169=nctemp174;
int nctemp176=0;
int nctemp178=ny;
int nctemp165=(nctemp169-nctemp167)*(nctemp178-nctemp176);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp165;nctempno+=blockDim.x*gridDim.x){
j=nctemp176+nctempno/(nctemp169-nctemp167);
i=nctemp167+nctempno%(nctemp169-nctemp167);
{
sum =0.0;
k =1;
int nctemp196 = l + 1;
int nctemp188 = (k < nctemp196);
while(nctemp188){
{
int nctemp212 = k - 1;
int nctemp207=nctemp212;
int nctemp222 = i - k;
int nctemp217=nctemp222;
nctemp217=j*A->d[0]+nctemp217;
float nctemp216= -A->a[nctemp217];
int nctemp234 = k - 1;
int nctemp235 = i + nctemp234;
int nctemp225=nctemp235;
nctemp225=j*A->d[0]+nctemp225;
float nctemp237 = nctemp216 + A->a[nctemp225];
float nctemp238 = w->a[nctemp207] * nctemp237;
float nctemp240 = nctemp238 + sum;
sum =nctemp240;
}
int nctemp249 = k + 1;
k =nctemp249;
int nctemp258 = l + 1;
int nctemp250 = (k < nctemp258);
nctemp188=nctemp250;
}
int nctemp262=i;
nctemp262=j*dA->d[0]+nctemp262;
float nctemp270 = sum / dx;
dA->a[nctemp262] =nctemp270;
}
}
int nctemp279 = nx - l;
int nctemp273=nctemp279;
int nctemp280=nx;
int nctemp283=0;
int nctemp285=ny;
int nctemp271=(nctemp280-nctemp273)*(nctemp285-nctemp283);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp271;nctempno+=blockDim.x*gridDim.x){
j=nctemp283+nctempno/(nctemp280-nctemp273);
i=nctemp273+nctempno%(nctemp280-nctemp273);
{
sum =0.0;
k =1;
int nctemp303 = l + 1;
int nctemp295 = (k < nctemp303);
while(nctemp295){
{
int nctemp319 = k - 1;
int nctemp314=nctemp319;
float nctemp313= -w->a[nctemp314];
int nctemp326 = i - k;
int nctemp321=nctemp326;
nctemp321=j*A->d[0]+nctemp321;
float nctemp328 = nctemp313 * A->a[nctemp321];
float nctemp330 = nctemp328 + sum;
sum =nctemp330;
}
int nctemp339 = k + 1;
k =nctemp339;
int nctemp348 = l + 1;
int nctemp340 = (k < nctemp348);
nctemp295=nctemp340;
}
k =1;
int nctemp364 = nx - i;
int nctemp366 = nctemp364 + 1;
int nctemp353 = (k < nctemp366);
while(nctemp353){
{
int nctemp382 = k - 1;
int nctemp377=nctemp382;
int nctemp393 = k - 1;
int nctemp394 = i + nctemp393;
int nctemp384=nctemp394;
nctemp384=j*A->d[0]+nctemp384;
float nctemp396 = w->a[nctemp377] * A->a[nctemp384];
float nctemp398 = nctemp396 + sum;
sum =nctemp398;
}
int nctemp407 = k + 1;
k =nctemp407;
int nctemp419 = nx - i;
int nctemp421 = nctemp419 + 1;
int nctemp408 = (k < nctemp421);
nctemp353=nctemp408;
}
int nctemp425=i;
nctemp425=j*dA->d[0]+nctemp425;
float nctemp433 = sum / dx;
dA->a[nctemp425] =nctemp433;
}
}
}
int DiffDxminus (float *A, float *dA, float *w,float dx,int l,int nx, int ny)
{
  kernel_DiffDxminus<<<NBLOCKS,NTHREADS>>>(A,dA,w,dx,l,nx,ny);
GpuError();
return(1);
}

__global__ void kernel_DiffDxplus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
__global__ void kernel_DiffDxplus (float *_A, float *_dA, float *_w,float dx,int l, int _nx, int _ny)
{
int nx;
int ny;
int i;
int j;
int k;
float sum;

nctempfloat2 A_;
nctempfloat2 *A; 
A_.d[0] = _nx;
A_.d[1] = _ny;
A_.a = _A;
A= &A_;

nctempfloat2 dA_;
nctempfloat2 *dA; 
dA_.d[0] = _nx;
dA_.d[1] = _ny;
dA_.a = _dA;
dA= &dA_;

nctempfloat1 w_;
nctempfloat1 *w; 
w_.d[0] = l;
w_.a = _w;
w= &w_;

int nctemp438=A->d[0];nx =nctemp438;
int nctemp446=A->d[1];ny =nctemp446;
int nctemp452=0;
int nctemp454=l;
int nctemp457=0;
int nctemp459=ny;
int nctemp450=(nctemp454-nctemp452)*(nctemp459-nctemp457);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp450;nctempno+=blockDim.x*gridDim.x){
j=nctemp457+nctempno/(nctemp454-nctemp452);
i=nctemp452+nctempno%(nctemp454-nctemp452);
{
sum =0.0;
k =1;
int nctemp477 = i + 2;
int nctemp469 = (k < nctemp477);
while(nctemp469){
{
int nctemp493 = k - 1;
int nctemp488=nctemp493;
float nctemp487= -w->a[nctemp488];
int nctemp504 = k - 1;
int nctemp505 = i - nctemp504;
int nctemp495=nctemp505;
nctemp495=j*A->d[0]+nctemp495;
float nctemp507 = nctemp487 * A->a[nctemp495];
float nctemp509 = nctemp507 + sum;
sum =nctemp509;
}
int nctemp518 = k + 1;
k =nctemp518;
int nctemp527 = i + 2;
int nctemp519 = (k < nctemp527);
nctemp469=nctemp519;
}
k =1;
int nctemp540 = l + 1;
int nctemp532 = (k < nctemp540);
while(nctemp532){
{
int nctemp556 = k - 1;
int nctemp551=nctemp556;
int nctemp563 = i + k;
int nctemp558=nctemp563;
nctemp558=j*A->d[0]+nctemp558;
float nctemp565 = w->a[nctemp551] * A->a[nctemp558];
float nctemp567 = nctemp565 + sum;
sum =nctemp567;
}
int nctemp576 = k + 1;
k =nctemp576;
int nctemp585 = l + 1;
int nctemp577 = (k < nctemp585);
nctemp532=nctemp577;
}
int nctemp589=i;
nctemp589=j*dA->d[0]+nctemp589;
float nctemp597 = sum / dx;
dA->a[nctemp589] =nctemp597;
}
}
int nctemp600=l;
int nctemp607 = nx - l;
int nctemp602=nctemp607;
int nctemp609=0;
int nctemp611=ny;
int nctemp598=(nctemp602-nctemp600)*(nctemp611-nctemp609);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp598;nctempno+=blockDim.x*gridDim.x){
j=nctemp609+nctempno/(nctemp602-nctemp600);
i=nctemp600+nctempno%(nctemp602-nctemp600);
{
sum =0.0;
k =1;
int nctemp629 = l + 1;
int nctemp621 = (k < nctemp629);
while(nctemp621){
{
int nctemp645 = k - 1;
int nctemp640=nctemp645;
int nctemp659 = k - 1;
int nctemp660 = i - nctemp659;
int nctemp650=nctemp660;
nctemp650=j*A->d[0]+nctemp650;
float nctemp649= -A->a[nctemp650];
int nctemp668 = i + k;
int nctemp663=nctemp668;
nctemp663=j*A->d[0]+nctemp663;
float nctemp670 = nctemp649 + A->a[nctemp663];
float nctemp671 = w->a[nctemp640] * nctemp670;
float nctemp673 = nctemp671 + sum;
sum =nctemp673;
}
int nctemp682 = k + 1;
k =nctemp682;
int nctemp691 = l + 1;
int nctemp683 = (k < nctemp691);
nctemp621=nctemp683;
}
int nctemp695=i;
nctemp695=j*dA->d[0]+nctemp695;
float nctemp703 = sum / dx;
dA->a[nctemp695] =nctemp703;
}
}
int nctemp712 = nx - l;
int nctemp706=nctemp712;
int nctemp713=nx;
int nctemp716=0;
int nctemp718=ny;
int nctemp704=(nctemp713-nctemp706)*(nctemp718-nctemp716);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp704;nctempno+=blockDim.x*gridDim.x){
j=nctemp716+nctempno/(nctemp713-nctemp706);
i=nctemp706+nctempno%(nctemp713-nctemp706);
{
sum =0.0;
k =1;
int nctemp736 = l + 1;
int nctemp728 = (k < nctemp736);
while(nctemp728){
{
int nctemp752 = k - 1;
int nctemp747=nctemp752;
float nctemp746= -w->a[nctemp747];
int nctemp763 = k - 1;
int nctemp764 = i - nctemp763;
int nctemp754=nctemp764;
nctemp754=j*A->d[0]+nctemp754;
float nctemp766 = nctemp746 * A->a[nctemp754];
float nctemp768 = nctemp766 + sum;
sum =nctemp768;
}
int nctemp777 = k + 1;
k =nctemp777;
int nctemp786 = l + 1;
int nctemp778 = (k < nctemp786);
nctemp728=nctemp778;
}
k =1;
int nctemp799 = nx - i;
int nctemp791 = (k < nctemp799);
while(nctemp791){
{
int nctemp815 = k - 1;
int nctemp810=nctemp815;
int nctemp822 = i + k;
int nctemp817=nctemp822;
nctemp817=j*A->d[0]+nctemp817;
float nctemp824 = w->a[nctemp810] * A->a[nctemp817];
float nctemp826 = nctemp824 + sum;
sum =nctemp826;
}
int nctemp835 = k + 1;
k =nctemp835;
int nctemp844 = nx - i;
int nctemp836 = (k < nctemp844);
nctemp791=nctemp836;
}
int nctemp848=i;
nctemp848=j*dA->d[0]+nctemp848;
float nctemp856 = sum / dx;
dA->a[nctemp848] =nctemp856;
}
}
}
int DiffDxplus (float *A, float *dA, float *w,float dx,int l, int nx, int ny)
{
  kernel_DiffDxplus<<<NBLOCKS,NTHREADS>>>(A,dA,w,dx,l,nx,ny);
GpuError();
return(1);
}
__global__ void kernel_DiffDyminus (float *_A, float *_dA, float *_w,float dx,int l, int nx, int ny);
__global__ void kernel_DiffDyminus (float *_A, float *_dA, float *_w,float dx,int l, int _nx, int _ny)
{
int nx;
int ny;
int i;
int j;
int k;

nctempfloat2 A_;
nctempfloat2 *A; 
A_.d[0] = _nx;
A_.d[1] = _ny;
A_.a = _A;
A= &A_;

nctempfloat2 dA_;
nctempfloat2 *dA; 
dA_.d[0] = _nx;
dA_.d[1] = _ny;
dA_.a = _dA;
dA= &dA_;

nctempfloat1 w_;
nctempfloat1 *w; 
w_.d[0] = l;
w_.a = _w;
w= &w_;
float sum;
int nctemp861=A->d[0];nx =nctemp861;
int nctemp869=A->d[1];ny =nctemp869;
int nctemp875=0;
int nctemp877=nx;
int nctemp880=0;
int nctemp882=l;
int nctemp873=(nctemp877-nctemp875)*(nctemp882-nctemp880);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp873;nctempno+=blockDim.x*gridDim.x){
j=nctemp880+nctempno/(nctemp877-nctemp875);
i=nctemp875+nctempno%(nctemp877-nctemp875);
{
sum =0.0;
k =1;
int nctemp900 = j + 1;
int nctemp892 = (k < nctemp900);
while(nctemp892){
{
int nctemp916 = k - 1;
int nctemp911=nctemp916;
float nctemp910= -w->a[nctemp911];
int nctemp918=i;
int nctemp924 = j - k;
nctemp918=nctemp924*A->d[0]+nctemp918;
float nctemp925 = nctemp910 * A->a[nctemp918];
float nctemp927 = nctemp925 + sum;
sum =nctemp927;
}
int nctemp936 = k + 1;
k =nctemp936;
int nctemp945 = j + 1;
int nctemp937 = (k < nctemp945);
nctemp892=nctemp937;
}
k =1;
int nctemp958 = l + 1;
int nctemp950 = (k < nctemp958);
while(nctemp950){
{
int nctemp974 = k - 1;
int nctemp969=nctemp974;
int nctemp976=i;
int nctemp986 = k - 1;
int nctemp987 = j + nctemp986;
nctemp976=nctemp987*A->d[0]+nctemp976;
float nctemp988 = w->a[nctemp969] * A->a[nctemp976];
float nctemp990 = nctemp988 + sum;
sum =nctemp990;
}
int nctemp999 = k + 1;
k =nctemp999;
int nctemp1008 = l + 1;
int nctemp1000 = (k < nctemp1008);
nctemp950=nctemp1000;
}
int nctemp1012=i;
nctemp1012=j*dA->d[0]+nctemp1012;
float nctemp1020 = sum / dx;
dA->a[nctemp1012] =nctemp1020;
}
}
int nctemp1023=0;
int nctemp1025=nx;
int nctemp1028=l;
int nctemp1035 = ny - l;
int nctemp1030=nctemp1035;
int nctemp1021=(nctemp1025-nctemp1023)*(nctemp1030-nctemp1028);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1021;nctempno+=blockDim.x*gridDim.x){
j=nctemp1028+nctempno/(nctemp1025-nctemp1023);
i=nctemp1023+nctempno%(nctemp1025-nctemp1023);
{
sum =0.0;
k =1;
int nctemp1052 = l + 1;
int nctemp1044 = (k < nctemp1052);
while(nctemp1044){
{
int nctemp1068 = k - 1;
int nctemp1063=nctemp1068;
int nctemp1073=i;
int nctemp1079 = j - k;
nctemp1073=nctemp1079*A->d[0]+nctemp1073;
float nctemp1072= -A->a[nctemp1073];
int nctemp1081=i;
int nctemp1091 = k - 1;
int nctemp1092 = j + nctemp1091;
nctemp1081=nctemp1092*A->d[0]+nctemp1081;
float nctemp1093 = nctemp1072 + A->a[nctemp1081];
float nctemp1094 = w->a[nctemp1063] * nctemp1093;
float nctemp1096 = nctemp1094 + sum;
sum =nctemp1096;
}
int nctemp1105 = k + 1;
k =nctemp1105;
int nctemp1114 = l + 1;
int nctemp1106 = (k < nctemp1114);
nctemp1044=nctemp1106;
}
int nctemp1118=i;
nctemp1118=j*dA->d[0]+nctemp1118;
float nctemp1126 = sum / dx;
dA->a[nctemp1118] =nctemp1126;
}
}
int nctemp1129=0;
int nctemp1131=nx;
int nctemp1140 = ny - l;
int nctemp1134=nctemp1140;
int nctemp1141=ny;
int nctemp1127=(nctemp1131-nctemp1129)*(nctemp1141-nctemp1134);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1127;nctempno+=blockDim.x*gridDim.x){
j=nctemp1134+nctempno/(nctemp1131-nctemp1129);
i=nctemp1129+nctempno%(nctemp1131-nctemp1129);
{
sum =0.0;
k =1;
int nctemp1159 = l + 1;
int nctemp1151 = (k < nctemp1159);
while(nctemp1151){
{
int nctemp1175 = k - 1;
int nctemp1170=nctemp1175;
float nctemp1169= -w->a[nctemp1170];
int nctemp1177=i;
int nctemp1183 = j - k;
nctemp1177=nctemp1183*A->d[0]+nctemp1177;
float nctemp1184 = nctemp1169 * A->a[nctemp1177];
float nctemp1186 = nctemp1184 + sum;
sum =nctemp1186;
}
int nctemp1195 = k + 1;
k =nctemp1195;
int nctemp1204 = l + 1;
int nctemp1196 = (k < nctemp1204);
nctemp1151=nctemp1196;
}
k =1;
int nctemp1220 = ny - j;
int nctemp1222 = nctemp1220 + 1;
int nctemp1209 = (k < nctemp1222);
while(nctemp1209){
{
int nctemp1238 = k - 1;
int nctemp1233=nctemp1238;
int nctemp1240=i;
int nctemp1250 = k - 1;
int nctemp1251 = j + nctemp1250;
nctemp1240=nctemp1251*A->d[0]+nctemp1240;
float nctemp1252 = w->a[nctemp1233] * A->a[nctemp1240];
float nctemp1254 = nctemp1252 + sum;
sum =nctemp1254;
}
int nctemp1263 = k + 1;
k =nctemp1263;
int nctemp1275 = ny - j;
int nctemp1277 = nctemp1275 + 1;
int nctemp1264 = (k < nctemp1277);
nctemp1209=nctemp1264;
}
int nctemp1281=i;
nctemp1281=j*dA->d[0]+nctemp1281;
float nctemp1289 = sum / dx;
dA->a[nctemp1281] =nctemp1289;
}
}
}
int DiffDyminus (float *A, float *dA, float *w, float dx,int l, int nx, int ny)
{
  kernel_DiffDyminus<<<NBLOCKS,NTHREADS>>>(A,dA,w,dx,l,nx,ny);
GpuError();
return(1);
}

__global__ void kernel_DiffDyplus (float *_A, float *_dA, float *_w,float dx,int l,int _nx,int _ny);
__global__ void kernel_DiffDyplus (float *_A, float *_dA, float *_w,float dx,int l,int _nx,int _ny)
{
int nx;
int ny;
int i;
int j;
int k;
float sum;

nctempfloat2 A_;
nctempfloat2 *A; 
A_.d[0] = _nx;
A_.d[1] = _ny;
A_.a = _A;
A= &A_;

nctempfloat2 dA_;
nctempfloat2 *dA; 
dA_.d[0] = _nx;
dA_.d[1] = _ny;
dA_.a = _dA;
dA= &dA_;

nctempfloat1 w_;
nctempfloat1 *w; 
w_.d[0] = l;
w_.a = _w;
w= &w_;

int nctemp1294=A->d[0];nx =nctemp1294;
int nctemp1302=A->d[1];ny =nctemp1302;
int nctemp1308=0;
int nctemp1310=nx;
int nctemp1313=0;
int nctemp1315=l;
int nctemp1306=(nctemp1310-nctemp1308)*(nctemp1315-nctemp1313);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1306;nctempno+=blockDim.x*gridDim.x){
j=nctemp1313+nctempno/(nctemp1310-nctemp1308);
i=nctemp1308+nctempno%(nctemp1310-nctemp1308);
{
sum =0.0;
k =1;
int nctemp1333 = j + 2;
int nctemp1325 = (k < nctemp1333);
while(nctemp1325){
{
int nctemp1349 = k - 1;
int nctemp1344=nctemp1349;
float nctemp1343= -w->a[nctemp1344];
int nctemp1351=i;
int nctemp1361 = k - 1;
int nctemp1362 = j - nctemp1361;
nctemp1351=nctemp1362*A->d[0]+nctemp1351;
float nctemp1363 = nctemp1343 * A->a[nctemp1351];
float nctemp1365 = nctemp1363 + sum;
sum =nctemp1365;
}
int nctemp1374 = k + 1;
k =nctemp1374;
int nctemp1383 = j + 2;
int nctemp1375 = (k < nctemp1383);
nctemp1325=nctemp1375;
}
k =1;
int nctemp1396 = l + 1;
int nctemp1388 = (k < nctemp1396);
while(nctemp1388){
{
int nctemp1412 = k - 1;
int nctemp1407=nctemp1412;
int nctemp1414=i;
int nctemp1420 = j + k;
nctemp1414=nctemp1420*A->d[0]+nctemp1414;
float nctemp1421 = w->a[nctemp1407] * A->a[nctemp1414];
float nctemp1423 = nctemp1421 + sum;
sum =nctemp1423;
}
int nctemp1432 = k + 1;
k =nctemp1432;
int nctemp1441 = l + 1;
int nctemp1433 = (k < nctemp1441);
nctemp1388=nctemp1433;
}
int nctemp1445=i;
nctemp1445=j*dA->d[0]+nctemp1445;
float nctemp1453 = sum / dx;
dA->a[nctemp1445] =nctemp1453;
}
}
int nctemp1456=0;
int nctemp1458=nx;
int nctemp1461=l;
int nctemp1468 = ny - l;
int nctemp1463=nctemp1468;
int nctemp1454=(nctemp1458-nctemp1456)*(nctemp1463-nctemp1461);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1454;nctempno+=blockDim.x*gridDim.x){
j=nctemp1461+nctempno/(nctemp1458-nctemp1456);
i=nctemp1456+nctempno%(nctemp1458-nctemp1456);
{
sum =0.0;
k =1;
int nctemp1485 = l + 1;
int nctemp1477 = (k < nctemp1485);
while(nctemp1477){
{
int nctemp1501 = k - 1;
int nctemp1496=nctemp1501;
int nctemp1506=i;
int nctemp1516 = k - 1;
int nctemp1517 = j - nctemp1516;
nctemp1506=nctemp1517*A->d[0]+nctemp1506;
float nctemp1505= -A->a[nctemp1506];
int nctemp1519=i;
int nctemp1525 = j + k;
nctemp1519=nctemp1525*A->d[0]+nctemp1519;
float nctemp1526 = nctemp1505 + A->a[nctemp1519];
float nctemp1527 = w->a[nctemp1496] * nctemp1526;
float nctemp1529 = nctemp1527 + sum;
sum =nctemp1529;
}
int nctemp1538 = k + 1;
k =nctemp1538;
int nctemp1547 = l + 1;
int nctemp1539 = (k < nctemp1547);
nctemp1477=nctemp1539;
}
int nctemp1551=i;
nctemp1551=j*dA->d[0]+nctemp1551;
float nctemp1559 = sum / dx;
dA->a[nctemp1551] =nctemp1559;
}
}
int nctemp1562=0;
int nctemp1564=nx;
int nctemp1573 = ny - l;
int nctemp1567=nctemp1573;
int nctemp1574=ny;
int nctemp1560=(nctemp1564-nctemp1562)*(nctemp1574-nctemp1567);
for(int nctempno=blockIdx.x*blockDim.x + threadIdx.x; nctempno<nctemp1560;nctempno+=blockDim.x*gridDim.x){
j=nctemp1567+nctempno/(nctemp1564-nctemp1562);
i=nctemp1562+nctempno%(nctemp1564-nctemp1562);
{
sum =0.0;
k =1;
int nctemp1592 = l + 1;
int nctemp1584 = (k < nctemp1592);
while(nctemp1584){
{
int nctemp1608 = k - 1;
int nctemp1603=nctemp1608;
float nctemp1602= -w->a[nctemp1603];
int nctemp1610=i;
int nctemp1620 = k - 1;
int nctemp1621 = j - nctemp1620;
nctemp1610=nctemp1621*A->d[0]+nctemp1610;
float nctemp1622 = nctemp1602 * A->a[nctemp1610];
float nctemp1624 = nctemp1622 + sum;
sum =nctemp1624;
}
int nctemp1633 = k + 1;
k =nctemp1633;
int nctemp1642 = l + 1;
int nctemp1634 = (k < nctemp1642);
nctemp1584=nctemp1634;
}
k =1;
int nctemp1655 = ny - j;
int nctemp1647 = (k < nctemp1655);
while(nctemp1647){
{
int nctemp1671 = k - 1;
int nctemp1666=nctemp1671;
int nctemp1673=i;
int nctemp1679 = j + k;
nctemp1673=nctemp1679*A->d[0]+nctemp1673;
float nctemp1680 = w->a[nctemp1666] * A->a[nctemp1673];
float nctemp1682 = nctemp1680 + sum;
sum =nctemp1682;
}
int nctemp1691 = k + 1;
k =nctemp1691;
int nctemp1700 = ny - j;
int nctemp1692 = (k < nctemp1700);
nctemp1647=nctemp1692;
}
int nctemp1704=i;
nctemp1704=j*dA->d[0]+nctemp1704;
float nctemp1712 = sum / dx;
dA->a[nctemp1704] =nctemp1712;
}
}
}
int DiffDyplus (float *A, float *dA, float *w,float dx, int l, int nx, int ny)
{
  kernel_DiffDyplus<<<NBLOCKS,NTHREADS>>>(A,dA,w,dx,l,nx,ny);
GpuError();
return(1);
}
