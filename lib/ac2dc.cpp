#include "balloc.h"
#include "Ac2dc.h"
#include "ac2d.h"
#include <stdio.h>
              
int Ac2dFwstepvxc(float *Pleft, float *Pright,   
              float *Vx,float *Rx, float *df,    
              float *Alftstag,float *Blftstag,float *Clftstag, 
              float *Arbbstag,float *Brbbstag,float *Crbbstag, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng) 
{
  // Add descriptors

  nctempfloat2 *dPleft;
  dPleft = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dPleft->a = Pleft;
  dPleft->d[0] = lpml;
  dPleft->d[1] = nz;

  nctempfloat2 *dPright;
  dPright = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dPright->a = Pright;
  dPright->d[0] = lpml; 
  dPright->d[1] = nz;

  nctempfloat2 *dVx;
  dVx = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVx->a = Vx;
  dVx->d[0] = nx;
  dVx->d[1] = nz;

  nctempfloat2 *ddf;
  ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  ddf->a = df;
  ddf->d[0] = nx;
  ddf->d[1] = nz;

  nctempfloat2 *dRx;
  dRx = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dRx->a = Rx;
  dRx->d[0] = nx;
  dRx->d[1] = nz;

  nctempfloat1 *dAlftstag;
  dAlftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dAlftstag->a = Alftstag;
  dAlftstag->d[0] = lpml;

  nctempfloat1 *dBlftstag;
  dBlftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBlftstag->a = Blftstag;
  dBlftstag->d[0] = lpml;

  nctempfloat1 *dClftstag;
  dClftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dClftstag->a = Clftstag;
  dClftstag->d[0] = lpml;

  nctempfloat1 *dArbbstag;
  dArbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dArbbstag->a = Arbbstag;
  dArbbstag->d[0] = lpml;

  nctempfloat1 *dBrbbstag;
  dBrbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBrbbstag->a = Brbbstag;
  dBrbbstag->d[0] = lpml;

  nctempfloat1 *dCrbbstag;
  dCrbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dCrbbstag->a = Crbbstag;
  dCrbbstag->d[0] = lpml;

  nctempint1 *dGetapplypml;
  dGetapplypml= (nctempint1 *)BallocNew(1,sizeof(nctempint1));
  dGetapplypml->a = Getapplypml;
  dGetapplypml->d[0] = ng;
  // Call gpu kernel
  Ac2dFwstepvx(dPleft, dPright,   
               dVx,dRx,ddf,    
               dAlftstag,dBlftstag,dClftstag, 
               dArbbstag,dBrbbstag,dCrbbstag, 
               dGetapplypml,ix0,iz0,nxo,nzo,dt); 
  // Call gpu kernel
  Ac2dFwstepvx2(dPleft, dPright,   
               dVx,dRx,ddf,    
               dAlftstag,dBlftstag,dClftstag, 
               dArbbstag,dBrbbstag,dCrbbstag, 
               dGetapplypml,ix0,iz0,nxo,nzo,dt); 

  // Remove descriptors
  BallocDelete(dPleft);
  BallocDelete(dPright);
  BallocDelete(dVx);
  BallocDelete(ddf);
  BallocDelete(dRx);
  BallocDelete(dAlftstag);
  BallocDelete(dBlftstag);
  BallocDelete(dClftstag);
  BallocDelete(dArbbstag);
  BallocDelete(dBrbbstag);
  BallocDelete(dCrbbstag);
  BallocDelete(dGetapplypml);
  return(1);
}

int Ac2dFwstepvzc(float *Ptop, float *Pbot,   
              float *Vz,float *Rz, float *df,    
              float *Alftstag,float *Blftstag,float *Clftstag, 
              float *Arbbstag,float *Brbbstag,float *Crbbstag, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng) 
{
  // Add descriptors

  nctempfloat2 *dPtop;
  dPtop = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dPtop->a = Ptop;
  dPtop->d[0] = nx;
  dPtop->d[1] = lpml;

  nctempfloat2 *dPbot;
  dPbot = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dPbot->a = Pbot;
  dPbot->d[0] = nx;
  dPbot->d[1] = lpml;

  nctempfloat2 *dVz;
  dVz = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVz->a = Vz;
  dVz->d[0] = nx;
  dVz->d[1] = nz;

  nctempfloat2 *ddf;
  ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  ddf->a = df;
  ddf->d[0] = nx;
  ddf->d[1] = nz;

  nctempfloat2 *dRz;
  dRz = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dRz->a = Rz;
  dRz->d[0] = nx;
  dRz->d[1] = nz;

  nctempfloat1 *dAlftstag;
  dAlftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dAlftstag->a = Alftstag;
  dAlftstag->d[0] = lpml;

  nctempfloat1 *dBlftstag;
  dBlftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBlftstag->a = Blftstag;
  dBlftstag->d[0] = lpml;

  nctempfloat1 *dClftstag;
  dClftstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dClftstag->a = Clftstag;
  dClftstag->d[0] = lpml;

  nctempfloat1 *dArbbstag;
  dArbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dArbbstag->a = Arbbstag;
  dArbbstag->d[0] = lpml;

  nctempfloat1 *dBrbbstag;
  dBrbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBrbbstag->a = Brbbstag;
  dBrbbstag->d[0] = lpml;

  nctempfloat1 *dCrbbstag;
  dCrbbstag= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dCrbbstag->a = Crbbstag;
  dCrbbstag->d[0] = lpml;

  nctempint1 *dGetapplypml;
  dGetapplypml= (nctempint1 *)BallocNew(1,sizeof(nctempint1));
  dGetapplypml->a = Getapplypml;
  dGetapplypml->d[0] = ng;


  // Call gpu kernel
  Ac2dFwstepvz(dPtop, dPbot,   
               dVz,dRz,ddf,    
               dAlftstag,dBlftstag,dClftstag, 
               dArbbstag,dBrbbstag,dCrbbstag, 
               dGetapplypml,ix0,iz0,nxo,nzo,dt); 
  // Call gpu kernel
  Ac2dFwstepvz2(dPtop, dPbot,   
               dVz,dRz,ddf,    
               dAlftstag,dBlftstag,dClftstag, 
               dArbbstag,dBrbbstag,dCrbbstag, 
               dGetapplypml,ix0,iz0,nxo,nzo,dt); 
              
  // Remove descriptors
  BallocDelete(dPtop);
  BallocDelete(dPbot);
  BallocDelete(dVz);
  BallocDelete(ddf);
  BallocDelete(dRz);
  BallocDelete(dAlftstag);
  BallocDelete(dBlftstag);
  BallocDelete(dClftstag);
  BallocDelete(dArbbstag);
  BallocDelete(dBrbbstag);
  BallocDelete(dCrbbstag);
  BallocDelete(dGetapplypml);

  return(1);
}

int Ac2dFwstepstressxc(float *Vxxleft, float *Vxxright,   
              float *P,float *L, float *df,    
              float *Alft,float *Blft,float *Clft, 
              float *Arbb,float *Brbb,float *Crbb, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng) 
{
  // Add descriptors

  nctempfloat2 *dVxxleft;
  dVxxleft = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVxxleft->a = Vxxleft;
  dVxxleft->d[0] = lpml;
  dVxxleft->d[1] = nz;

  nctempfloat2 *dVxxright;
  dVxxright = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVxxright->a = Vxxright;
  dVxxright->d[0] = lpml; 
  dVxxright->d[1] = nz;

  nctempfloat2 *dP;
  dP = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dP->a = P;
  dP->d[0] = nx;
  dP->d[1] = nz;

  nctempfloat2 *ddf;
  ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  ddf->a = df;
  ddf->d[0] = nx;
  ddf->d[1] = nz;

  nctempfloat2 *dL;
  dL = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dL->a = L;
  dL->d[0] = nx;
  dL->d[1] = nz;

  nctempfloat1 *dAlft;
  dAlft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dAlft->a = Alft;
  dAlft->d[0] = lpml;

  nctempfloat1 *dBlft;
  dBlft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBlft->a = Blft;
  dBlft->d[0] = lpml;

  nctempfloat1 *dClft;
  dClft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dClft->a = Clft;
  dClft->d[0] = lpml;

  nctempfloat1 *dArbb;
  dArbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dArbb->a = Arbb;
  dArbb->d[0] = lpml;

  nctempfloat1 *dBrbb;
  dBrbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBrbb->a = Brbb;
  dBrbb->d[0] = lpml;

  nctempfloat1 *dCrbb;
  dCrbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dCrbb->a = Crbb;
  dCrbb->d[0] = lpml;

  nctempint1 *dGetapplypml;
  dGetapplypml= (nctempint1 *)BallocNew(1,sizeof(nctempint1));
  dGetapplypml->a = Getapplypml;
  dGetapplypml->d[0] = ng;

  // Call gpu kernel
  Ac2dFwstepstressx(dVxxleft, dVxxright,
              dP,dL,ddf,
              dAlft,dBlft,dClft,
              dArbb,dBrbb,dCrbb,
              dGetapplypml,ix0,iz0,nxo,nzo,dt);
  Ac2dFwstepstressx2(dVxxleft, dVxxright,
              dP,dL,ddf,
              dAlft,dBlft,dClft,
              dArbb,dBrbb,dCrbb,
              dGetapplypml,ix0,iz0,nxo,nzo,dt);

  // Remove descriptors
  BallocDelete(dVxxleft);
  BallocDelete(dVxxright);
  BallocDelete(dP);
  BallocDelete(ddf);
  BallocDelete(dL);
  BallocDelete(dAlft);
  BallocDelete(dBlft);
  BallocDelete(dClft);
  BallocDelete(dArbb);
  BallocDelete(dBrbb);
  BallocDelete(dCrbb);
  BallocDelete(dGetapplypml);
  return(1);
}

int Ac2dFwstepstresszc(float *Vzztop, float *Vzzbot,   
              float *P,float *L, float *df,    
              float *Alft,float *Blft,float *Clft, 
              float *Arbb,float *Brbb,float *Crbb, 
              int   *Getapplypml,int ix0,int iz0,int nxo,int nzo,
              float dt,int nx, int nz,int lpml,int ng) 
{
  // Add descriptors

  nctempfloat2 *dVzztop;
  dVzztop = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVzztop->a = Vzztop;
  dVzztop->d[0] = nx;
  dVzztop->d[1] = lpml;

  nctempfloat2 *dVzzbot;
  dVzzbot = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVzzbot->a = Vzzbot;
  dVzzbot->d[0] = nx; 
  dVzzbot->d[1] = lpml;

  nctempfloat2 *dP;
  dP = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dP->a = P;
  dP->d[0] = nx;
  dP->d[1] = nz;

  nctempfloat2 *ddf;
  ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  ddf->a = df;
  ddf->d[0] = nx;
  ddf->d[1] = nz;

  nctempfloat2 *dL;
  dL = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dL->a = L;
  dL->d[0] = nx;
  dL->d[1] = nz;

  nctempfloat1 *dAlft;
  dAlft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dAlft->a = Alft;
  dAlft->d[0] = lpml;

  nctempfloat1 *dBlft;
  dBlft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBlft->a = Blft;
  dBlft->d[0] = lpml;

  nctempfloat1 *dClft;
  dClft= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dClft->a = Clft;
  dClft->d[0] = lpml;

  nctempfloat1 *dArbb;
  dArbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dArbb->a = Arbb;
  dArbb->d[0] = lpml;

  nctempfloat1 *dBrbb;
  dBrbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dBrbb->a = Brbb;
  dBrbb->d[0] = lpml;

  nctempfloat1 *dCrbb;
  dCrbb= (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  dCrbb->a = Crbb;
  dCrbb->d[0] = lpml;

  nctempint1 *dGetapplypml;
  dGetapplypml= (nctempint1 *)BallocNew(1,sizeof(nctempint1));
  dGetapplypml->a = Getapplypml;
  dGetapplypml->d[0] = ng;

  // Call gpu kernel
  Ac2dFwstepstressz(dVzztop, dVzzbot,
              dP,dL,ddf,
              dAlft,dBlft,dClft,
              dArbb,dBrbb,dCrbb,
              dGetapplypml,ix0,iz0,nxo,nzo,dt);
  Ac2dFwstepstressz2(dVzztop, dVzzbot,
              dP,dL,ddf,
              dAlft,dBlft,dClft,
              dArbb,dBrbb,dCrbb,
              dGetapplypml,ix0,iz0,nxo,nzo,dt);

  // Remove descriptors
  BallocDelete(dVzztop);
  BallocDelete(dVzzbot);
  BallocDelete(dP);
  BallocDelete(ddf);
  BallocDelete(dL);
  BallocDelete(dAlft);
  BallocDelete(dBlft);
  BallocDelete(dClft);
  BallocDelete(dArbb);
  BallocDelete(dBrbb);
  BallocDelete(dCrbb);
  BallocDelete(dGetapplypml);
  return(1);
}

// Xcorr performs cross correlations for gradient computation
int Ac2dXcorrc(float *Vp,  int nx, int nz, int padr, int pads, 
              float *Rho, float *Rx, float *Rz,int nxs, int nzs, 
	      int nxr, int nzr, float *wsp, float *wrx, float *wrz, 
              float *vpgraddata, float *rhograddata, 
	      float dx, float dz, int srcilumset, float * srcilumdata)  
{

  // Add descriptors

  nctempfloat2 *dVp;
  dVp = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dVp->a = Vp;
  dVp->d[0] = nx;
  dVp->d[1] = nz;

  nctempfloat2 *dRho;
  dRho = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dRho->a = Rho;
  dRho->d[0] = nx; 
  dRho->d[1] = nz;

  nctempfloat2 *dRx;
  dRx = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dRx->a = Rx;
  dRx->d[0] = nxr;
  dRx->d[1] = nzr;

  nctempfloat2 *dRz;
  dRz = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dRz->a = Rz;
  dRz->d[0] = nxr;
  dRz->d[1] = nzr;

  nctempfloat2 *dwsp;
  dwsp = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dwsp->a = wsp;
  dwsp->d[0] = nxs;
  dwsp->d[1] = nzs;

  nctempfloat2 *dwrx;
  dwrx = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dwrx->a = wrx;
  dwrx->d[0] = nxr;
  dwrx->d[1] = nzr;

  nctempfloat2 *dwrz;
  dwrz = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dwrz->a = wrz;
  dwrz->d[0] = nxr;
  dwrz->d[1] = nzr;

  nctempfloat2 *dvpgraddata;
  dvpgraddata = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dvpgraddata->a = vpgraddata;
  dvpgraddata->d[0] = nx;
  dvpgraddata->d[1] = nz;

  nctempfloat2 *drhograddata;
  drhograddata = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  drhograddata->a = rhograddata;
  drhograddata->d[0] = nx;
  drhograddata->d[1] = nz;

  nctempfloat2 *dsrcilumdata;
  dsrcilumdata = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  dsrcilumdata->a = srcilumdata;
  dsrcilumdata->d[0] = nx;
  dsrcilumdata->d[1] = nz;

  // Call gpu kernel
   
  Ac2dXcorr(dVp, padr, pads, 
              dRho, dRx,  dRz, 
	      dwsp, dwrx, dwrz, 
              dvpgraddata, drhograddata, 
	      dx, dz, srcilumset, dsrcilumdata);  

  // Remove descriptors
  BallocDelete(dVp);
  BallocDelete(dRho);
  BallocDelete(dRx);
  BallocDelete(dRz);
  BallocDelete(dwsp);
  BallocDelete(dwrx);
  BallocDelete(dwrz);
  BallocDelete(dvpgraddata);
  BallocDelete(drhograddata);
  BallocDelete(dsrcilumdata);

  return(1);
}

// Memcpy copies data
int Ac2dMemcpyc(void *s,  void *t, int n)
{

  // Add descriptors

  nctempchar1 *ds;
  ds = (nctempchar1 *)BallocNew(1,sizeof(nctempchar1));
  ds->a = (char *)s;
  ds->d[0] = n;

  nctempchar1 *dt;
  dt = (nctempchar1 *)BallocNew(1,sizeof(nctempchar1));
  dt->a = (char *)t;
  dt->d[0] = n;

  Ac2dMemcpy(ds,dt);
  
  BallocDelete(ds);
  BallocDelete(dt);

  return(1);
}
