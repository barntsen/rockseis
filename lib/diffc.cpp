#include "Diffc.h"
#include "diff.h"
#include "balloc.h"
#include <stdio.h>
              
int DiffDxminusc(float *f, float *df, float *coeffs, 
               float dx, int order, int nx, int nz)   
{
  // Add descriptors

  nctempfloat2 *Df;
  Df = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Df->a = f;
  Df->d[0] = nx;
  Df->d[1] = nz;

  nctempfloat2 *Ddf;
  Ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Ddf->a = df;
  Ddf->d[0] = nx;
  Ddf->d[1] = nz;

  nctempfloat1 *Dcoeffs;
  Dcoeffs = (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat2));
  Dcoeffs->a = coeffs;
  Dcoeffs->d[0] = order;

  // Call gpu kernel
  DiffDxminus(Df,Ddf,Dcoeffs,dx,order);

  // Remove descriptors
  BallocDelete(Df);
  BallocDelete(Ddf);
  BallocDelete(Dcoeffs);
  return(1);
}

int DiffDxplusc(float *f, float *df, float *coeffs, 
               float dx, int order, int nx, int nz)   
{
  // Add descriptors

  nctempfloat2 *Df;
  Df = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Df->a = f;
  Df->d[0] = nx;
  Df->d[1] = nz;

  nctempfloat2 *Ddf;
  Ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Ddf->a = df;
  Ddf->d[0] = nx;
  Ddf->d[1] = nz;

  nctempfloat1 *Dcoeffs;
  Dcoeffs = (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  Dcoeffs->a = coeffs;
  Dcoeffs->d[0] = order;

  // Call gpu kernel
  DiffDxplus(Df,Ddf,Dcoeffs,dx,order);

  // Remove descriptors
  BallocDelete(Df);
  BallocDelete(Ddf);
  BallocDelete(Dcoeffs);
  return(1);
}

int DiffDyminusc(float *f, float *df, float *coeffs, 
               float dx, int order, int nx, int nz)   
{
  // Add descriptors

  nctempfloat2 *Df;
  Df = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Df->a = f;
  Df->d[0] = nx;
  Df->d[1] = nz;

  nctempfloat2 *Ddf;
  Ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Ddf->a = df;
  Ddf->d[0] = nx;
  Ddf->d[1] = nz;

  nctempfloat1 *Dcoeffs;
  Dcoeffs = (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  Dcoeffs->a = coeffs;
  Dcoeffs->d[0] = order;

  // Call gpu kernel
  DiffDyminus(Df,Ddf,Dcoeffs,dx,order);

  // Remove descriptors
  BallocDelete(Df);
  BallocDelete(Ddf);
  BallocDelete(Dcoeffs);
  return(1);
}

int DiffDyplusc(float *f, float *df, float *coeffs, 
               float dx, int order, int nx, int nz)   
{
  // Add descriptors

  nctempfloat2 *Df;
  Df = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Df->a = f;
  Df->d[0] = nx;
  Df->d[1] = nz;

  nctempfloat2 *Ddf;
  Ddf = (nctempfloat2 *)BallocNew(1,sizeof(nctempfloat2));
  Ddf->a = df;
  Ddf->d[0] = nx;
  Ddf->d[1] = nz;

  nctempfloat1 *Dcoeffs;
  Dcoeffs = (nctempfloat1 *)BallocNew(1,sizeof(nctempfloat1));
  Dcoeffs->a = coeffs;
  Dcoeffs->d[0] = order;

  // Call gpu kernel
  DiffDyplus(Df,Ddf,Dcoeffs,dx,order);

  // Remove descriptors
  BallocDelete(Df);
  BallocDelete(Ddf);
  BallocDelete(Dcoeffs);
  return(1);
}
