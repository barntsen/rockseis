//==============================================================================
// Differentiators
//==============================================================================
//include <libe.i>
include "diff.i"

//------------------------------------------------------------------------------
// Dxminus computes the backward derivative in the x-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+(k-1)dx,j]-A[(i-kdx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
int DiffDxminus(float [*,*] A, float [*,*] dA, float [*] w, float dx, int l){
  int nx, ny;
  int i,j,k;
  float sum;

  nx = len(A,0);
  ny = len(A,1);

  //
  // Left border (1 <i < l+1)
  //

  parallel(i=0:l,j=0:ny)
  {
    sum=0.0;
    for(k=1; k<i+1; k=k+1){
      sum = -w[k-1]*A[i-k,j] + sum; 
    }
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*A[i+(k-1),j] +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Outside border area 
  //
  parallel(i=l:nx-l,j=0:ny)
  {
    sum=0.0;
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*(-A[i-k,j]+A[i+(k-1),j]) +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Right border 
  //
  parallel(i=nx-l:nx,j=0:ny)
  {
    sum = 0.0;
    for(k=1; k<l+1; k=k+1){
      sum = -w[k-1]*A[i-k,j] + sum;
    }

    for(k=1; k<(nx-i+1); k=k+1){
      sum = w[k-1]*A[i+(k-1),j] + sum;
    }
    dA[i,j] = sum/dx;
  }
}
//------------------------------------------------------------------------------
// Dxplus computes the forward derivative in the x-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+kdx,j]-A[(i-(k-1)dx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
int DiffDxplus(float [*,*] A, float [*,*] dA, float [*] w, float dx, int l){
  int nx, ny;
  int i,j,k;
  float sum;

  nx = len(A,0);
  ny = len(A,1);

  //
  // Left border (1 <i < l+1)
  //

  // Left border

  parallel(i=0:l,j=0:ny)
  {
    sum=0.0;
    for(k=1; k<i+2; k=k+1){
      sum = -w[k-1]*A[i-(k-1),j] + sum; 
    }
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*A[i+k,j] +sum; 
    }
    dA[i,j] = sum/dx;
  } 
  //
  // Between left and right border
  //
  parallel(i=l:nx-l,j=0:ny)
  {
    sum=0.0;
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*(-A[i-(k-1),j]+A[i+k,j]) +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Right border 
  //
  parallel(i=nx-l:nx,j=0:ny)
  {
    sum = 0.0;
    for(k=1; k<l+1; k=k+1){
      sum = -w[k-1]*A[i-(k-1),j] + sum;
    }

    for(k=1; k<nx-i; k=k+1){
      sum = w[k-1]*A[i+k,j] + sum;
    }
    dA[i,j] = sum/dx;
  }
}
//------------------------------------------------------------------------------
// Dyminus computes the backward derivative in the y-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i,j+(k-1)dx]-A[i,j-kdx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
int DiffDyminus(float [*,*] A, float [*,*] dA, float [*] w, float dx, int l){
  int nx, ny;
  int i,j,k;
  float sum;

  nx = len(A,0);
  ny = len(A,1);

  //
  // Top border (1 <i < l+1)
  //

  // Left border 

  parallel(i=0:nx,j=0:l)
  {
    sum=0.0;
    for(k=1; k<j+1; k=k+1){
      sum = -w[k-1]*A[i,j-k] + sum; 
    }
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*A[i,j+(k-1)] +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Outside border area 
  //
  parallel(i=0:nx,j=l:ny-l)
  {
    sum=0.0;
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*(-A[i,j-k]+A[i,j+(k-1)]) +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Right border 
  //
  parallel(i=0:nx,j=ny-l:ny)
  {
    sum = 0.0;
    for(k=1; k<l+1; k=k+1){
      sum = -w[k-1]*A[i,j-k] + sum;
   }

    for(k=1; k<(ny-j+1); k=k+1){
      sum = w[k-1]*A[i,j+(k-1)] + sum;
    }

    dA[i,j] = sum/dx;
  }
}
//------------------------------------------------------------------------------
// Dyplus computes the forward derivative in the x-direction
//
// Arguments:
//  Diff: Diff object 
//  float  A       : Input 2D array
//  float dx       : Sampling interval
//  float dA       : Output array 
//
//  The output array, dA, contains the derivative for each point computed
//  as:
//  dA[i,j] = (1/dx) sum_{k=1}^l w[k](A[i+kdx,j]-A[(i-(k-1)dx,j]
//
//  w[k] are weights and l is the length of the differentiator.
//  (see DiffNew for the definitions of these)
//------------------------------------------------------------------------------
int DiffDyplus(float [*,*] A, float [*,*] dA, float [*] w, float dx, int l){
  int nx, ny;
  int i,j,k;
  float sum;

  nx = len(A,0);
  ny = len(A,1);

  //
  // Left border (1 <i < l+1)
  //

  // Left border
  parallel(i=0:nx,j=0:l)
  {
    sum=0.0;
    for(k=1; k<j+2; k=k+1){
      sum = -w[k-1]*A[i,j-(k-1)] + sum; 
    }
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*A[i,j+k] +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Between left and right border
  //
  parallel(i=0:nx,j=l:ny-l)
  {
    sum=0.0;
    for(k=1; k<l+1; k=k+1){
      sum = w[k-1]*(-A[i,j-(k-1)]+A[i,j+k]) +sum; 
    }
    dA[i,j] = sum/dx;
  } 

  //
  // Right border 
  //
  parallel(i=0:nx,j=ny-l:ny)
  {
    sum = 0.0;
    for(k=1; k<l+1; k=k+1){
      sum = -w[k-1]*A[i,j-(k-1)] + sum;
    }

    for(k=1; k<ny-j; k=k+1){
      sum = w[k-1]*A[i,j+k] + sum;
    }
    dA[i,j] = sum/dx;
  }
}
