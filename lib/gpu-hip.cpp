//
//Gpu contains support functions for cuda
//The run time gpu support is written in c and consists of
//a few calls to the Cuda run time api.
//
//

#include "hip/hip_runtime.h"
#include "gpu.h"

extern "C" {
#include <stdio.h>
}

//
//GpuNew allocates unified memory on host and gpu
//
void * GpuNew(int n){
    void *f;
    char *a;
    hipError_t cerr;
    cerr = hipMallocManaged(&f, (size_t)n);

    if(cerr != hipSuccess){
      fprintf(stderr,"GpuAlloc:%s\n ", hipGetErrorString(cerr)) ;
      exit(1);
    }

    cerr=hipDeviceSynchronize();
    if(cerr != hipSuccess){
      fprintf(stderr,"GpuAlloc:%s\n ", hipGetErrorString(cerr)) ;
      exit(1);
    }

    //Zero the allocated array
    a = (char*)f;
    for(int i=0; i<n; i++){
      a[i] = 0;
    }

    return(f);
}
//
//GpuDelete deletes unified memory on host and gpu
//
void GpuDelete(void *f){
    hipError_t cerr;

    cerr=hipFree(f);
    if(cerr != hipSuccess){
       fprintf(stderr,"GpuDelete:%s\n ", hipGetErrorString(cerr)) ;
       exit(1);
    }
    cerr=hipDeviceSynchronize();
    if(cerr != hipSuccess){
       fprintf(stderr,"GpuDelete:%s\n ", hipGetErrorString(cerr)) ;
       exit(1);
    }
}
//
// GpuError checks for gpu errors and sync
//
void GpuError(){
    hipDeviceSynchronize();
    hipError_t cerr;
    cerr = hipGetLastError();
    if(cerr != hipSuccess){
        fprintf(stderr,"GpuError: %s\n",hipGetErrorString(cerr));
        exit(1);
    }
}
