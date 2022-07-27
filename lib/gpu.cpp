//
//Gpu contains support functions for cuda
//The run time gpu support is written in c and consists of
//a few calls to the Cuda run time api.
//
//

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
    cudaError_t cerr;
    cerr = cudaMallocManaged(&f, (size_t)n);

    if(cerr != cudaSuccess){
      fprintf(stderr,"GpuAlloc:%s\n ", cudaGetErrorString(cerr)) ;
      exit(1);
    }

    cerr=cudaDeviceSynchronize();
    if(cerr != cudaSuccess){
      fprintf(stderr,"GpuAlloc:%s\n ", cudaGetErrorString(cerr)) ;
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
    cudaError_t cerr;

    cerr=cudaFree(f);
    if(cerr != cudaSuccess){
       fprintf(stderr,"GpuDelete:%s\n ", cudaGetErrorString(cerr)) ;
       exit(1);
    }
    cerr=cudaDeviceSynchronize();
    if(cerr != cudaSuccess){
       fprintf(stderr,"GpuDelete:%s\n ", cudaGetErrorString(cerr)) ;
       exit(1);
    }
}
//
// GpuError checks for gpu errors and sync
//
void GpuError(){
    cudaDeviceSynchronize();
    cudaError_t cerr;
    cerr = cudaGetLastError();
    if(cerr != cudaSuccess){
        fprintf(stderr,"GpuError: %s\n",cudaGetErrorString(cerr));
        exit(1);
    }
}
