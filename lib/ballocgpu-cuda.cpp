//
//Balloc contains support functions for gpu memory allocation.
//
#include "balloc.h"
#include "gpu.h"

extern "C" {
#include <stdio.h>
}

//
//BallocNew allocates memory on host and gpu
//
void * BallocNew(int n, int size){

    return(GpuNew(n*size));
}
//
//BallocDelete free up memory
//
void BallocDelete(void * a){

    GpuDelete(a);

    return;
}
