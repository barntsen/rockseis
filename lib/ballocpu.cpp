//
//Balloc contains support functions for memory allocation.
//
#include "balloc.h"
#include <cstdlib>

extern "C" {
#include <stdio.h>
}

//
//BallocNew allocates memory on host and gpu
//
void * BallocNew(int n, int size){

    return(calloc(n,size));
}
//
//BallocDelete free up memory
//
void BallocDelete(void * a){

    free(a);

    return;
}
