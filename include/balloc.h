//
//Balloc contains support functions for memory allocation.
//

//
//BallocNew allocates memory on host and gpu
//
void * BallocNew(int n, int size);

//
//BallocDelete free up memory
//
void BallocDelete(void * a);

