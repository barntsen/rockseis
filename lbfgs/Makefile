# Simple makefile for lbfgs
all: liblbfgs.a liblbfgs.so

lbfgs.o : lbfgs.c

#Dynamic library
liblbfgs.so : lbfgs.c
		gcc -fPIC -c lbfgs.c
		#gcc -shared -W1,-soname,liblbfgs.so -o liblbfgs.so lbfgs.o
		gcc -shared -o liblbfgs.so lbfgs.o

#Static library
liblbfgs.a: lbfgs.o 
	ar -cvq liblbfgs.a lbfgs.o

#Clean
clean: 
	rm *.o *.a *.so
