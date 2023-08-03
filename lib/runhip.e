/*
%============================================================
\section{Run.c -- The C run time library}
%============================================================
The run time library is written in C.
Most of the routines in the library
are wrappers to unix system calls and math functions.
\begin{verbatim}
*/
#include "hip/hip_runtime.h"
extern "C" {
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>
#include<time.h>
#include<stdlib.h>
#include<string.h>
#include<stdio.h>
}
#define OK   1
#define ERR  0
#define EOF (-1)


#define PERMS 0666
#define MAXRANK 1

typedef struct nctempchar1 { int d[MAXRANK]; char *a;} nctempchar1; 
struct MainArg {nctempchar1 *arg;};
struct nctempMainArg1 {int d[MAXRANK]; struct MainArg *a; };
int Main (struct nctempMainArg1 *MainArgs);
/*
\end{verbatim}
%============================================================
\section{Main -- the main function}
%============================================================
\begin{verbatim}
*/
int main(int argc, char ** argv)
{
  struct nctempMainArg1 *cmlargs;
  int i;
  int rval;

  cmlargs = (struct nctempMainArg1*)malloc(sizeof(struct nctempMainArg1));
  cmlargs->a = (struct MainArg*) malloc(argc*sizeof(struct MainArg));
  cmlargs->d[0] = argc;
  for(i=0; i<argc; i=i+1){
    cmlargs->a[i].arg = (struct nctempchar1*)malloc(sizeof(struct nctempchar1));
    cmlargs->a[i].arg->a = argv[i]; 
    cmlargs->a[i].arg->d[0] = strlen(argv[i])+1;
  }
  rval = Main(cmlargs);
  if(rval == ERR) 
    rval = -1;
  else
    rval = 0;

  return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunMalloc -- allocate memory}
%============================================================
\begin{verbatim}
*/
char * RunMalloc(int nb)
{
    return((char *)malloc(nb));
}
/*
\end{verbatim}
%============================================================
\section{RunFree -- free memory}
%============================================================
\begin{verbatim}
*/
int RunFree(char* p)
{
    free(p);
    return(OK);
}
/*
\end{verbatim}
%============================================================
\section{RunClock -- measure elapsed time}
%============================================================
\begin{verbatim}
*/
float RunClock()
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (float)((double)tp.tv_sec + (double)tp.tv_nsec*1.0e-9) ;
}
/*
\end{verbatim}
%============================================================
\section{RunCreat -- create a file}
%============================================================
\begin{verbatim}
*/
int RunCreate(nctempchar1* name)
{
  int rval;
  
  rval = creat(name->a, PERMS);
  if(rval == -1)
    return(ERR);
  else
    return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunOpen -- open a file}
%============================================================
\begin{verbatim}
*/
int RunOpen(nctempchar1 *name, nctempchar1 *mode)
{
  int rval;
  int flag;

  if(*(mode->a) == 'r')
    flag=O_RDONLY;
  else if(*(mode->a) == 'w')
    flag = O_WRONLY;
  else if(*(mode->a) == 'a')
    flag = O_RDWR | O_APPEND;
  else
    return(ERR);

  rval = open(name->a, flag);
  if(rval == -1)
    return(ERR);
  else
    return (rval);
}
/*
\end{verbatim}
%============================================================
\section{RunClose -- close a file}
%============================================================
\begin{verbatim}
*/
int RunClose(int fd)
{
  int rval;
  
  rval = close(fd);
  if(rval == -1)
    return(ERR);
  else
    return(OK);
}
/*
\end{verbatim}
%============================================================
\section{RunRead -- read from a file}
%============================================================
{\tt RunRead} reads in {\tt lbuff} characters into the
{\tt buffer} array from a file with descriptor {\tt fd}.
The return value is the number of characters actually read.
If an error has occured {\tt ERR} will be returned.
If the end of a file is reached, {\tt EOF} will
be returned.
The {\tt read} routine is a standard UNIX system call.
\begin{verbatim}
*/
int RunRead(int fd, int lbuff, nctempchar1 *buffer)
{
  int rval;
  rval = (int)read(fd, (void *)buffer->a, (size_t)lbuff);
  if(rval == 0)
    rval=EOF;
  else if(rval == -1)
    rval = ERR;
  return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunWrite -- write to a file}
%============================================================
{\tt RunWrite} writes {\tt lbuff} from the {\tt buffer} array
into a file with file descriptor {\tt fd}.
{\tt bufferdesc} is an integer array containg the length
of the {\tt buffer} array.
The return value is the number of characters actually written.
{\tt ERR} is returned whenever an error has occured.
The {\tt write} routine is a standard UNIX system call.
\begin{verbatim}
*/
int RunWrite(int fd, int lbuff, nctempchar1 *buffer)
{
  int rval;
  if(lbuff==0)return(0);
  rval = (int)write(fd, (void *)buffer->a, (size_t)lbuff);
  if(rval == -1)rval=ERR;
  return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunSeek -- Seek a file}
%============================================================
{\tt RunSeek} sets the position of the file pointer to {\tt pos} 
bytes relative to the beginning of the file if flag equals 0, or
realtive to the current position if {\tt flag} equals 1 or relative 
to the end of the file if {\tt flag} equals 2.
\begin{verbatim}
*/
int RunSeek(int fd, int flag, int pos)
{
  int rval;
  int fflag;

  if(flag == 0)
    fflag=SEEK_SET;
  else if(flag == 1)
    fflag=SEEK_CUR;
  else if(flag == 2)
    fflag=SEEK_END;
  else
    return(-1);

  rval = (int)lseek(fd, (off_t) pos, fflag);
  if(rval == -1)return(-1);
  return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunGetenv -- Get environment variable}
%============================================================
\begin{verbatim}
*/
nctempchar1 *RunGetenv(nctempchar1* str)
{
  nctempchar1* rval=(nctempchar1*)malloc(sizeof(nctempchar1));
  rval->a=getenv(str->a); 
  rval->d[0] = strlen(rval->a);
  return(rval);
}
/*
\end{verbatim}
%============================================================
\section{RunExit -- clean up and exit}
%============================================================
{\tt RunExit} will attempt to close all excisting files
and then exit.
*/
int RunExit()
{
  exit(-1);
  return(OK);
}
/*
\end{verbatim}
%==============================================================
\section{GpuNew -- Allocate unified memory on host and gpu}
%==============================================================
\begin{verbatim}
*/
void *GpuNew(int n);
void GpuDelete(void *f);
void GpuErrorCheck(char *s);
void * GpuNew(int n){
  void *f;
  hipError_t cerr;
  cerr = hipMallocManaged(&f, n);
  if(cerr != hipSuccess){
    fprintf(stderr,"GpuAlloc:%s\n ", hipGetErrorString(cerr)) ;
    exit(1);
  }
  return(f);
}
/*
\end{verbatim}
%==============================================================
\section{GpuDelete -- Delete unified memory on host and gpu}
%==============================================================
\begin{verbatim}
*/
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
/*
\end{verbatim}
%==============================================================
\section{GpuError -- Check for gpu errors}
%==============================================================
\begin{verbatim}
*/
void GpuError(){
  hipDeviceSynchronize();
  hipError_t cerr;
  cerr = hipGetLastError();
  if(cerr != hipSuccess){
    fprintf(stderr,"%s\n",hipGetErrorString(cerr));
    exit(1);
  }
}
/*
\end{verbatim}
*/

