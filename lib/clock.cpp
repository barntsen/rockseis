/* Author: B. Arntsen <borge.arntsen@ntnu.no>
 * Date: November 15 2017
 */

// Libraries
#include <ctime>
#include "clock.h"
#include<time.h>
#include<sys/time.h>
#include<sys/types.h>

namespace rockseis {
// =========== CLOCK CLASS =========== //

template<typename T>
Clock<T>::Clock() {
tend=0;
tstart=tend;
telapsed=0;
wtstart=0;
wtelapsed=0;
}

template<typename T>
Clock<T>::~Clock() {
      // Nothing
}
template<typename T>
void Clock<T>::start() {

  // Cpu time
  tstart=((long double)std::clock())/CLOCKS_PER_SEC;

  // Wall clock time
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  wtstart=(T)tp.tv_sec + (T)tp.tv_nsec*1.0e-9;

}

template<typename T>
void Clock<T>::stop() {

  //Cpu clock time
  tend=((long double)std::clock())/CLOCKS_PER_SEC;
  telapsed = telapsed +tend-tstart;

  // Wall clock time
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  wtend=(T)tp.tv_sec + (T)tp.tv_nsec*1.0e-09;
  wtelapsed=wtelapsed+wtend-wtstart;
}

template<typename T>
void Clock<T>::print(){

//std::cout <<" Cpu  time (seconds) : " << telapsed << "\n";
std::cout <<" Wall time (seconds) : " << wtelapsed << "\n";
}


// =========== TEMPLATE CLASS INSTANTIATION =========== //
template class Clock<double>;
template class Clock<float>;

}
