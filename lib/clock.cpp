/* Author: B. Arntsen <borge.arntsen@ntnu.no>
 * Date: November 15 2017
 */

// Libraries
#include <ctime>
#include "clock.h"
#include<time.h>
#include<sys/time.h>

namespace rockseis {
// =========== CLOCK CLASS =========== //

template<typename T>
Clock<T>::Clock() {

tend=((long double)std::clock())/CLOCKS_PER_SEC;
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
struct timeval time;

tstart=tend;
gettimeofday(&time,NULL);
wtstart=(T)time.tv_sec + (T)time.tv_usec*0.000001;

}

template<typename T>
void Clock<T>::stop() {
struct timeval time;

tend=((long double)std::clock())/CLOCKS_PER_SEC;
telapsed += tend-tstart;
tstart=tend;
gettimeofday(&time,NULL);
wtend=(T)time.tv_sec + (T)time.tv_usec*0.000001;
wtelapsed=wtend-wtstart+wtelapsed;
}

template<typename T>
T Clock<T>::elapsed() {

tend=((long double)std::clock())/CLOCKS_PER_SEC;
return(telapsed+tend-tstart);
}

template<typename T>
void Clock<T>::print(){

std::cout <<" Cpu  time (seconds) : " << telapsed << "\n";
std::cout <<" Wall time (seconds) : " << wtelapsed << "\n";
}


// =========== TEMPLATE CLASS INSTANTIATION =========== //
template class Clock<double>;
template class Clock<float>;

}
