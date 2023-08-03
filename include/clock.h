/* Author: B. Arntsen <borge.arntsen@ntnu.no>
 * Date: November 20 2017
 * Small library for measuring cpu time
 *
 */

#ifndef CLOCK_H
#define CLOCK_H
#include<iostream>

namespace rockseis{
// =========== CLOCK CLASS =========== //

template<typename T>
class Clock {
public:
    Clock();                               ///< Default constructor
    ~Clock();                              ///< Destructor

    // Functions
    T  wtime();                         ///< Get Wall clock time

    // Functions
    void start();                      ///< Start clock
    void stop();                       ///< End clock
    void print();                      ///< Print results

private:
    // Variables
    T tstart, tend, telapsed;  // Cpu time
    T wtstart,wtend,wtelapsed; // Wall clock time
};

}
#endif //CLOCK_H
