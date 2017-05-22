#ifndef REVOLVE_H
#define REVOLVE_H

// Include statements
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <time.h>
#include "file.h"
#include "waves.h"

#define REVOLVE_OK 1
#define REVOLVE_ERR 0

/* Constants */
#define CHECKUP 512 
#define REPSUP 64 
#define MAXINT 2147483647


namespace rockseis {

typedef enum { advance, takeshot, restore, firsturn, youturn, terminate, error} revolve_action;
/*^*/

typedef struct 
 {
  int advances;
  int takeshots;
  int commands;
 } revolve_nums;
/*^*/

// =============== ABSTRACT REVOLVE CLASS =============== //
/** The abstract chkpointing class
 *
 */

template<typename T>
class Revolve {
public:
    Revolve();	///< Constructor
    Revolve(unsigned int _nt, unsigned int _nsnaps);	///< Constructor
    int adjust(int nsteps); ///< Automatically compute optimal number of snaps
    double expense(int steps, int snaps); ///< Get ratio between full and optimal checkpointing forward modelling steps
    int numforw(int steps, int snaps); ///< Get total number of forward modellings
    revolve_action revolve(); ///< Get what to do
    int getCapo() { return capo; } ///< Get capo
    int getInfo() { return info; }
    int getFine() { return fine; }
    int getCheck() { return check; }
    void readCheckpoint();
    void writeCheckpoint();
    
    // Revolve functions
    ~Revolve();	///< Destructor

private:
    int maxrange(int ss, int tt);
	int check;
    int capo;
    int fine;
    int steps;
    int snaps;
    int nfw;
	int info;
    bool incore;
    T *chp;
    std::shared_ptr<File> Fc;
    std::string filename;
    revolve_nums numbers;
    int  turn;
    int reps; 
    int range;
    int ch[CHECKUP];
    int oldsnaps;
    int oldfine;
};

}
#endif //REVOLVE_H
