/***************************************************************
 * Name:      Picks.h
 * Purpose:   Defines Picks Class
 * Author:    Wiktor Weibull (wiktor@rockseis.no)
 * Created:   2016-01-25
 * Copyright: Wiktor Weibull ()
 * License:
 **************************************************************/

#ifndef PICKS_H
#define PICKS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <wx/filefn.h>
#include "utils.h"

#define PICK_HORIZONTAL 0
#define PICK_VERTICAL 1

#define PICKS_OK 0
#define PICKS_ERR 1

class Picks
{
	private:
	    int nt;
	    float dt;
	    float ot;
		int ncmp;
		int maxpicks;
        rockseis::Point2D<float> *picks;
		int *npicks;
		float *vrms;
		bool changed;
        std::string filename;
        int type;


	public:
	    Picks() { /* Do nothing */ }
	    Picks(int maxp, int nc, int n, float d, float o, int type);
	    rockseis::Point2D<float> *Getpicks() { return picks; }
	    int *Getnpicks() { return npicks; }
	    int Getmaxpicks() { return maxpicks; }
	    bool Checkforpicks();
	    float *Getvrms() { return vrms; }
	    void Addpick(int n, float x, float y);
	    void Removepick(int n, float x, float y);
	    void Interp(int cmp_number);
	    void Project();
	    int Savepicks(std::string filename);
	    int Loadpicks(std::string filename);
	    void Importpicks(char *filename);
	    bool Getchanged() { return changed; }
	    void Setchanged(bool update) { changed = update; }
	    void Clearpicks();
	    ~Picks();
};


#endif // PICKS_H

