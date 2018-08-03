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

class Picks
{
	private:
	    float dt;
	    int nt;
		int ncmp;
		int maxpicks;
        rockseis::Point2D<float> *picks;
		int *npicks;
		float *vrms;
		bool changed;


	public:
	    Picks() { /* Do nothing */ }
	    Picks(int maxp, int nc, int n, float d);
	    rockseis::Point2D<float> *Getpicks() { return picks; }
	    int *Getnpicks() { return npicks; }
	    int Getmaxpicks() { return maxpicks; }
	    bool Checkforpicks();
	    float *Getvrms() { return vrms; }
	    void Addpick(int n, float x, float y);
	    void Removepick(int n, float x, float y);
	    void Interp(int cmp_number);
	    void Savepicks(char *folder);
	    void Loadpicks(char *folder);
	    void Importpicks(char *filename);
	    bool Getchanged() { return changed; }
	    void Setchanged(bool update) { changed = update; }
	    void Clearpicks();
	    ~Picks();
};


#endif // PICKS_H

