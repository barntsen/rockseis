/***************************************************************
 * Name:      Picks.cpp
 * Purpose:   Code for Picks Class functions
 * Author:    Wiktor Weibull (wiktor@rockseis.no)
 * Created:   2016-01-25
 * Copyright: Wiktor Weibull ()
 * License:
 **************************************************************/

#include "picks.h"
#define SQ(x) ((x) * (x))

// Member function implementations
Picks::Picks(int maxp, int nc, int n, float d)
{
    dt=d;
    nt=n;
    ncmp=nc;
    maxpicks=maxp;

    //Allocating picks array
    picks = (rockseis::Point2D<float> *) calloc(ncmp*maxpicks, sizeof(rockseis::Point2D<float>));

    // Allocating npicks array
    npicks = (int *) calloc(ncmp, sizeof(int));

    //Allocating VRMS
    vrms = (float *) calloc(nt, sizeof(float));
    changed=0;
}

void Picks::Addpick(int n, float x, float y)
{
    rockseis::Point2D<float> evpos;
    evpos.x=x;
    evpos.y=y;
    float dy, dist, mindist;
    int i, minindex;

    if(npicks[n]==0){
        picks[n*maxpicks].x=evpos.x;
        picks[n*maxpicks].y=evpos.y;
        npicks[n]++;
    }else{
        for(i=0; i<npicks[n]; i++){
            dy=picks[n*maxpicks + i].y - evpos.y;
            dist=sqrtf(dy*dy);
            if(i==0){
              mindist=dist;
              minindex=0;
            }
            if(dist<mindist){
                mindist=dist;
                minindex=i;
            }
        }
        if(evpos.y < picks[n*maxpicks + minindex].y){
                npicks[n]++;
            for(i=npicks[n]-1; i>minindex; i--){
                picks[n*maxpicks + i].x=picks[n*maxpicks + i-1].x;
                picks[n*maxpicks + i].y=picks[n*maxpicks + i-1].y;
            }
            picks[n*maxpicks + minindex].x=evpos.x;
            picks[n*maxpicks + minindex].y=evpos.y;
        }
        if(evpos.y > picks[n*maxpicks + minindex].y){
                npicks[n]++;
            for(i=npicks[n]-1; i>minindex+1; i--){
                picks[n*maxpicks + i].x=picks[n*maxpicks + i-1].x;
                picks[n*maxpicks + i].y=picks[n*maxpicks + i-1].y;
            }
            picks[n*maxpicks + minindex+1].x=evpos.x;
            picks[n*maxpicks + minindex+1].y=evpos.y;
        }
    }
    changed=1;
}

void Picks::Removepick(int n, float x, float y)
{
    rockseis::Point2D<float> evpos;
    evpos.x=x;
    evpos.y=y;
    float dy, dist, mindist;
    int i, minindex;
    minindex=0;
    dist=0;
    mindist=0;

    if(npicks[n]>0){
        for(i=0; i<npicks[n]; i++){
            dy=picks[n*maxpicks + i].y - evpos.y;
            dist=sqrt(dy*dy);
            if(i==0) mindist=dist;
            if(dist<mindist){
                mindist=dist;
                minindex=i;
            }
        }
        npicks[n]--;
        for(i=minindex; i<npicks[n]; i++){
            picks[n*maxpicks + i].x=picks[n*maxpicks + i+1].x;
            picks[n*maxpicks + i].y=picks[n*maxpicks + i+1].y;
        }
    }
    changed=1;
}

int seekindex(const float *buffer, const int n, const float value)
/*< Find interpolation position through binary search >*/
{
	int Mid;
	int High;
	int Low;

	High=n-1;
	Low=0;

        if(value<=buffer[Low]){
		return (Low);
	}
        if(value>=buffer[High]){
		return (High);
	}

	while(1){
		Mid=(High+Low)/2;
		if(value==buffer[Mid]){
			return (Mid);
		}
		if(value>buffer[Mid]){
			Low=Mid;
		}
		else{
			High=Mid;
		}
		if((High-Low)<2){
			return (Low);
		}
	}
}

void Picks::Interp(int cmp_number)
{
    int i,it;
    float ti,d,p0,p1;
    int n = npicks[cmp_number];
    if(n>0){
        float *t;
        t = (float *) calloc((n+1), sizeof(float));

        for(i=0; i < n; i++){
            t[i]=picks[maxpicks*cmp_number + i].y;
        }
        for(i=0; i<nt; i++)
        {
            ti=i*dt;
            if(ti <= t[0] || ti > t[n-1])
            {
                if(ti <= t[0]) vrms[i]=picks[maxpicks*cmp_number].x;
                if(ti > t[n-1]) vrms[i]=picks[maxpicks*cmp_number+n-1].x;
            }else{
                it=seekindex(t, n, ti);
                d=(ti-t[it])/(t[it+1]-t[it]);
                p0=picks[maxpicks*cmp_number + it].x;
                p1=picks[maxpicks*cmp_number + it+1].x;
                vrms[i]=p0*(1.0-d) + p1*d;
            }
        }
    }else{
        for(i=0; i<nt; i++){
            vrms[i]=0.0;
        }
    }
}

void Picks::Savepicks(char *folder)
{
    // No need to save
    changed=0;
}

void Picks::Loadpicks(char *folder)
{
    FILE *fp;

    char filename[256];
    int i,j;

    // Get Semblance data
    snprintf(filename, 256, "%s/Picks/picksdata.bin", folder);
    if(wxFileExists(filename)){
        fp=fopen(filename, "rb");
        int nump;
        float x,y;
        for (i=0; i < ncmp; i++){
            fread(&nump, sizeof(int), 1, fp);
            npicks[i]=nump;
            for(j=0; j<nump; j++){
                fread(&x, sizeof(float), 1, fp);
                fread(&y, sizeof(float), 1, fp);
                picks[i*maxpicks +j].x = x;
                picks[i*maxpicks +j].y = y;
            }
        }
        fclose(fp);
    }
}

void Picks::Importpicks(char *filename)
{
    FILE *fp;

    int i,j;

    // Get Picks data
    if(wxFileExists(filename)){
        fp=fopen(filename, "rb");
        int nump;
        float x,y;
        for (i=0; i < ncmp; i++){
            fread(&nump, sizeof(int), 1, fp);
            npicks[i]=nump;
            for(j=0; j<nump; j++){
                fread(&x, sizeof(float), 1, fp);
                fread(&y, sizeof(float), 1, fp);
                picks[i*maxpicks +j].x = x;
                picks[i*maxpicks +j].y = y;
            }
        }
        fclose(fp);
    }
}

void Picks::Clearpicks()
{
    int i;
    for (i=0; i<ncmp; i++){
        npicks[i]=0;
    }
}

bool Picks::Checkforpicks()
{
    int i;
    for (i=0; i<ncmp; i++){
        if(npicks[i]) return true;
    }
    return false;
}

Picks::~Picks()
{
    //dtor
    free(picks);
    free(npicks);
    free(vrms);
}
