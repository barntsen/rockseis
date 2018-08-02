/***************************************************************
 * Name:      zoom.cpp
 * Purpose:   Code for Zoom Class functions
 * Author:    Wiktor Weibull (wiktor@rockseis.no)
 * Created:   2016-01-25
 * Copyright: Wiktor Weibull ()
 * License:
 **************************************************************/

#include "zoom.h"


// Member function implementations
Zoom::Zoom(float xa, float xb, float ya, float yb, int i, int n, int j, int m)
{
    int a;
    for(a=0; a<5; a++){
        box[a].x=0;
        box[a].y=0;
    }
    zpos1.x=0;
    zpos1.y=0;
    zpos2.x=0;
    zpos2.y=0;
    y0=ya; y1=yb;
    x0=xa; x1=xb;
    xfactor=1.0;
    yfactor=1.0;
    ix0=i;
    nx=n;
    iy0=j;
    ny=m;
    zooming=0;

}

Zoom::Zoom(float xa, float xb, float ya, float yb,float za, float zb, int i, int n, int j, int m, int k, int l)
/* Overloaded function for 3D plotting version */
{
    int a;
    for(a=0; a<5; a++){
        box[a].x=0;
        box[a].y=0;
    }
    zpos1.x=0;
    zpos1.y=0;
    zpos2.x=0;
    zpos2.y=0;
    y0=ya; y1=yb;
    x0=xa; x1=xb;
    z0=za; z1=zb;
    xfactor=1.0;
    yfactor=1.0;
    ix0=i;
    nx=n;
    iy0=j;
    ny=m;
    iz0=k;
    nz=l;
    zooming=0;

}


Zoom::~Zoom()
{
    // Nothing to do
}


