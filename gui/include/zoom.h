/***************************************************************
 * Name:      zoom.h
 * Purpose:   Defines Zoom Class
 * Author:    Wiktor Weibull (wiktor@rockseis.no)
 * Created:   2016-01-25
 * Copyright: Wiktor Weibull ()
 * License:
 **************************************************************/

#ifndef ZOOM_H
#define ZOOM_H
#include <wx/gdicmn.h>

class Zoom
{
	private:
		wxPoint zpos1,zpos2; // topleft and bottom right corners
		wxPoint box[5];  // Zoom box to draw
		float x0, x1;   // Smallest and largest xs
		float y0, y1;   // Smallest and largest ys
		float z0, z1;   // Smallest and largest zs
		float xfactor, yfactor, zfactor; // Zoom factors
		int iy0;
		int ny;
		int ix0;
		int nx;
		int iz0;
		int nz;
		bool zooming;


	public:
		Zoom() { /* Default constructor */ }
		Zoom(float ya, float yb, float xa, float xb, int i, int n, int j, int m);
		Zoom(float ya, float yb, float xa, float xb, float za, float zb, int i, int n, int j, int m, int k, int l);
		~Zoom();
		void Setx0(float v) { x0 = v; }
		void Setx1(float v) { x1 = v; }
		void Sety0(float t) { y0 = t; }
		void Sety1(float t) { y1 = t; }
		void Setz0(float v) { z0 = v; }
		void Setz1(float v) { z1 = v; }
		float Getx0() { return x0; }
		float Getx1() { return x1; }
		float Gety0() { return y0; }
		float Gety1() { return y1; }
		float Getz0() { return z0; }
		float Getz1() { return z1; }
		wxPoint *Getbox() { return &box[0]; }
		void Setbox(wxPoint pos, int i) { box[i].x = pos.x; box[i].y = pos.y;}
		wxPoint Getzpos1() { return zpos1; }
		void Setzpos1(wxPoint pos) { zpos1=pos; }
		wxPoint Getzpos2() { return zpos2; }
		void Setzpos2(wxPoint pos) { zpos2=pos; }
		void Setxfactor(float factor) { xfactor = factor; }
		float Getxfactor() { return xfactor; }
		void Setyfactor(float factor) { yfactor = factor; }
		float Getyfactor() { return yfactor; }
		void Setzfactor(float factor) { zfactor = factor; }
		float Getzfactor() { return zfactor; }
		int Getix0() { return ix0; }
		void Setix0(int ix) { ix0 = ix; }
		int Getnx() { return nx; }
		void Setnx(int n) {  nx = n; }
		int Getny() { return ny; }
		void Setny(int n) { ny = n; }
		int Getiy0() { return iy0; }
		void Setiy0(int iy) { iy0 = iy; }
		int Getnz() { return nz; }
		void Setnz(int n) { nz = n; }
		int Getiz0() { return iz0; }
		void Setiz0(int iz) { iz0 = iz; }
		void Setzooming(bool zoom) { zooming=zoom; }
		bool Getzooming() { return zooming; }
};
#endif // ZOOM_H
