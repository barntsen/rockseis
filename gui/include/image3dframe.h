#ifndef IMAGE3DFRAME_H
#define IMAGE3DFRAME_H

//(*Headers(Image3dframe)
#include <wx/wxprec.h>
#include <wx/app.h>
#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif

#include <wx/textdlg.h>
#include <wx/toolbar.h>
#include <wx/sizer.h>
#include <wx/menu.h>
#include <wx/panel.h>
#include <wx/scrolwin.h>
#include <wx/statusbr.h>
#include <wx/frame.h>
#include "wx/bitmap.h"
#include "wx/dcbuffer.h"
#include "wx/msgdlg.h"
#include "wx/filedlg.h"
#include <wx/numdlg.h>
#include <wx/string.h>
#include <wx/intl.h>
#include <wx/image.h>
#include <wx/artprov.h>
#include "zoom.h"
#include <algorithm>
#include <vector>
#include "colors.h"

#define FONTSIZE1 6 
#define FONTSIZE2 9 
#define NCOLORS 5
#define MAXPOS 9999

class Image3dframe: public wxFrame
{
	public:
        Image3dframe(size_t _n1, float _d1, float _o1, size_t _n2, float _d2, float _o2,  size_t _n3, float _d3, float _o3, float *_imagedata, wxWindow* parent, wxWindowID id=wxID_ANY,const wxPoint& pos=wxDefaultPosition,const wxSize& size=wxDefaultSize);
		virtual ~Image3dframe();

		wxStatusBar* StatusBar1;
		wxPanel* Xaxis;
		wxPanel* Yaxis_vert;
		wxPanel* Yaxis_hor;
		wxPanel* Zaxis;
		wxPanel* Cornerup;
		wxPanel* Cornerdown;
		wxPanel* Slider1panel;
		wxPanel* Slider2panel;
		wxPanel* Slider3panel;
		wxPanel* Gap1;
		wxPanel* Gap2;
		wxPanel* Gap3;
		wxPanel* Gap4;
		wxPanel* Xplotwindow;
		wxPanel* Zplotwindow;
		wxPanel* Yplotwindow;
		wxSlider* Slider1;
		wxSlider* Slider2;
		wxSlider* Slider3;

        void setToolbarset(bool val) { toolbarset = val; } 
        bool getToolbarset() { return toolbarset; } 
        void createToolbar();
        void destroyToolbar();

        void setMenubarset(bool val) { menubarset = val; } 
        bool getMenubarset() { return menubarset; } 
        void createMenubar();
        void LoadXimage(int zi1, int zn1, int zi2, int zn2);
        void LoadYimage(int zi1, int zn1, int zi2, int zn2);
        void LoadZimage(int zi1, int zn1, int zi2, int zn2);
        bool getXimageAlloc() { return ximage_allocated; }
        void setXimageAlloc(bool val) { ximage_allocated = val; }
        bool getYimageAlloc() { return yimage_allocated; }
        void setYimageAlloc(bool val) { yimage_allocated = val; }
        bool getZimageAlloc() { return zimage_allocated; }
        void setZimageAlloc(bool val) { zimage_allocated = val; }
        void setCmpnumber(int val) { cmpnumber = val; }
        int getCmpnumber() { return cmpnumber; }
        void setDcmp(int val) { dcmp = val; }
        int getDcmp() { return dcmp; }
        void setMaxcmp(int val) { maxcmp = val; }
        int getMaxcmp() { return maxcmp; }
        void setStatusbar(int cmp, float x, float y, float val);
        void ComputeClip();
        void setDisplaycrosshair(bool val) { displaycrosshair = val; }
        void setGetcrosshair(bool val) { getcrosshair = val; }
        void setN1(size_t val) { n1 = val; }
        void setD1(float val) { d1 = val; }
        void setO1(float val) { o1 = val; }

        void setN2(size_t val) { n2 = val; }
        void setD2(float val) { d2 = val; }
        void setO2(float val) { o2 = val; }

        void setN3(size_t val) { n3 = val; }
        void setD3(float val) { d3 = val; }
        void setO3(float val) { o3 = val; }

        size_t getN1() { return n1; }
        float getD1() { return d1; }
        float getO1() { return o1; }

        size_t getN2() { return n2; }
        float getD2() { return d2; }
        float getO2() { return o2; }

        size_t getN3() { return n3; }
        float getD3() { return d3; }
        float getO3() { return o3; }
        void setImagedata (float *data) { imagedata = data; }
        void setXzoom(float o1, float max1, float o2, float max2, int ix0, int nx, int iy0, int ny);
        void setYzoom(float o1, float max1, float o2, float max2, int ix0, int nx, int iy0, int ny);
        void setZzoom(float o1, float max1, float o2, float max2, int ix0, int nx, int iy0, int ny);

        float *getCrosshair_pt() { return &crosshair_pt[0]; }

        void createPicks();
        void createPicks(int dir, long n, float d, float o);

        void PlotXcrosshair(wxDC &dc, int w, int h);
        void PlotYcrosshair(wxDC &dc, int w, int h);
        void PlotZcrosshair(wxDC &dc, int w, int h);
    protected:

		//(*Identifiers(StackVel3dframe)
		static const long ID_PANEL1;
		static const long ID_PANEL2;
		static const long ID_PANEL3;
		static const long ID_PANEL4;
		static const long ID_PANEL5;
		static const long ID_PANEL6;
		static const long ID_PANEL7;
		static const long ID_PANEL8;
		static const long ID_PANEL9;
		static const long ID_PANEL10;
		static const long ID_PANEL11;
		static const long ID_PANEL12;
		static const long ID_PANEL13;
		static const long ID_PANEL14;
		static const long ID_PANEL15;
		static const long ID_PANEL16;
		static const long ID_STATUSBAR1;
		static const long ID_SLIDER1;
		static const long ID_SLIDER2;
		static const long ID_SLIDER3;
		//*)


	private:
		//(*Handlers(Image2dframe)
		void OnZaxisPaint(wxPaintEvent& event);
		void OnXaxisPaint(wxPaintEvent& event);
		void OnYaxis_vertPaint(wxPaintEvent& event);
		void OnYaxis_horPaint(wxPaintEvent& event);
		void OnXplotwindowPaint(wxPaintEvent& event);
		void OnYplotwindowPaint(wxPaintEvent& event);
		void OnZplotwindowPaint(wxPaintEvent& event);
		void OnSlider1CmdScroll(wxScrollEvent& event);
		void OnSlider2CmdScroll(wxScrollEvent& event);
		void OnSlider3CmdScroll(wxScrollEvent& event);
		void OnKeyUp(wxKeyEvent& event);
        void OnXplotwindowMouseMove(wxMouseEvent& event);
        void OnYplotwindowMouseMove(wxMouseEvent& event);
        void OnZplotwindowMouseMove(wxMouseEvent& event);
		void OnXplotwindowLeftUp(wxMouseEvent& event);
		void OnXplotwindowLeftDown(wxMouseEvent& event);
		void OnYplotwindowLeftUp(wxMouseEvent& event);
		void OnYplotwindowLeftDown(wxMouseEvent& event);
		void OnZplotwindowLeftUp(wxMouseEvent& event);
		void OnZplotwindowLeftDown(wxMouseEvent& event);

        //Variables
        size_t n1;
        size_t n2;
        size_t n3;
        float d1;
        float d2;
        float d3;
        float o1;
        float o2;
        float o3;
        float *imagedata;
        float pclip[100];
        unsigned char RGB[768]; 
        int iminclip, imaxclip;
        wxImage *ximage;
        wxImage *yimage;
        wxImage *zimage;
        bool ximage_allocated;
        bool yimage_allocated;
        bool zimage_allocated;
        Zoom *xslicezoom;
        Zoom *yslicezoom;
        Zoom *zslicezoom;
        int color; 
        bool toolbarset;
        bool menubarset;
        int cmpnumber, dcmp, maxcmp;
        float crosshair_pt[2];
        bool displaycrosshair;
        bool getcrosshair;

		DECLARE_EVENT_TABLE()
};


#endif
