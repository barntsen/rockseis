#ifndef IMAGE2DFRAME_H
#define IMAGE2DFRAME_H

//(*Headers(Image2dframe)
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
#include "picks.h"
#include <algorithm>
#include <vector>

#define FONTSIZE1 6 
#define FONTSIZE2 9 
#define NCOLORS 5
#define MAXPOS 9999


class Image2dframe: public wxFrame
{
	public:
		Image2dframe(size_t n1, float d1, float o1, size_t n2, float d2, float o2, float *imagedata, wxWindow* parent,wxWindowID id=wxID_ANY,const wxPoint& pos=wxDefaultPosition,const wxSize& size=wxDefaultSize);
		virtual ~Image2dframe();
        void setToolbarset(bool val) { toolbarset = val; } 
        bool getToolbarset() { return toolbarset; } 
        void createToolbar();
        void destroyToolbar();
        void LoadImage(int zi1, int zn1, int zi2, int zn2);
        void setCmpnumber(int val) { cmpnumber = val; }
        int getCmpnumber() { return cmpnumber; }
        void setDcmp(int val) { dcmp = val; }
        int getDcmp() { return dcmp; }
        void setMaxcmp(int val) { maxcmp = val; }
        int getMaxcmp() { return maxcmp; }
        void setStatusbar(int cmp, float x, float y, float val);
        void ComputeClip();
        void setDisplaypicks(int val) { displaypicks = val; }
        void setN1(size_t val) { n1 = val; }
        void setD1(float val) { d1 = val; }
        void setO1(float val) { o1 = val; }

        void setN2(size_t val) { n2 = val; }
        void setD2(float val) { d2 = val; }
        void setO2(float val) { o2 = val; }

        size_t getN1() { return n1; }
        float getD1() { return d1; }
        float getO1() { return o1; }

        size_t getN2() { return n2; }
        float getD2() { return d2; }
        float getO2() { return o2; }
        void setImagedata (float *data) { imagedata = data; }
        void setZoom(float o1, float max1, float o2, float max2, int ix0, int nx, int iy0, int ny);

        void createPicks();

		//(*Declarations(Image2dframe)
		wxPanel* Imagewindow;
		wxStatusBar* StatusBar1;
		wxPanel* LeftCorner;
		wxPanel* Zaxis;
		wxPanel* Xaxis;
		wxPanel* RightCorner;

		wxToolBar* ToolBar1;
		wxToolBarToolBase* ToolBarItem1;
		wxToolBarToolBase* ToolBarItem2;
		wxToolBarToolBase* ToolBarItem3;
		wxToolBarToolBase* ToolBarItem4;
		wxToolBarToolBase* ToolBarItem5;
		wxToolBarToolBase* ToolBarItem6;
		wxToolBarToolBase* ToolBarItem7;
		//*)

	protected:

		//(*Identifiers(Image2dframe)
		static const long ID_PANEL1;
		static const long ID_PANEL2;
		static const long ID_PANEL3;
		static const long ID_PANEL4;
		static const long ID_IMAGEWINDOW1;
		static const long ID_STATUSBAR1;
		static const long idToolNext;
		static const long idToolcmpint;
		static const long idToolPrev;
		static const long idToolpick;
		static const long idToolzoom;
		static const long idToolSave;
		static const long idToolLoad;
		static const long ID_TOOLBAR1;
		static const long ID_LISTBOX1;
		//*)

	private:
		//(*Handlers(Image2dframe)
		void OnImagewindowPaint(wxPaintEvent& event);
		void OnZaxisPaint(wxPaintEvent& event);
		void OnXaxisPaint(wxPaintEvent& event);
		void OnImagewindowMouseMove(wxMouseEvent& event);
		void OnImagewindowLeftUp(wxMouseEvent& event);
		void OnImagewindowRightUp(wxMouseEvent& event);
		void OnImagewindowLeftDown(wxMouseEvent& event);
		void OnImagewindowKeyUp(wxKeyEvent& event);
		void OnNextClicked(wxCommandEvent& event);
		void OnPreviousClicked(wxCommandEvent& event);
		void OnCMPinterval(wxCommandEvent& event);
		void OnClose(wxCloseEvent& event);
		void OnImagewindowEraseBackground(wxEraseEvent& event);
        bool getImage2dAlloc() { return image2d_allocated; }
        void setImage2dAlloc(bool val) { image2d_allocated = val; }
        void getRgb(int color);
        void Plotpicks(wxDC &dc, int w, int h);
        void OnPick(wxCommandEvent& event);
        void OnZoom(wxCommandEvent& event);

        //Variables
        size_t n1;
        size_t n2;
        float d1;
        float d2;
        float o1;
        float o2;
        float *imagedata;
        bool image2d_allocated;
        float pclip[100];
        unsigned char RGB[768]; 
        int iminclip, imaxclip;
        wxImage *image2d;
        Zoom *zoom;
        int color; 
        bool toolbarset;
        int cmpnumber, dcmp, maxcmp;
        bool displaypicks;

        // Picking variables
        std::vector<Picks*> picks;
        int nlayers;
        int layer;
        wxPoint pos[MAXPOS];

		//*)

		DECLARE_EVENT_TABLE()
};

wxDECLARE_EVENT(SelectCmp, wxCommandEvent);

#endif
