#include "image3dframe.h"

//(*IdInit(Image3dframe)
const long Image3dframe::ID_PANEL1 = wxNewId();
const long Image3dframe::ID_PANEL2 = wxNewId();
const long Image3dframe::ID_PANEL3 = wxNewId();
const long Image3dframe::ID_PANEL4 = wxNewId();
const long Image3dframe::ID_PANEL5 = wxNewId();
const long Image3dframe::ID_PANEL6 = wxNewId();
const long Image3dframe::ID_PANEL7 = wxNewId();
const long Image3dframe::ID_PANEL8 = wxNewId();
const long Image3dframe::ID_PANEL9 = wxNewId();
const long Image3dframe::ID_PANEL10 = wxNewId();
const long Image3dframe::ID_PANEL11 = wxNewId();
const long Image3dframe::ID_PANEL12 = wxNewId();
const long Image3dframe::ID_PANEL13 = wxNewId();
const long Image3dframe::ID_PANEL14 = wxNewId();
const long Image3dframe::ID_PANEL15 = wxNewId();
const long Image3dframe::ID_PANEL16 = wxNewId();
const long Image3dframe::ID_STATUSBAR1 = wxNewId();
const long Image3dframe::ID_SLIDER1 = wxNewId();
const long Image3dframe::ID_SLIDER2 = wxNewId();
const long Image3dframe::ID_SLIDER3 = wxNewId();
//*)

BEGIN_EVENT_TABLE(Image3dframe,wxFrame)
	//(*EventTable(Image3dframe)
	//*)
END_EVENT_TABLE()

Image3dframe::Image3dframe(size_t _n1, float _d1, float _o1, size_t _n2, float _d2, float _o2,  size_t _n3, float _d3, float _o3, float *_imagedata, wxWindow* parent, wxWindowID id,const wxPoint& pos,const wxSize& size)
{
    //Setting up variables
    n1 = _n1;
    d1 = _d1;
    o1 = _o1;

    n2 = _n2;
    d2 = _d2;
    o2 = _o2;

    n3 = _n3;
    d3 = _d3;
    o3 = _o3;

    imagedata = _imagedata;
    ximage_allocated = false;
    yimage_allocated = false;
    zimage_allocated = false;
    toolbarset = false;
    menubarset = false;

    cmpnumber = 0;
    dcmp = 1;
    maxcmp = 1;

	//(*Initialize(Image3dframe)
	wxFlexGridSizer* FlexGridSizer1;

	Create(parent, wxID_ANY, _("3D Stack/Velocity viewer"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE, _T("wxID_ANY"));
	FlexGridSizer1 = new wxFlexGridSizer(4, 4, 0, 0);
	FlexGridSizer1->AddGrowableCol(1);
	FlexGridSizer1->AddGrowableRow(0);
	FlexGridSizer1->AddGrowableCol(3);
	FlexGridSizer1->AddGrowableRow(2);
	Yaxis_vert = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxSize(50,300), wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1->Add(Yaxis_vert, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Zplotwindow = new wxPanel(this, ID_PANEL2, wxDefaultPosition, wxSize(300,300), wxTAB_TRAVERSAL, _T("ID_PANEL2"));
	FlexGridSizer1->Add(Zplotwindow, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Slider1panel = new wxPanel(this, ID_PANEL3, wxDefaultPosition, wxSize(50,300), wxTAB_TRAVERSAL, _T("ID_PANEL3"));
	Slider1 = new wxSlider(Slider1panel, ID_SLIDER1, 0, 0, 10, wxPoint(1,50), wxSize(50,200), wxSL_VERTICAL|wxSL_LABELS, wxDefaultValidator, _T("ID_SLIDER1"));
	FlexGridSizer1->Add(Slider1panel, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Cornerup = new wxPanel(this, ID_PANEL4, wxDefaultPosition, wxSize(300,300), wxTAB_TRAVERSAL, _T("ID_PANEL4"));
	FlexGridSizer1->Add(Cornerup, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Gap1 = new wxPanel(this, ID_PANEL5, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL5"));
	FlexGridSizer1->Add(Gap1, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
    
	Slider2panel = new wxPanel(this, ID_PANEL6, wxDefaultPosition, wxSize(300,50), wxTAB_TRAVERSAL, _T("ID_PANEL6"));
	Slider2 = new wxSlider(Slider2panel, ID_SLIDER2, 0, 0, 10, wxPoint(50,1), wxSize(200,50), wxSL_HORIZONTAL|wxSL_LABELS, wxDefaultValidator, _T("ID_SLIDER2"));
	FlexGridSizer1->Add(Slider2panel, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Gap2 = new wxPanel(this, ID_PANEL7, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL7"));
	FlexGridSizer1->Add(Gap2, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Slider3panel = new wxPanel(this, ID_PANEL8, wxDefaultPosition, wxSize(300,50), wxTAB_TRAVERSAL, _T("ID_PANEL8"));
	Slider3 = new wxSlider(Slider3panel, ID_SLIDER2, 0, 0, 10, wxPoint(50,1), wxSize(200,50), wxSL_HORIZONTAL|wxSL_LABELS, wxDefaultValidator, _T("ID_SLIDER3"));
	FlexGridSizer1->Add(Slider3panel, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Zaxis = new wxPanel(this, ID_PANEL9, wxDefaultPosition, wxSize(50,300), wxTAB_TRAVERSAL, _T("ID_PANEL9"));
	FlexGridSizer1->Add(Zaxis, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Xplotwindow = new wxPanel(this, ID_PANEL10, wxDefaultPosition, wxSize(300,300), wxTAB_TRAVERSAL, _T("ID_PANEL10"));
	FlexGridSizer1->Add(Xplotwindow, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Gap3 = new wxPanel(this, ID_PANEL11, wxDefaultPosition, wxSize(50,300), wxTAB_TRAVERSAL, _T("ID_PANEL11"));
	FlexGridSizer1->Add(Gap3, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Yplotwindow = new wxPanel(this, ID_PANEL12, wxDefaultPosition, wxSize(300,300), wxTAB_TRAVERSAL, _T("ID_PANEL12"));
	FlexGridSizer1->Add(Yplotwindow, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Cornerdown = new wxPanel(this, ID_PANEL13, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL13"));
	FlexGridSizer1->Add(Cornerdown, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Xaxis = new wxPanel(this, ID_PANEL14, wxDefaultPosition, wxSize(300,50), wxTAB_TRAVERSAL, _T("ID_PANEL14"));
	FlexGridSizer1->Add(Xaxis, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Gap4 = new wxPanel(this, ID_PANEL15, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL15"));
	FlexGridSizer1->Add(Gap4, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);

	Yaxis_hor = new wxPanel(this, ID_PANEL16, wxDefaultPosition, wxSize(300,50), wxTAB_TRAVERSAL, _T("ID_PANEL16"));
	FlexGridSizer1->Add(Yaxis_hor, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	SetSizer(FlexGridSizer1);
	StatusBar1 = new wxStatusBar(this, ID_STATUSBAR1, 0, _T("ID_STATUSBAR1"));
	int __wxStatusBarWidths_1[1] = { -10 };
	int __wxStatusBarStyles_1[1] = { wxSB_NORMAL };
	StatusBar1->SetFieldsCount(1,__wxStatusBarWidths_1);
	StatusBar1->SetStatusStyles(1,__wxStatusBarStyles_1);
	SetStatusBar(StatusBar1);

	Xaxis->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnXaxisPaint,0,this);
	Yaxis_vert->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnYaxis_vertPaint,0,this);
	Yaxis_hor->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnYaxis_horPaint,0,this);
	Zaxis->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnZaxisPaint,0,this);
	Xplotwindow->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnXplotwindowPaint,0,this);
	Yplotwindow->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnYplotwindowPaint,0,this);
	Zplotwindow->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image3dframe::OnZplotwindowPaint,0,this);
    Slider1->Connect(wxEVT_SLIDER, wxScrollEventHandler(Image3dframe::OnSlider1CmdScroll),NULL, this);
    Slider2->Connect(wxEVT_SLIDER, wxScrollEventHandler(Image3dframe::OnSlider2CmdScroll),NULL, this);
    Slider3->Connect(wxEVT_SLIDER, wxScrollEventHandler(Image3dframe::OnSlider3CmdScroll),NULL, this);
	Xplotwindow->Connect(wxEVT_KEY_UP,(wxObjectEventFunction)&Image3dframe::OnKeyUp,0,this);
	Yplotwindow->Connect(wxEVT_KEY_UP,(wxObjectEventFunction)&Image3dframe::OnKeyUp,0,this);
	Zplotwindow->Connect(wxEVT_KEY_UP,(wxObjectEventFunction)&Image3dframe::OnKeyUp,0,this);

    // Compute Maxmimum and minimum values
    this->ComputeClip();

    // Setup color pallete
    color = 1;
    getRgb(color, &RGB[0]);

    // Set zoom 
    xslicezoom = new Zoom(o1, (n1-1)*d1 + o1, o3, (n3-1)*d3 + o3, 0, n1, 0, n3);
    yslicezoom = new Zoom(o2, (n2-1)*d2 + o2, o3, (n3-1)*d3 + o3, 0, n2, 0, n3);
    zslicezoom = new Zoom(o1, (n1-1)*d1 + o1, o2, (n2-1)*d2 + o2, 0, n1, 0, n2);

    int X, Y, Z;
    X = (int) n1/4;
    Y = (int) n2/4;
    Z = (int) n3/4;
    Slider1->SetRange(0,n3-1);
    Slider1->SetValue(Z);
    Slider1->SetTickFreq(Z);

    Slider2->SetRange(0,n1-1);
    Slider2->SetValue(X);
    Slider2->SetTickFreq(X);

    Slider3->SetRange(0,n2-1);
    Slider3->SetValue(Y);
    Slider3->SetTickFreq(Y);

    // Create image 
    this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
    this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
    this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());

}

void Image3dframe::OnXplotwindowPaint(wxPaintEvent& event)
{
    int w,h;
    wxPoint pos;
    wxPaintDC dc( Xplotwindow );

    Xplotwindow->GetSize(&w,&h);
    wxImage image_scaled=this->ximage->Scale(w,h);
    wxBitmap bitmap(image_scaled);
    dc.DrawBitmap(bitmap, 0, 0, false);
}

void Image3dframe::OnYplotwindowPaint(wxPaintEvent& event)
{
    int w,h;
    wxPoint pos;
    wxPaintDC dc( Yplotwindow );

    Xplotwindow->GetSize(&w,&h);
    wxImage image_scaled=this->yimage->Scale(w,h);
    wxBitmap bitmap(image_scaled);
    dc.DrawBitmap(bitmap, 0, 0, false);
}

void Image3dframe::OnZplotwindowPaint(wxPaintEvent& event)
{
    int w,h;
    wxPoint pos;
    wxPaintDC dc( Zplotwindow );

    Xplotwindow->GetSize(&w,&h);
    wxImage image_scaled=this->zimage->Scale(w,h);
    wxBitmap bitmap(image_scaled);
    dc.DrawBitmap(bitmap, 0, 0, false);
}

void Image3dframe::OnXaxisPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Xaxis );
    int w,h;

    float t0, t1;
    t0=xslicezoom->Getx0();
    t1=xslicezoom->Getx1();

    Xaxis->GetSize(&w,&h);

    float a = (t1-t0)/(w-1);
    float t;
    int tl;
    int i,n;
    float d;

    n=24;
    d=(t1-t0)/(n-1);
    char label[32];
    wxFont font(FONTSIZE1, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font);
    wxPoint pt;
    //Small ticks
    for(i=0; i<25; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        dc.SetPen(wxPen(*wxBLACK, 1, wxSOLID));
        dc.DrawLine(tl, h-4, tl, h-1);
    }
    wxCoord ht,wt;

    //Big ticks and text
    n=4;
    d=(t1-t0)/(n-1);
    for(i=1; i<4; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);

        snprintf(label, 32, "%.1f", t);
        dc.GetTextExtent(_(label), &wt, &ht);
        pt.y=h-5-ht;
        pt.x=tl-wt/2;
        dc.DrawText(_(label), pt);
        dc.SetPen(wxPen(*wxBLACK, 2, wxSOLID));
        dc.DrawLine(tl, h-5, tl, h-1);
    }
    wxFont font2(FONTSIZE2, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font2);
    dc.GetTextExtent(wxT("X"), &wt, &ht);
    dc.DrawRotatedText(wxT("X"), w/2-wt/2, 0, 0);
}

void Image3dframe::OnZaxisPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Zaxis );

    float t0, t1;
    t0=xslicezoom->Gety0();
    t1=xslicezoom->Gety1();
    int w,h;
    Zaxis->GetSize(&w,&h);

    float a = (t1-t0)/(h-1);
    float t;
    int tl;
    int i,n;
    float d;
    
    n=51;
    d=(t1-t0)/(n-1);
    char label[32];
    wxFont font(FONTSIZE1, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font);
    wxPoint pt;
    //Small ticks
    for(i=0; i<n; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        dc.SetPen(wxPen(*wxBLACK, 1, wxSOLID));
        dc.DrawLine(w-4, tl, w-1, tl);
    }
    
    n=11;
    d=(t1-t0)/(n-1);
    //Big ticks and text
    for(i=1; i<n-1; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        pt.x=w-30;
        pt.y=tl-5;
        snprintf(label, 32, "%.2f", t);
        dc.DrawText(_(label), pt);
        dc.SetPen(wxPen(*wxBLACK, 2, wxSOLID));
        dc.DrawLine(w-5, tl, w-1, tl);
    }
    wxCoord ht,wt;
    wxFont font2(FONTSIZE2, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font2);
    dc.GetTextExtent(wxT("Z"), &wt, &ht);
    dc.DrawRotatedText(wxT("Z"), 0, h/2+ht, 90);

}

void Image3dframe::OnYaxis_vertPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Yaxis_vert );

    float t0, t1;
    t0=yslicezoom->Getx0();
    t1=yslicezoom->Getx1();
    int w,h;
    Yaxis_vert->GetSize(&w,&h);

    float a = (t1-t0)/(h-1);
    float t;
    int tl;
    int i,n;
    float d;
    
    n=51;
    d=(t1-t0)/(n-1);
    char label[32];
    wxFont font(FONTSIZE1, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font);
    wxPoint pt;
    //Small ticks
    for(i=0; i<n; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        dc.SetPen(wxPen(*wxBLACK, 1, wxSOLID));
        dc.DrawLine(w-4, tl, w-1, tl);
    }
    
    n=11;
    d=(t1-t0)/(n-1);
    //Big ticks and text
    for(i=1; i<n-1; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        pt.x=w-30;
        pt.y=tl-5;
        snprintf(label, 32, "%.2f", t);
        dc.DrawText(_(label), pt);
        dc.SetPen(wxPen(*wxBLACK, 2, wxSOLID));
        dc.DrawLine(w-5, tl, w-1, tl);
    }
    wxCoord ht,wt;
    wxFont font2(FONTSIZE2, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font2);
    dc.GetTextExtent(wxT("Y"), &wt, &ht);
    dc.DrawRotatedText(wxT("Y"), 0, h/2+ht, 90);

}

void Image3dframe::OnYaxis_horPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Yaxis_hor );
    int w,h;

    float t0, t1;
    t0=yslicezoom->Getx0();
    t1=yslicezoom->Getx1();

    Yaxis_hor->GetSize(&w,&h);

    float a = (t1-t0)/(w-1);
    float t;
    int tl;
    int i,n;
    float d;

    n=24;
    d=(t1-t0)/(n-1);
    char label[32];
    wxFont font(FONTSIZE1, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font);
    wxPoint pt;
    //Small ticks
    for(i=0; i<25; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);
        dc.SetPen(wxPen(*wxBLACK, 1, wxSOLID));
        dc.DrawLine(tl, h-4, tl, h-1);
    }
    wxCoord ht,wt;

    //Big ticks and text
    n=4;
    d=(t1-t0)/(n-1);
    for(i=1; i<4; i++){
        t= t0 + i*d;
        tl = (int) ((t-t0)/a);

        snprintf(label, 32, "%.1f", t);
        dc.GetTextExtent(_(label), &wt, &ht);
        pt.y=h-5-ht;
        pt.x=tl-wt/2;
        dc.DrawText(_(label), pt);
        dc.SetPen(wxPen(*wxBLACK, 2, wxSOLID));
        dc.DrawLine(tl, h-5, tl, h-1);
    }
    wxFont font2(FONTSIZE2, wxFONTFAMILY_SWISS, wxNORMAL, wxBOLD);
    dc.SetFont(font2);
    dc.GetTextExtent(wxT("Y"), &wt, &ht);
    dc.DrawRotatedText(wxT("Y"), w/2-wt/2, 0, 0);

}

void Image3dframe::LoadXimage(int zi1, int zn1, int zi2, int zn2)
{
    if(getXimageAlloc()){
        delete(ximage);
    }
    ximage = new wxImage(zn1,zn2);
    setXimageAlloc(true);

    int Y = Slider3->GetValue();

    unsigned char *idata=ximage->GetData();
    int i,j;
    int val;
    float minclip = pclip[iminclip];
    float maxclip = pclip[imaxclip];

    if(minclip==maxclip) {
        minclip=0.0;
        maxclip=1.0;
    }

    for(j=0; j<zn2; j++){
        for (i=0; i<zn1; i++){
            val = (int) rintf(255*(imagedata[(j+zi2)*n1*n2 + Y*n1 + (i+zi1)]/(maxclip-minclip) - minclip/(maxclip-minclip)));
            //Clip data
            if( val < 0 )  val=0;
            if( val >255 ) val=255;
            idata[3*j*zn1 +3*i] = RGB[val*3 + 0];
            idata[3*j*zn1 +3*i +1] = RGB[val*3 +1];
            idata[3*j*zn1 +3*i +2] = RGB[val*3 +2];
        }
    }
}

void Image3dframe::LoadYimage(int zi1, int zn1, int zi2, int zn2)
{
    if(getYimageAlloc()){
        delete(yimage);
    }
    yimage = new wxImage(zn1,zn2);
    setYimageAlloc(true);

    int X = Slider2->GetValue();

    unsigned char *idata=yimage->GetData();
    int i,j;
    int val;
    float minclip = pclip[iminclip];
    float maxclip = pclip[imaxclip];

    if(minclip==maxclip) {
        minclip=0.0;
        maxclip=1.0;
    }

    for(j=0; j<zn2; j++){
        for (i=0; i<zn1; i++){
            val = (int) rintf(255*(imagedata[(j+zi2)*n1*n2 + (i+zi1)*n1 + X]/(maxclip-minclip) - minclip/(maxclip-minclip)));
            //Clip data
            if( val < 0 )  val=0;
            if( val >255 ) val=255;
            idata[3*j*zn1 +3*i] = RGB[val*3 + 0];
            idata[3*j*zn1 +3*i +1] = RGB[val*3 +1];
            idata[3*j*zn1 +3*i +2] = RGB[val*3 +2];
        }
    }
}

void Image3dframe::LoadZimage(int zi1, int zn1, int zi2, int zn2)
{
    if(getZimageAlloc()){
        delete(zimage);
    }
    zimage = new wxImage(zn1,zn2);
    setZimageAlloc(true);
    
    int Z = Slider1->GetValue();

    unsigned char *idata=zimage->GetData();
    int i,j;
    int val;
    float minclip = pclip[iminclip];
    float maxclip = pclip[imaxclip];

    if(minclip==maxclip) {
        minclip=0.0;
        maxclip=1.0;
    }

    for(j=0; j<zn2; j++){
        for (i=0; i<zn1; i++){
            val = (int) rintf(255*(imagedata[Z*n1*n2 + (j+zi2)*n1 + (i+zi1)]/(maxclip-minclip) - minclip/(maxclip-minclip)));
            //Clip data
            if( val < 0 )  val=0;
            if( val >255 ) val=255;
            idata[3*j*zn1 +3*i] = RGB[val*3 + 0];
            idata[3*j*zn1 +3*i +1] = RGB[val*3 +1];
            idata[3*j*zn1 +3*i +2] = RGB[val*3 +2];
        }
    }
}


void Image3dframe::ComputeClip()
{
    float *wrk;
    int ind;
    wrk = (float *) calloc(n1*n2*n3,sizeof(float)); 
    for(size_t i=0; i < n1*n2*n3; i++){
       wrk[i] = imagedata[i];
    }
    std::sort(wrk, wrk+(n1*n2*n3));
    if(n1*n2*n3 < 100){
        for(size_t i=0; i < n1*n2*n3; i++){
            pclip[i] = wrk[i];
        }
        iminclip=0; 
        imaxclip=(n1*n2*n3)-1;
    }else{
        for(size_t i=0; i < 100; i++){
            ind = (int) i*((n1*n2*n3)-1)/99;
            pclip[i] = wrk[ind];
        }
        iminclip=0; 
        imaxclip=99;
    }
    free(wrk);
}

void Image3dframe::OnSlider1CmdScroll(wxScrollEvent& event)
{
    this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());
    Refresh();
}

void Image3dframe::OnSlider2CmdScroll(wxScrollEvent& event)
{
    this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
    Refresh();
}

void Image3dframe::OnSlider3CmdScroll(wxScrollEvent& event)
{
    this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
    Refresh();
}

void Image3dframe::OnKeyUp(wxKeyEvent& event)
{

    if(event.GetKeyCode() == wxKeyCode('X')){
        // Increase clip
        iminclip++;
        imaxclip--;
        if(iminclip == imaxclip){ 
            iminclip--; 
            imaxclip++;
        }
        this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
        this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
        this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());
        Refresh();
    }
    if(event.GetKeyCode() == wxKeyCode('Z')){
        // Decrease clip
        iminclip--;
        imaxclip++;
        if(iminclip < 0){ 
            iminclip = 0; 
        }

        if(imaxclip > 99){ 
            imaxclip = 99; 
        }
        this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
        this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
        this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('E')){
        color--; 
        color = color%NCOLORS;
        if(color < 0) color+=NCOLORS;
        getRgb(color, &RGB[0]);
        this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
        this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
        this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('R')){
        color++; 
        color = color%NCOLORS;
        getRgb(color, &RGB[0]);
        this->LoadXimage(xslicezoom->Getix0(), xslicezoom->Getnx(), xslicezoom->Getiy0(), xslicezoom->Getny());
        this->LoadYimage(yslicezoom->Getix0(), yslicezoom->Getnx(), yslicezoom->Getiy0(), yslicezoom->Getny());
        this->LoadZimage(zslicezoom->Getix0(), zslicezoom->Getnx(), zslicezoom->Getiy0(), zslicezoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('Q')){
        Close();
    }
}

Image3dframe::~Image3dframe()
{
    //(*Destroy(Image3dframe)
    //*)
}


