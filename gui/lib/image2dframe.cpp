#include "image2dframe.h"
#define MAXPOS 1024


//(*IdInit(Image2dframe)
const long Image2dframe::ID_PANEL1 = wxNewId();
const long Image2dframe::ID_PANEL2 = wxNewId();
const long Image2dframe::ID_PANEL3 = wxNewId();
const long Image2dframe::ID_PANEL4 = wxNewId();
const long Image2dframe::ID_SCROLLEDWINDOW1 = wxNewId();
const long Image2dframe::ID_STATUSBAR1 = wxNewId();
//*)

BEGIN_EVENT_TABLE(Image2dframe,wxFrame)
	//(*EventTable(Image2dframe)
	//*)
END_EVENT_TABLE()

Image2dframe::Image2dframe(size_t _n1, float _d1, float _o1, size_t _n2, float _d2, float _o2, float *_imagedata, wxWindow* parent,wxWindowID id,const wxPoint& pos,const wxSize& size)
{
    //Setting up variables
    n1 = _n1;
    d1 = _d1;
    o1 = _o1;

    n2 = _n2;
    d2 = _d2;
    o2 = _o2;

    imagedata = _imagedata;
    image2d_allocated = false;

	//(*Initialize(Image2dframe)
	wxFlexGridSizer* FlexGridSizer1;
	wxBoxSizer* BoxSizer1;

	Create(parent, id, _("2D image window"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE, _T("id"));
	SetClientSize(wxDefaultSize);
	Move(wxDefaultPosition);
	BoxSizer1 = new wxBoxSizer(wxHORIZONTAL);
	FlexGridSizer1 = new wxFlexGridSizer(2, 3, 0, 0);
	FlexGridSizer1->AddGrowableCol(1);
	FlexGridSizer1->AddGrowableRow(1);
	LeftCorner = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1->Add(LeftCorner, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Xaxis = new wxPanel(this, ID_PANEL2, wxDefaultPosition, wxSize(600,50), wxTAB_TRAVERSAL, _T("ID_PANEL2"));
	FlexGridSizer1->Add(Xaxis, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	RightCorner = new wxPanel(this, ID_PANEL3, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL3"));
	FlexGridSizer1->Add(RightCorner, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Zaxis = new wxPanel(this, ID_PANEL4, wxDefaultPosition, wxSize(50,600), wxTAB_TRAVERSAL, _T("ID_PANEL4"));
	FlexGridSizer1->Add(Zaxis, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	Imagewindow = new wxScrolledWindow(this, ID_SCROLLEDWINDOW1, wxDefaultPosition, wxSize(600,600), wxVSCROLL|wxHSCROLL|wxFULL_REPAINT_ON_RESIZE, _T("ID_SCROLLEDWINDOW1"));
	FlexGridSizer1->Add(Imagewindow, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	BoxSizer1->Add(FlexGridSizer1, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 1);
	SetSizer(BoxSizer1);
    StatusBar1 = new wxStatusBar(this, ID_STATUSBAR1, wxSIMPLE_BORDER, _T("ID_STATUSBAR1"));
    int __wxStatusBarWidths_1[1] = { -10 };
    int __wxStatusBarStyles_1[1] = { wxSB_NORMAL };
    StatusBar1->SetFieldsCount(1,__wxStatusBarWidths_1);
    StatusBar1->SetStatusStyles(1,__wxStatusBarStyles_1);
    SetStatusBar(StatusBar1);

    BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);

	Xaxis->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image2dframe::OnXaxisPaint,0,this);
	Zaxis->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image2dframe::OnZaxisPaint,0,this);
	Imagewindow->Connect(wxEVT_PAINT,(wxObjectEventFunction)&Image2dframe::OnImagewindowPaint,0,this);
	Imagewindow->Connect(wxEVT_ERASE_BACKGROUND,(wxObjectEventFunction)&Image2dframe::OnImagewindowEraseBackground,0,this);
	Imagewindow->Connect(wxEVT_KEY_UP,(wxObjectEventFunction)&Image2dframe::OnImagewindowKeyUp,0,this);
	Imagewindow->Connect(wxEVT_LEFT_DOWN,(wxObjectEventFunction)&Image2dframe::OnImagewindowLeftDown,0,this);
	Imagewindow->Connect(wxEVT_LEFT_UP,(wxObjectEventFunction)&Image2dframe::OnImagewindowLeftUp,0,this);
	Imagewindow->Connect(wxEVT_RIGHT_UP,(wxObjectEventFunction)&Image2dframe::OnImagewindowRightUp,0,this);
	Imagewindow->Connect(wxEVT_MOTION,(wxObjectEventFunction)&Image2dframe::OnImagewindowMouseMove,0,this);
	Connect(wxID_ANY,wxEVT_CLOSE_WINDOW,(wxObjectEventFunction)&Image2dframe::OnClose);
	//*)
    
    // Compute Maxmimum and minimum values
    this->ComputeClip();

    zoom = new Zoom(o1, (n1-1)*d1 + o1, o2, (n2-1)*d2 + o2, 0, n1, 0, n2);

    // Create image 
    this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
}

Image2dframe::~Image2dframe()
{
	//(*Destroy(Image2dframe)
	//*)
    if(getImage2dAlloc()){
        delete(image2d);
    }

}

void Image2dframe::getRgb(int val, unsigned char *RGB)
{
    float pi=3.141592653589793;
    RGB[0]=rintf(255*sqrt((cos(val*pi/2/255))*(cos(val*pi/2/255))));
    RGB[1]=rintf(255*sqrt((cos(pi/2 + val*pi/255))*(cos(pi/2 + val*pi/255))));
    RGB[2]=rintf(255*sqrt((cos(pi/2 + val*pi/2/255))*(cos(pi/2 + val*pi/2/255))));
}

void Image2dframe::LoadImage(int zi1, int zn1, int zi2, int zn2)
{
    if(getImage2dAlloc()){
        delete(image2d);
    }
    image2d = new wxImage(zn1,zn2);
    setImage2dAlloc(true);

    unsigned char *idata=image2d->GetData();
    int i,j;
    int val;
    float minclip = pclip[iminclip];
    float maxclip = pclip[imaxclip];

    if(minclip==maxclip) {
        minclip=0.0;
        maxclip=1.0;
    }

    unsigned char RGB[3];
    for(j=0; j<zn2; j++){
        for (i=0; i<zn1; i++){
            val=rintf(255*(imagedata[(j+zi2)*n1 + (i+zi1)]/(maxclip-minclip) - minclip/(maxclip-minclip)));
            //Clip data
            if( val < 0 )  val=0;
            if( val >255 ) val=255;
            this->getRgb(val, &RGB[0]);
            //Convert to unsigned char
            idata[3*j*zn1 +3*i]=(unsigned char) RGB[2];
            idata[3*j*zn1 +3*i +1]=(unsigned char) RGB[1];
            idata[3*j*zn1 +3*i +2]=(unsigned char) RGB[0];
        }
    }
}

void Image2dframe::ComputeClip()
{
    float *wrk;
    int ind;
    wrk = (float *) calloc(n1*n2,sizeof(float)); 
    for(size_t i=0; i < n1*n2; i++){
       wrk[i] = imagedata[i];
    }
    std::sort(wrk, wrk+(n1*n2));
    if(n1*n2 < 100){
        for(size_t i=0; i < n1*n2; i++){
            pclip[i] = wrk[i];
        }
        iminclip=0; 
        imaxclip=(n1*n2)-1;
    }else{
        for(size_t i=0; i < 100; i++){
            ind = (int) i*((n1*n2)-1)/99;
            pclip[i] = wrk[ind];
        }
        iminclip=0; 
        imaxclip=99;
    }
    free(wrk);
}

void Image2dframe::OnImagewindowPaint(wxPaintEvent& event)
{
    int w,h;
    wxPoint pos;
    wxPaintDC dc( Imagewindow );

    Imagewindow->GetSize(&w,&h);
    wxImage image_scaled=this->image2d->Scale(w,h);
    wxBitmap bitmap(image_scaled);
    dc.DrawBitmap(bitmap, 0, 0, false);

    //Plot Zoom box
    if(zoom->Getzooming()){
        wxPen myWhitePen(*wxWHITE,2,wxSOLID);
        dc.SetPen(myWhitePen);
        dc.DrawLines( 5, zoom->Getbox() );
    }

}

void Image2dframe::OnZaxisPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Zaxis );

    float t0, t1;
    t0=zoom->Gety0();
    t1=zoom->Gety1();
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

void Image2dframe::OnXaxisPaint(wxPaintEvent& event)
{
    wxPaintDC dc( Xaxis );
    int w,h;

    float t0, t1;
    t0=zoom->Getx0();
    t1=zoom->Getx1();

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

void Image2dframe::OnImagewindowMouseMove(wxMouseEvent& event)
{

    wxPaintDC dc( Imagewindow );
    wxPoint pos;
    wxPoint box[5];
    wxPoint zpos1=zoom->Getzpos1();
    pos=event.GetLogicalPosition(dc);


    float x0, x1, y0, y1;
    x0=zoom->Getx0();
    x1=zoom->Getx1();
    y0=zoom->Gety0();
    y1=zoom->Gety1();

    int w,h;
    Imagewindow->GetSize(&w,&h);
    float ax, ay;
    ay=(y1 - y0)/(h-1);
    ax=(x1 - x0)/(w-1);

    // Compute current position
    float x,y;
    x=ax*pos.x + x0;
    y=ay*pos.y + y0;

    // Get current value
    int ix,iy;
    float val;
    ix = (int) (x-x0)/d1;
    iy = (int) (y-y0)/d2;
    val = imagedata[iy*n1 + ix];

    char label[48];
    snprintf(label, 48, "X: %.2f, Z: %.2f, Val: %f", x, y, val);
    StatusBar1->SetStatusText(_(label));

    if(event.LeftIsDown()){
        box[0].x=zpos1.x;
        box[0].y=zpos1.y;
        zoom->Setbox(box[0], 0);
        box[1].x=zpos1.x;
        box[1].y=pos.y;
        zoom->Setbox(box[1], 1);
        box[2].x=pos.x;
        box[2].y=pos.y;
        zoom->Setbox(box[2], 2);
        box[3].x=pos.x;
        box[3].y=zpos1.y;
        zoom->Setbox(box[3], 3);
        box[4].x=zpos1.x;
        box[4].y=zpos1.y;
        zoom->Setbox(box[4], 4);
        Refresh();
    }
    event.Skip();
}

void Image2dframe::OnImagewindowLeftUp(wxMouseEvent& event)
{
    wxPoint pos;
    wxPaintDC dc( Imagewindow );
    pos=event.GetLogicalPosition(dc);
    int w,h;
    Imagewindow->GetSize(&w,&h);

    float x0, x1, y0, y1;
    x0=zoom->Getx0();
    x1=zoom->Getx1();
    y0=zoom->Gety0();
    y1=zoom->Gety1();

    
    float ax, ay;
    ay=(y1 - y0)/(h-1);
    ax=(x1 - x0)/(w-1);

    wxPoint zpos1=zoom->Getzpos1();
    wxPoint temp;
    if(zpos1.y > pos.y){
        temp.y=zpos1.y;
        zpos1.y=pos.y;
        pos.y=temp.y;
    }
    if(zpos1.x > pos.x){
        temp.x=zpos1.x;
        zpos1.x=pos.x;
        pos.x=temp.x;
    }
    if(zpos1.x<0) zpos1.x=0;
    if(zpos1.y<0) zpos1.y=0;
    if(zpos1.x>w-1) zpos1.x=w-1;
    if(zpos1.y>h-1) zpos1.y=h-1;
    if(pos.x<0) pos.x=0;
    if(pos.y<0) pos.y=0;
    if(pos.x>w-1) pos.x=w-1;
    if(pos.y>h-1) pos.y=h-1;
    zoom->Setzpos1(zpos1);
    zoom->Setzpos2(pos);
    if(zpos1.x != pos.x && zpos1.y != pos.y){
        //Compute limits
        y1=ay*pos.y + y0;
        y0=ay*zpos1.y + y0;
        x1=ax*pos.x + x0;
        x0=ax*zpos1.x + x0;
        zoom->Sety0(y0);
        zoom->Setiy0(int ((y0-o2)/d2));
        zoom->Setny(int ((y1-y0)/d2 + 1 ));
        zoom->Sety1(y1);
        zoom->Setx0(x0);
        zoom->Setix0(int ((x0-o1)/d1));
        zoom->Setnx(int ((x1-x0)/d1 + 1 ));
        zoom->Setx1(x1);
    }else{
        zoom->Sety0(0.0);
        zoom->Sety1((n2-1)*d2 + o2);
        zoom->Setx0(o1);
        zoom->Setx1((n1-1)*d1 + o1);
        zoom->Setiy0(0);
        zoom->Setny(n2);
        zoom->Setix0(0);
        zoom->Setnx(n1);
    }
    this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
    zoom->Setzooming(false);
    Refresh();
}

void Image2dframe::OnImagewindowRightUp(wxMouseEvent& event)
{
}

void Image2dframe::OnImagewindowLeftDown(wxMouseEvent& event)
{
    wxPaintDC dc( Imagewindow );
    wxPoint zpos1=event.GetLogicalPosition(dc);
    zoom->Setzpos1(zpos1);
    zoom->Setzooming(true);
}

void Image2dframe::OnImagewindowKeyUp(wxKeyEvent& event)
{
    if(event.GetKeyCode() == 88){
        // Increase clip
        iminclip++;
        imaxclip--;
        if(iminclip == imaxclip){ 
            iminclip--; 
            imaxclip++;
        }
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }
    if(event.GetKeyCode() == 90){
        // Decrease clip
        iminclip--;
        imaxclip++;
        if(iminclip < 0){ 
            iminclip = 0; 
        }

        if(imaxclip > 99){ 
            imaxclip = 99; 
        }
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }

}

void Image2dframe::OnClose(wxCloseEvent& event)
{
    event.Skip();  // Quit
}

void Image2dframe::OnImagewindowEraseBackground(wxEraseEvent& event)
{
    event.Skip();  // Quit
}
