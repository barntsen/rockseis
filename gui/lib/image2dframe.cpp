#include "image2dframe.h"

//helper functions
enum wxbuildinfoformat {
    short_f, long_f };

wxString wxbuildinfo(wxbuildinfoformat format)
{
    wxString wxbuild(wxVERSION_STRING);

    if (format == long_f )
    {
#if defined(__WXMSW__)
        wxbuild << _T("-Windows");
#elif defined(__UNIX__)
        wxbuild << _T("-Linux");
#endif

#if wxUSE_UNICODE
        wxbuild << _T("-Unicode build");
#else
        wxbuild << _T("-ANSI build");
#endif // wxUSE_UNICODE
    }

    return wxbuild;
}

//(*IdInit(Image2dframe)
const long Image2dframe::ID_PANEL1 = wxNewId();
const long Image2dframe::ID_PANEL2 = wxNewId();
const long Image2dframe::ID_PANEL3 = wxNewId();
const long Image2dframe::ID_PANEL4 = wxNewId();
const long Image2dframe::ID_IMAGEWINDOW1 = wxNewId();
const long Image2dframe::ID_STATUSBAR1 = wxNewId();
const long Image2dframe::idToolNext = wxNewId();
const long Image2dframe::idToolcmpint = wxNewId();
const long Image2dframe::idToolPrev = wxNewId();
const long Image2dframe::idToolpick = wxNewId();
const long Image2dframe::idToolzoom = wxNewId();
const long Image2dframe::idToolSave = wxNewId();
const long Image2dframe::idToolLoad = wxNewId();
const long Image2dframe::ID_TOOLBAR1 = wxNewId();
const long Image2dframe::ID_LISTBOX1 = wxNewId();
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
    toolbarset = false;

    cmpnumber = 0;
    dcmp = 1;
    maxcmp = 1;

    nlayers=0;
    layer=0;

    displaycrosshair = false;
    getcrosshair = false;

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
	Imagewindow = new wxPanel(this, ID_IMAGEWINDOW1, wxDefaultPosition, wxSize(600,600), wxWANTS_CHARS, _T("ID_IMAGEWINDOW1"));
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
	Connect(idToolpick,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnPick);
	Connect(idToolzoom,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnZoom);
	//*)
    
    // Compute Maxmimum and minimum values
    this->ComputeClip();

    // Setup color pallete
    color = 0;
    this->getRgb(color);
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

    for(j=0; j<zn2; j++){
        for (i=0; i<zn1; i++){
            val = (int) rintf(255*(imagedata[(j+zi2)*n1 + (i+zi1)]/(maxclip-minclip) - minclip/(maxclip-minclip)));
            //Clip data
            if( val < 0 )  val=0;
            if( val >255 ) val=255;
            idata[3*j*zn1 +3*i] = RGB[val*3 + 0];
            idata[3*j*zn1 +3*i +1] = RGB[val*3 +1];
            idata[3*j*zn1 +3*i +2] = RGB[val*3 +2];
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

    if(nlayers>0){
        //Plot picks
        Plotpicks(dc, w, h);
    }


    if(displaycrosshair){
        //Plot crosshair
        Plotcrosshair(dc, w, h);
    }

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
    Imagewindow->SetFocus();

    bool dozoom = false;
    if(!toolbarset){
        dozoom = true;
    }
    if(toolbarset){
        if(ToolBarItem4->IsToggled()){
            dozoom = true;
        }
    }

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
    size_t ix,iy;
    float val;
    ix = (size_t) ((x-o1)/d1);
    iy = (size_t) ((y-o2)/d2);
    if(ix < n1 && iy < n2){
	    val = imagedata[iy*n1 + ix];
    }else{
	    val = 0.0;
    }

    char label[48];
    if(getToolbarset()){
        snprintf(label, 48, "CMP: %d, X: %.2f, Z: %.2f, Val: %f", this->getCmpnumber(), x, y, val);
    }else{
        snprintf(label, 48, "CMP: %zu, X: %.2f, Z: %.2f, Val: %f", ix, x, y, val);
    }
    StatusBar1->SetStatusText(_(label));

    if(event.LeftIsDown() && dozoom){
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

    if(getcrosshair){
        if(getToolbarset()){
            crosshair_pt[0] = this->getCmpnumber();
            crosshair_pt[1] = y;
        }else{
            crosshair_pt[0] = x;
            crosshair_pt[1] = y;
        }
        wxCommandEvent parevent(Crosshair, GetId());
        parevent.SetEventObject(this);
        parevent.SetClientData((void*) &crosshair_pt[0]);
        // Send event to App
        ProcessWindowEvent(parevent);
    }
    event.Skip(true);
}

void Image2dframe::OnImagewindowLeftUp(wxMouseEvent& event)
{
    wxPoint pos;
    wxPaintDC dc( Imagewindow );
    pos=event.GetLogicalPosition(dc);
    int w,h;
    Imagewindow->GetSize(&w,&h);
    bool dozoom = false;


    float x0, x1, y0, y1;
    x0=zoom->Getx0();
    x1=zoom->Getx1();
    y0=zoom->Gety0();
    y1=zoom->Gety1();
    
    float ax, ay;
    ay=(y1 - y0)/(h-1);
    ax=(x1 - x0)/(w-1);

    float x,y;
    x=ax*pos.x + x0;
    y=ay*pos.y + y0;

    if(toolbarset){
        if(ToolBarItem3->IsToggled() && nlayers > 0){
            picks[layer]->Addpick(this->getCmpnumber(), x, y);
            Refresh();
        }
    }
    if(!toolbarset){
        dozoom = true;
    }
    if(toolbarset){
        if(ToolBarItem4->IsToggled()){
            dozoom = true;
        }
    }


    if(dozoom){
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
        zoom->Setzooming(false);
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }
}

void Image2dframe::OnImagewindowRightUp(wxMouseEvent& event)
{
    if(toolbarset){
        if(ToolBarItem3->IsToggled() && nlayers > 0){
            int w,h;
            float x0, x1, y0, y1;
            float ax;
            float ay;
            float x,y;

            wxPaintDC dc( Imagewindow );
            wxPoint pos;
            pos=event.GetLogicalPosition(dc);
            Imagewindow->GetSize(&w,&h);

            x0=zoom->Getx0();
            x1=zoom->Getx1();
            y0=zoom->Gety0();
            y1=zoom->Gety1();

            ay=(y1 - y0)/(h-1);
            ax=(x1 - x0)/(w-1);

            x=ax*pos.x + x0;
            y=ay*pos.y + y0;

            picks[layer]->Removepick(this->getCmpnumber(), x, y);
            Refresh();
        }
    }
}

void Image2dframe::OnImagewindowLeftDown(wxMouseEvent& event)
{
    wxPaintDC dc( Imagewindow );
    wxPoint zpos1=event.GetLogicalPosition(dc);
    bool dozoom = false;
    if(!toolbarset){
        dozoom = true;
    }
    if(toolbarset){
        if(ToolBarItem4->IsToggled()){
            dozoom = true;
        }
    }
    if(dozoom){
        zoom->Setzpos1(zpos1);
        zoom->Setzooming(true);
    }
}

void Image2dframe::OnImagewindowKeyUp(wxKeyEvent& event)
{

    if(event.GetKeyCode() == wxKeyCode('X')){
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
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('E')){
        color--; 
        color = color%NCOLORS;
        if(color < 0) color+=NCOLORS;
        this->getRgb(color);
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('R')){
        color++; 
        color = color%NCOLORS;
        this->getRgb(color);
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('Q')){
        Close();
    }
}

void Image2dframe::createToolbar()
{
    if(!this->getToolbarset()){
        ToolBar1 = new wxToolBar(this, ID_TOOLBAR1, wxDefaultPosition, wxDefaultSize, wxTB_VERTICAL|wxTB_HORIZONTAL|wxNO_BORDER, _T("ID_TOOLBAR1"));
        ToolBarItem1 = ToolBar1->AddTool(idToolNext, _("Next"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_GO_FORWARD")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, _("Jump to next CMP position"), wxEmptyString);
        ToolBarItem7 = ToolBar1->AddTool(idToolcmpint, _("CMP interval"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_REDO")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, _("Change CMP interval"), wxEmptyString);
        ToolBarItem2 = ToolBar1->AddTool(idToolPrev, _("Previous"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_GO_BACK")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, _("Jump to previous CMP position "), wxEmptyString);
        ToolBarItem3 = ToolBar1->AddTool(idToolpick, _("Pick"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_ADD_BOOKMARK")),wxART_TOOLBAR), wxNullBitmap, wxITEM_CHECK, _("Pick velocities"), wxEmptyString);
        ToolBarItem4 = ToolBar1->AddTool(idToolzoom, _("Zoom"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_FIND")),wxART_TOOLBAR), wxNullBitmap, wxITEM_CHECK, wxEmptyString, wxEmptyString);
        ToolBarItem5 = ToolBar1->AddTool(idToolSave, _("Save picks"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_FILE_SAVE")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, _("Save picks"), wxEmptyString);
        ToolBarItem6 = ToolBar1->AddTool(idToolLoad, _("Load picks"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_FILE_OPEN")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, _("Load picks"), wxEmptyString);
        SetToolBar(ToolBar1);
        ToolBar1->Realize();
        this->setToolbarset(true);
        Layout();

        Connect(idToolNext,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnNextClicked);
        Connect(idToolPrev,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnPreviousClicked);
        Connect(idToolcmpint,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnCMPinterval);
	    Connect(idToolSave,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnSavepicks);
    }
}

void Image2dframe::destroyToolbar()
{
    if(this->getToolbarset()){
        delete GetToolBar();
        SetToolBar(0);
        this->setToolbarset(false);
        Layout();
    }
}

void Image2dframe::OnClose(wxCloseEvent& event)
{
    event.Skip(true);
    Destroy();
}

void Image2dframe::OnImagewindowEraseBackground(wxEraseEvent& event)
{
}

void Image2dframe::getRgb(int which)
{
    float n = 64.;
    float u = 0.0;
    float r,g,b;

    for (int i=0; i<256; i++)
    {
        RGB[i*3] = (unsigned char) 0;
        RGB[i*3 + 1] = (unsigned char) 0;
        RGB[i*3 + 2] = (unsigned char) 0;
    }

    switch(which){
        case 0:
            for (int i=0; i<256; i++)
            {
                RGB[i*3] = (unsigned char) i;
                RGB[i*3 + 1] = (unsigned char) i;
                RGB[i*3 + 2] = (unsigned char) i;
            }
            break;
        case 1:
            for (int i=0; i<256; i++)
            {
                if(i<n) u =  (i)/n;
                if(i >=n && i < 2*n ) u = 1;
                if(i >=2*n && i < 3*n ) u = 1 - (i-2*n)/n;
                if(i >= 3*n) u = 0;
                if(i < 192) {
                    g = n/2. + i; 
                }else{
                    g = 0.0;
                }
                r = g + n;
                b = g - n;

                if(r >= 0 && r < 256) RGB[((int) r)*3] = (unsigned char) (255*u);
                if(g >= 0 && g < 256) RGB[((int) g)*3 + 1] = (unsigned char) (255*u);
                if(b >= 0 && b < 256) RGB[((int) b)*3 + 2] = (unsigned char) (255*u);
            }
            break;
        case 2:
            RGB[0*3] = 53;
            RGB[0*3 +1] = 42;
            RGB[0*3 +2] = 135;
            RGB[1*3] = 53;
            RGB[1*3 +1] = 44;
            RGB[1*3 +2] = 138;
            RGB[2*3] = 54;
            RGB[2*3 +1] = 45;
            RGB[2*3 +2] = 141;
            RGB[3*3] = 54;
            RGB[3*3 +1] = 47;
            RGB[3*3 +2] = 144;
            RGB[4*3] = 54;
            RGB[4*3 +1] = 48;
            RGB[4*3 +2] = 147;
            RGB[5*3] = 54;
            RGB[5*3 +1] = 50;
            RGB[5*3 +2] = 150;
            RGB[6*3] = 54;
            RGB[6*3 +1] = 51;
            RGB[6*3 +2] = 153;
            RGB[7*3] = 54;
            RGB[7*3 +1] = 53;
            RGB[7*3 +2] = 156;
            RGB[8*3] = 54;
            RGB[8*3 +1] = 54;
            RGB[8*3 +2] = 160;
            RGB[9*3] = 54;
            RGB[9*3 +1] = 56;
            RGB[9*3 +2] = 163;
            RGB[10*3] = 54;
            RGB[10*3 +1] = 57;
            RGB[10*3 +2] = 166;
            RGB[11*3] = 54;
            RGB[11*3 +1] = 59;
            RGB[11*3 +2] = 169;
            RGB[12*3] = 53;
            RGB[12*3 +1] = 61;
            RGB[12*3 +2] = 172;
            RGB[13*3] = 53;
            RGB[13*3 +1] = 62;
            RGB[13*3 +2] = 175;
            RGB[14*3] = 52;
            RGB[14*3 +1] = 64;
            RGB[14*3 +2] = 179;
            RGB[15*3] = 51;
            RGB[15*3 +1] = 65;
            RGB[15*3 +2] = 182;
            RGB[16*3] = 50;
            RGB[16*3 +1] = 67;
            RGB[16*3 +2] = 185;
            RGB[17*3] = 49;
            RGB[17*3 +1] = 69;
            RGB[17*3 +2] = 188;
            RGB[18*3] = 48;
            RGB[18*3 +1] = 70;
            RGB[18*3 +2] = 191;
            RGB[19*3] = 46;
            RGB[19*3 +1] = 72;
            RGB[19*3 +2] = 195;
            RGB[20*3] = 44;
            RGB[20*3 +1] = 74;
            RGB[20*3 +2] = 198;
            RGB[21*3] = 42;
            RGB[21*3 +1] = 76;
            RGB[21*3 +2] = 201;
            RGB[22*3] = 39;
            RGB[22*3 +1] = 78;
            RGB[22*3 +2] = 204;
            RGB[23*3] = 36;
            RGB[23*3 +1] = 80;
            RGB[23*3 +2] = 208;
            RGB[24*3] = 33;
            RGB[24*3 +1] = 82;
            RGB[24*3 +2] = 211;
            RGB[25*3] = 29;
            RGB[25*3 +1] = 84;
            RGB[25*3 +2] = 214;
            RGB[26*3] = 25;
            RGB[26*3 +1] = 87;
            RGB[26*3 +2] = 217;
            RGB[27*3] = 21;
            RGB[27*3 +1] = 89;
            RGB[27*3 +2] = 219;
            RGB[28*3] = 16;
            RGB[28*3 +1] = 91;
            RGB[28*3 +2] = 221;
            RGB[29*3] = 12;
            RGB[29*3 +1] = 93;
            RGB[29*3 +2] = 222;
            RGB[30*3] = 8;
            RGB[30*3 +1] = 95;
            RGB[30*3 +2] = 224;
            RGB[31*3] = 5;
            RGB[31*3 +1] = 97;
            RGB[31*3 +2] = 224;
            RGB[32*3] = 3;
            RGB[32*3 +1] = 98;
            RGB[32*3 +2] = 225;
            RGB[33*3] = 2;
            RGB[33*3 +1] = 100;
            RGB[33*3 +2] = 225;
            RGB[34*3] = 2;
            RGB[34*3 +1] = 101;
            RGB[34*3 +2] = 225;
            RGB[35*3] = 1;
            RGB[35*3 +1] = 102;
            RGB[35*3 +2] = 225;
            RGB[36*3] = 1;
            RGB[36*3 +1] = 104;
            RGB[36*3 +2] = 225;
            RGB[37*3] = 2;
            RGB[37*3 +1] = 105;
            RGB[37*3 +2] = 225;
            RGB[38*3] = 2;
            RGB[38*3 +1] = 106;
            RGB[38*3 +2] = 225;
            RGB[39*3] = 3;
            RGB[39*3 +1] = 107;
            RGB[39*3 +2] = 225;
            RGB[40*3] = 4;
            RGB[40*3 +1] = 108;
            RGB[40*3 +2] = 224;
            RGB[41*3] = 5;
            RGB[41*3 +1] = 109;
            RGB[41*3 +2] = 224;
            RGB[42*3] = 6;
            RGB[42*3 +1] = 110;
            RGB[42*3 +2] = 223;
            RGB[43*3] = 7;
            RGB[43*3 +1] = 111;
            RGB[43*3 +2] = 223;
            RGB[44*3] = 8;
            RGB[44*3 +1] = 112;
            RGB[44*3 +2] = 223;
            RGB[45*3] = 9;
            RGB[45*3 +1] = 113;
            RGB[45*3 +2] = 222;
            RGB[46*3] = 10;
            RGB[46*3 +1] = 114;
            RGB[46*3 +2] = 222;
            RGB[47*3] = 11;
            RGB[47*3 +1] = 115;
            RGB[47*3 +2] = 221;
            RGB[48*3] = 12;
            RGB[48*3 +1] = 116;
            RGB[48*3 +2] = 221;
            RGB[49*3] = 13;
            RGB[49*3 +1] = 117;
            RGB[49*3 +2] = 220;
            RGB[50*3] = 14;
            RGB[50*3 +1] = 118;
            RGB[50*3 +2] = 220;
            RGB[51*3] = 15;
            RGB[51*3 +1] = 119;
            RGB[51*3 +2] = 219;
            RGB[52*3] = 16;
            RGB[52*3 +1] = 120;
            RGB[52*3 +2] = 218;
            RGB[53*3] = 16;
            RGB[53*3 +1] = 121;
            RGB[53*3 +2] = 218;
            RGB[54*3] = 17;
            RGB[54*3 +1] = 122;
            RGB[54*3 +2] = 217;
            RGB[55*3] = 18;
            RGB[55*3 +1] = 123;
            RGB[55*3 +2] = 217;
            RGB[56*3] = 18;
            RGB[56*3 +1] = 124;
            RGB[56*3 +2] = 216;
            RGB[57*3] = 19;
            RGB[57*3 +1] = 125;
            RGB[57*3 +2] = 216;
            RGB[58*3] = 19;
            RGB[58*3 +1] = 126;
            RGB[58*3 +2] = 215;
            RGB[59*3] = 19;
            RGB[59*3 +1] = 127;
            RGB[59*3 +2] = 215;
            RGB[60*3] = 20;
            RGB[60*3 +1] = 128;
            RGB[60*3 +2] = 214;
            RGB[61*3] = 20;
            RGB[61*3 +1] = 129;
            RGB[61*3 +2] = 214;
            RGB[62*3] = 20;
            RGB[62*3 +1] = 130;
            RGB[62*3 +2] = 213;
            RGB[63*3] = 20;
            RGB[63*3 +1] = 131;
            RGB[63*3 +2] = 213;
            RGB[64*3] = 20;
            RGB[64*3 +1] = 132;
            RGB[64*3 +2] = 212;
            RGB[65*3] = 20;
            RGB[65*3 +1] = 133;
            RGB[65*3 +2] = 212;
            RGB[66*3] = 20;
            RGB[66*3 +1] = 134;
            RGB[66*3 +2] = 211;
            RGB[67*3] = 20;
            RGB[67*3 +1] = 135;
            RGB[67*3 +2] = 211;
            RGB[68*3] = 19;
            RGB[68*3 +1] = 136;
            RGB[68*3 +2] = 211;
            RGB[69*3] = 19;
            RGB[69*3 +1] = 137;
            RGB[69*3 +2] = 211;
            RGB[70*3] = 18;
            RGB[70*3 +1] = 138;
            RGB[70*3 +2] = 210;
            RGB[71*3] = 18;
            RGB[71*3 +1] = 140;
            RGB[71*3 +2] = 210;
            RGB[72*3] = 17;
            RGB[72*3 +1] = 141;
            RGB[72*3 +2] = 210;
            RGB[73*3] = 16;
            RGB[73*3 +1] = 142;
            RGB[73*3 +2] = 210;
            RGB[74*3] = 15;
            RGB[74*3 +1] = 143;
            RGB[74*3 +2] = 210;
            RGB[75*3] = 14;
            RGB[75*3 +1] = 145;
            RGB[75*3 +2] = 210;
            RGB[76*3] = 13;
            RGB[76*3 +1] = 146;
            RGB[76*3 +2] = 210;
            RGB[77*3] = 12;
            RGB[77*3 +1] = 147;
            RGB[77*3 +2] = 210;
            RGB[78*3] = 11;
            RGB[78*3 +1] = 149;
            RGB[78*3 +2] = 210;
            RGB[79*3] = 10;
            RGB[79*3 +1] = 150;
            RGB[79*3 +2] = 210;
            RGB[80*3] = 9;
            RGB[80*3 +1] = 151;
            RGB[80*3 +2] = 209;
            RGB[81*3] = 9;
            RGB[81*3 +1] = 152;
            RGB[81*3 +2] = 209;
            RGB[82*3] = 8;
            RGB[82*3 +1] = 153;
            RGB[82*3 +2] = 209;
            RGB[83*3] = 8;
            RGB[83*3 +1] = 154;
            RGB[83*3 +2] = 208;
            RGB[84*3] = 7;
            RGB[84*3 +1] = 155;
            RGB[84*3 +2] = 208;
            RGB[85*3] = 7;
            RGB[85*3 +1] = 156;
            RGB[85*3 +2] = 207;
            RGB[86*3] = 7;
            RGB[86*3 +1] = 157;
            RGB[86*3 +2] = 207;
            RGB[87*3] = 6;
            RGB[87*3 +1] = 158;
            RGB[87*3 +2] = 206;
            RGB[88*3] = 6;
            RGB[88*3 +1] = 159;
            RGB[88*3 +2] = 206;
            RGB[89*3] = 6;
            RGB[89*3 +1] = 160;
            RGB[89*3 +2] = 205;
            RGB[90*3] = 6;
            RGB[90*3 +1] = 161;
            RGB[90*3 +2] = 204;
            RGB[91*3] = 6;
            RGB[91*3 +1] = 162;
            RGB[91*3 +2] = 203;
            RGB[92*3] = 6;
            RGB[92*3 +1] = 163;
            RGB[92*3 +2] = 203;
            RGB[93*3] = 6;
            RGB[93*3 +1] = 164;
            RGB[93*3 +2] = 202;
            RGB[94*3] = 6;
            RGB[94*3 +1] = 164;
            RGB[94*3 +2] = 201;
            RGB[95*3] = 6;
            RGB[95*3 +1] = 165;
            RGB[95*3 +2] = 200;
            RGB[96*3] = 6;
            RGB[96*3 +1] = 166;
            RGB[96*3 +2] = 199;
            RGB[97*3] = 6;
            RGB[97*3 +1] = 167;
            RGB[97*3 +2] = 198;
            RGB[98*3] = 6;
            RGB[98*3 +1] = 167;
            RGB[98*3 +2] = 197;
            RGB[99*3] = 6;
            RGB[99*3 +1] = 168;
            RGB[99*3 +2] = 196;
            RGB[100*3] = 6;
            RGB[100*3 +1] = 169;
            RGB[100*3 +2] = 195;
            RGB[101*3] = 7;
            RGB[101*3 +1] = 169;
            RGB[101*3 +2] = 194;
            RGB[102*3] = 7;
            RGB[102*3 +1] = 170;
            RGB[102*3 +2] = 193;
            RGB[103*3] = 8;
            RGB[103*3 +1] = 171;
            RGB[103*3 +2] = 192;
            RGB[104*3] = 9;
            RGB[104*3 +1] = 171;
            RGB[104*3 +2] = 191;
            RGB[105*3] = 10;
            RGB[105*3 +1] = 172;
            RGB[105*3 +2] = 190;
            RGB[106*3] = 11;
            RGB[106*3 +1] = 172;
            RGB[106*3 +2] = 189;
            RGB[107*3] = 12;
            RGB[107*3 +1] = 173;
            RGB[107*3 +2] = 188;
            RGB[108*3] = 13;
            RGB[108*3 +1] = 174;
            RGB[108*3 +2] = 186;
            RGB[109*3] = 15;
            RGB[109*3 +1] = 174;
            RGB[109*3 +2] = 185;
            RGB[110*3] = 16;
            RGB[110*3 +1] = 175;
            RGB[110*3 +2] = 184;
            RGB[111*3] = 18;
            RGB[111*3 +1] = 175;
            RGB[111*3 +2] = 183;
            RGB[112*3] = 19;
            RGB[112*3 +1] = 176;
            RGB[112*3 +2] = 182;
            RGB[113*3] = 21;
            RGB[113*3 +1] = 176;
            RGB[113*3 +2] = 180;
            RGB[114*3] = 23;
            RGB[114*3 +1] = 177;
            RGB[114*3 +2] = 179;
            RGB[115*3] = 24;
            RGB[115*3 +1] = 178;
            RGB[115*3 +2] = 178;
            RGB[116*3] = 26;
            RGB[116*3 +1] = 178;
            RGB[116*3 +2] = 177;
            RGB[117*3] = 28;
            RGB[117*3 +1] = 179;
            RGB[117*3 +2] = 175;
            RGB[118*3] = 30;
            RGB[118*3 +1] = 179;
            RGB[118*3 +2] = 174;
            RGB[119*3] = 32;
            RGB[119*3 +1] = 180;
            RGB[119*3 +2] = 173;
            RGB[120*3] = 34;
            RGB[120*3 +1] = 180;
            RGB[120*3 +2] = 171;
            RGB[121*3] = 36;
            RGB[121*3 +1] = 181;
            RGB[121*3 +2] = 170;
            RGB[122*3] = 38;
            RGB[122*3 +1] = 181;
            RGB[122*3 +2] = 169;
            RGB[123*3] = 40;
            RGB[123*3 +1] = 182;
            RGB[123*3 +2] = 167;
            RGB[124*3] = 43;
            RGB[124*3 +1] = 182;
            RGB[124*3 +2] = 166;
            RGB[125*3] = 45;
            RGB[125*3 +1] = 183;
            RGB[125*3 +2] = 165;
            RGB[126*3] = 47;
            RGB[126*3 +1] = 183;
            RGB[126*3 +2] = 163;
            RGB[127*3] = 49;
            RGB[127*3 +1] = 184;
            RGB[127*3 +2] = 162;
            RGB[128*3] = 52;
            RGB[128*3 +1] = 184;
            RGB[128*3 +2] = 160;
            RGB[129*3] = 54;
            RGB[129*3 +1] = 185;
            RGB[129*3 +2] = 159;
            RGB[130*3] = 57;
            RGB[130*3 +1] = 185;
            RGB[130*3 +2] = 157;
            RGB[131*3] = 59;
            RGB[131*3 +1] = 186;
            RGB[131*3 +2] = 156;
            RGB[132*3] = 62;
            RGB[132*3 +1] = 186;
            RGB[132*3 +2] = 154;
            RGB[133*3] = 64;
            RGB[133*3 +1] = 186;
            RGB[133*3 +2] = 153;
            RGB[134*3] = 67;
            RGB[134*3 +1] = 187;
            RGB[134*3 +2] = 151;
            RGB[135*3] = 70;
            RGB[135*3 +1] = 187;
            RGB[135*3 +2] = 150;
            RGB[136*3] = 73;
            RGB[136*3 +1] = 188;
            RGB[136*3 +2] = 148;
            RGB[137*3] = 75;
            RGB[137*3 +1] = 188;
            RGB[137*3 +2] = 147;
            RGB[138*3] = 78;
            RGB[138*3 +1] = 188;
            RGB[138*3 +2] = 145;
            RGB[139*3] = 81;
            RGB[139*3 +1] = 189;
            RGB[139*3 +2] = 144;
            RGB[140*3] = 84;
            RGB[140*3 +1] = 189;
            RGB[140*3 +2] = 142;
            RGB[141*3] = 87;
            RGB[141*3 +1] = 189;
            RGB[141*3 +2] = 141;
            RGB[142*3] = 90;
            RGB[142*3 +1] = 189;
            RGB[142*3 +2] = 139;
            RGB[143*3] = 93;
            RGB[143*3 +1] = 190;
            RGB[143*3 +2] = 138;
            RGB[144*3] = 96;
            RGB[144*3 +1] = 190;
            RGB[144*3 +2] = 136;
            RGB[145*3] = 99;
            RGB[145*3 +1] = 190;
            RGB[145*3 +2] = 135;
            RGB[146*3] = 102;
            RGB[146*3 +1] = 190;
            RGB[146*3 +2] = 133;
            RGB[147*3] = 105;
            RGB[147*3 +1] = 190;
            RGB[147*3 +2] = 132;
            RGB[148*3] = 108;
            RGB[148*3 +1] = 191;
            RGB[148*3 +2] = 131;
            RGB[149*3] = 111;
            RGB[149*3 +1] = 191;
            RGB[149*3 +2] = 129;
            RGB[150*3] = 113;
            RGB[150*3 +1] = 191;
            RGB[150*3 +2] = 128;
            RGB[151*3] = 116;
            RGB[151*3 +1] = 191;
            RGB[151*3 +2] = 127;
            RGB[152*3] = 119;
            RGB[152*3 +1] = 191;
            RGB[152*3 +2] = 126;
            RGB[153*3] = 122;
            RGB[153*3 +1] = 191;
            RGB[153*3 +2] = 124;
            RGB[154*3] = 125;
            RGB[154*3 +1] = 191;
            RGB[154*3 +2] = 123;
            RGB[155*3] = 128;
            RGB[155*3 +1] = 191;
            RGB[155*3 +2] = 122;
            RGB[156*3] = 130;
            RGB[156*3 +1] = 191;
            RGB[156*3 +2] = 121;
            RGB[157*3] = 133;
            RGB[157*3 +1] = 191;
            RGB[157*3 +2] = 120;
            RGB[158*3] = 136;
            RGB[158*3 +1] = 191;
            RGB[158*3 +2] = 119;
            RGB[159*3] = 138;
            RGB[159*3 +1] = 191;
            RGB[159*3 +2] = 118;
            RGB[160*3] = 141;
            RGB[160*3 +1] = 191;
            RGB[160*3 +2] = 117;
            RGB[161*3] = 143;
            RGB[161*3 +1] = 191;
            RGB[161*3 +2] = 116;
            RGB[162*3] = 146;
            RGB[162*3 +1] = 191;
            RGB[162*3 +2] = 114;
            RGB[163*3] = 148;
            RGB[163*3 +1] = 191;
            RGB[163*3 +2] = 114;
            RGB[164*3] = 151;
            RGB[164*3 +1] = 191;
            RGB[164*3 +2] = 113;
            RGB[165*3] = 153;
            RGB[165*3 +1] = 191;
            RGB[165*3 +2] = 112;
            RGB[166*3] = 156;
            RGB[166*3 +1] = 191;
            RGB[166*3 +2] = 111;
            RGB[167*3] = 158;
            RGB[167*3 +1] = 190;
            RGB[167*3 +2] = 110;
            RGB[168*3] = 160;
            RGB[168*3 +1] = 190;
            RGB[168*3 +2] = 109;
            RGB[169*3] = 163;
            RGB[169*3 +1] = 190;
            RGB[169*3 +2] = 108;
            RGB[170*3] = 165;
            RGB[170*3 +1] = 190;
            RGB[170*3 +2] = 107;
            RGB[171*3] = 167;
            RGB[171*3 +1] = 190;
            RGB[171*3 +2] = 106;
            RGB[172*3] = 170;
            RGB[172*3 +1] = 190;
            RGB[172*3 +2] = 105;
            RGB[173*3] = 172;
            RGB[173*3 +1] = 190;
            RGB[173*3 +2] = 104;
            RGB[174*3] = 174;
            RGB[174*3 +1] = 190;
            RGB[174*3 +2] = 103;
            RGB[175*3] = 176;
            RGB[175*3 +1] = 189;
            RGB[175*3 +2] = 102;
            RGB[176*3] = 179;
            RGB[176*3 +1] = 189;
            RGB[176*3 +2] = 101;
            RGB[177*3] = 181;
            RGB[177*3 +1] = 189;
            RGB[177*3 +2] = 101;
            RGB[178*3] = 183;
            RGB[178*3 +1] = 189;
            RGB[178*3 +2] = 100;
            RGB[179*3] = 185;
            RGB[179*3 +1] = 189;
            RGB[179*3 +2] = 99;
            RGB[180*3] = 187;
            RGB[180*3 +1] = 189;
            RGB[180*3 +2] = 98;
            RGB[181*3] = 189;
            RGB[181*3 +1] = 188;
            RGB[181*3 +2] = 97;
            RGB[182*3] = 192;
            RGB[182*3 +1] = 188;
            RGB[182*3 +2] = 96;
            RGB[183*3] = 194;
            RGB[183*3 +1] = 188;
            RGB[183*3 +2] = 95;
            RGB[184*3] = 196;
            RGB[184*3 +1] = 188;
            RGB[184*3 +2] = 95;
            RGB[185*3] = 198;
            RGB[185*3 +1] = 188;
            RGB[185*3 +2] = 94;
            RGB[186*3] = 200;
            RGB[186*3 +1] = 188;
            RGB[186*3 +2] = 93;
            RGB[187*3] = 202;
            RGB[187*3 +1] = 187;
            RGB[187*3 +2] = 92;
            RGB[188*3] = 204;
            RGB[188*3 +1] = 187;
            RGB[188*3 +2] = 91;
            RGB[189*3] = 206;
            RGB[189*3 +1] = 187;
            RGB[189*3 +2] = 90;
            RGB[190*3] = 208;
            RGB[190*3 +1] = 187;
            RGB[190*3 +2] = 89;
            RGB[191*3] = 210;
            RGB[191*3 +1] = 187;
            RGB[191*3 +2] = 89;
            RGB[192*3] = 212;
            RGB[192*3 +1] = 187;
            RGB[192*3 +2] = 88;
            RGB[193*3] = 214;
            RGB[193*3 +1] = 186;
            RGB[193*3 +2] = 87;
            RGB[194*3] = 216;
            RGB[194*3 +1] = 186;
            RGB[194*3 +2] = 86;
            RGB[195*3] = 218;
            RGB[195*3 +1] = 186;
            RGB[195*3 +2] = 85;
            RGB[196*3] = 220;
            RGB[196*3 +1] = 186;
            RGB[196*3 +2] = 84;
            RGB[197*3] = 222;
            RGB[197*3 +1] = 186;
            RGB[197*3 +2] = 83;
            RGB[198*3] = 224;
            RGB[198*3 +1] = 186;
            RGB[198*3 +2] = 82;
            RGB[199*3] = 226;
            RGB[199*3 +1] = 185;
            RGB[199*3 +2] = 81;
            RGB[200*3] = 228;
            RGB[200*3 +1] = 185;
            RGB[200*3 +2] = 80;
            RGB[201*3] = 230;
            RGB[201*3 +1] = 185;
            RGB[201*3 +2] = 79;
            RGB[202*3] = 232;
            RGB[202*3 +1] = 185;
            RGB[202*3 +2] = 78;
            RGB[203*3] = 234;
            RGB[203*3 +1] = 185;
            RGB[203*3 +2] = 77;
            RGB[204*3] = 236;
            RGB[204*3 +1] = 185;
            RGB[204*3 +2] = 76;
            RGB[205*3] = 238;
            RGB[205*3 +1] = 185;
            RGB[205*3 +2] = 75;
            RGB[206*3] = 240;
            RGB[206*3 +1] = 185;
            RGB[206*3 +2] = 74;
            RGB[207*3] = 242;
            RGB[207*3 +1] = 185;
            RGB[207*3 +2] = 73;
            RGB[208*3] = 244;
            RGB[208*3 +1] = 185;
            RGB[208*3 +2] = 72;
            RGB[209*3] = 246;
            RGB[209*3 +1] = 186;
            RGB[209*3 +2] = 70;
            RGB[210*3] = 248;
            RGB[210*3 +1] = 186;
            RGB[210*3 +2] = 69;
            RGB[211*3] = 249;
            RGB[211*3 +1] = 187;
            RGB[211*3 +2] = 67;
            RGB[212*3] = 251;
            RGB[212*3 +1] = 188;
            RGB[212*3 +2] = 66;
            RGB[213*3] = 252;
            RGB[213*3 +1] = 188;
            RGB[213*3 +2] = 64;
            RGB[214*3] = 253;
            RGB[214*3 +1] = 189;
            RGB[214*3 +2] = 62;
            RGB[215*3] = 254;
            RGB[215*3 +1] = 191;
            RGB[215*3 +2] = 61;
            RGB[216*3] = 254;
            RGB[216*3 +1] = 192;
            RGB[216*3 +2] = 59;
            RGB[217*3] = 255;
            RGB[217*3 +1] = 193;
            RGB[217*3 +2] = 57;
            RGB[218*3] = 255;
            RGB[218*3 +1] = 194;
            RGB[218*3 +2] = 56;
            RGB[219*3] = 255;
            RGB[219*3 +1] = 196;
            RGB[219*3 +2] = 55;
            RGB[220*3] = 255;
            RGB[220*3 +1] = 197;
            RGB[220*3 +2] = 53;
            RGB[221*3] = 254;
            RGB[221*3 +1] = 198;
            RGB[221*3 +2] = 52;
            RGB[222*3] = 254;
            RGB[222*3 +1] = 200;
            RGB[222*3 +2] = 51;
            RGB[223*3] = 254;
            RGB[223*3 +1] = 201;
            RGB[223*3 +2] = 50;
            RGB[224*3] = 253;
            RGB[224*3 +1] = 202;
            RGB[224*3 +2] = 49;
            RGB[225*3] = 253;
            RGB[225*3 +1] = 204;
            RGB[225*3 +2] = 48;
            RGB[226*3] = 252;
            RGB[226*3 +1] = 205;
            RGB[226*3 +2] = 46;
            RGB[227*3] = 252;
            RGB[227*3 +1] = 206;
            RGB[227*3 +2] = 45;
            RGB[228*3] = 251;
            RGB[228*3 +1] = 207;
            RGB[228*3 +2] = 44;
            RGB[229*3] = 251;
            RGB[229*3 +1] = 209;
            RGB[229*3 +2] = 43;
            RGB[230*3] = 250;
            RGB[230*3 +1] = 210;
            RGB[230*3 +2] = 42;
            RGB[231*3] = 249;
            RGB[231*3 +1] = 211;
            RGB[231*3 +2] = 41;
            RGB[232*3] = 249;
            RGB[232*3 +1] = 213;
            RGB[232*3 +2] = 40;
            RGB[233*3] = 248;
            RGB[233*3 +1] = 214;
            RGB[233*3 +2] = 39;
            RGB[234*3] = 248;
            RGB[234*3 +1] = 215;
            RGB[234*3 +2] = 38;
            RGB[235*3] = 247;
            RGB[235*3 +1] = 217;
            RGB[235*3 +2] = 37;
            RGB[236*3] = 247;
            RGB[236*3 +1] = 218;
            RGB[236*3 +2] = 36;
            RGB[237*3] = 246;
            RGB[237*3 +1] = 219;
            RGB[237*3 +2] = 35;
            RGB[238*3] = 246;
            RGB[238*3 +1] = 221;
            RGB[238*3 +2] = 34;
            RGB[239*3] = 245;
            RGB[239*3 +1] = 222;
            RGB[239*3 +2] = 33;
            RGB[240*3] = 245;
            RGB[240*3 +1] = 224;
            RGB[240*3 +2] = 32;
            RGB[241*3] = 245;
            RGB[241*3 +1] = 225;
            RGB[241*3 +2] = 31;
            RGB[242*3] = 245;
            RGB[242*3 +1] = 227;
            RGB[242*3 +2] = 30;
            RGB[243*3] = 244;
            RGB[243*3 +1] = 228;
            RGB[243*3 +2] = 29;
            RGB[244*3] = 244;
            RGB[244*3 +1] = 230;
            RGB[244*3 +2] = 28;
            RGB[245*3] = 244;
            RGB[245*3 +1] = 232;
            RGB[245*3 +2] = 26;
            RGB[246*3] = 245;
            RGB[246*3 +1] = 233;
            RGB[246*3 +2] = 25;
            RGB[247*3] = 245;
            RGB[247*3 +1] = 235;
            RGB[247*3 +2] = 24;
            RGB[248*3] = 245;
            RGB[248*3 +1] = 237;
            RGB[248*3 +2] = 23;
            RGB[249*3] = 245;
            RGB[249*3 +1] = 239;
            RGB[249*3 +2] = 22;
            RGB[250*3] = 246;
            RGB[250*3 +1] = 241;
            RGB[250*3 +2] = 20;
            RGB[251*3] = 246;
            RGB[251*3 +1] = 243;
            RGB[251*3 +2] = 19;
            RGB[252*3] = 247;
            RGB[252*3 +1] = 245;
            RGB[252*3 +2] = 18;
            RGB[253*3] = 248;
            RGB[253*3 +1] = 247;
            RGB[253*3 +2] = 17;
            RGB[254*3] = 248;
            RGB[254*3 +1] = 249;
            RGB[254*3 +2] = 15;
            RGB[255*3] = 249;
            RGB[255*3 +1] = 251;
            RGB[255*3 +2] = 14;
            break;
        case 3:
            RGB[0*3] = 3;
            RGB[0*3 +1] = 3;
            RGB[0*3 +2] = 3;
            RGB[1*3] = 9;
            RGB[1*3 +1] = 2;
            RGB[1*3 +2] = 7;
            RGB[2*3] = 12;
            RGB[2*3 +1] = 3;
            RGB[2*3 +2] = 10;
            RGB[3*3] = 15;
            RGB[3*3 +1] = 4;
            RGB[3*3 +2] = 13;
            RGB[4*3] = 17;
            RGB[4*3 +1] = 5;
            RGB[4*3 +2] = 15;
            RGB[5*3] = 20;
            RGB[5*3 +1] = 5;
            RGB[5*3 +2] = 17;
            RGB[6*3] = 21;
            RGB[6*3 +1] = 6;
            RGB[6*3 +2] = 19;
            RGB[7*3] = 23;
            RGB[7*3 +1] = 7;
            RGB[7*3 +2] = 20;
            RGB[8*3] = 25;
            RGB[8*3 +1] = 7;
            RGB[8*3 +2] = 21;
            RGB[9*3] = 26;
            RGB[9*3 +1] = 8;
            RGB[9*3 +2] = 23;
            RGB[10*3] = 27;
            RGB[10*3 +1] = 9;
            RGB[10*3 +2] = 24;
            RGB[11*3] = 29;
            RGB[11*3 +1] = 10;
            RGB[11*3 +2] = 25;
            RGB[12*3] = 30;
            RGB[12*3 +1] = 11;
            RGB[12*3 +2] = 26;
            RGB[13*3] = 31;
            RGB[13*3 +1] = 12;
            RGB[13*3 +2] = 27;
            RGB[14*3] = 32;
            RGB[14*3 +1] = 12;
            RGB[14*3 +2] = 28;
            RGB[15*3] = 34;
            RGB[15*3 +1] = 13;
            RGB[15*3 +2] = 30;
            RGB[16*3] = 35;
            RGB[16*3 +1] = 13;
            RGB[16*3 +2] = 31;
            RGB[17*3] = 36;
            RGB[17*3 +1] = 14;
            RGB[17*3 +2] = 31;
            RGB[18*3] = 38;
            RGB[18*3 +1] = 14;
            RGB[18*3 +2] = 33;
            RGB[19*3] = 38;
            RGB[19*3 +1] = 15;
            RGB[19*3 +2] = 33;
            RGB[20*3] = 39;
            RGB[20*3 +1] = 16;
            RGB[20*3 +2] = 34;
            RGB[21*3] = 41;
            RGB[21*3 +1] = 16;
            RGB[21*3 +2] = 35;
            RGB[22*3] = 41;
            RGB[22*3 +1] = 17;
            RGB[22*3 +2] = 37;
            RGB[23*3] = 43;
            RGB[23*3 +1] = 18;
            RGB[23*3 +2] = 38;
            RGB[24*3] = 44;
            RGB[24*3 +1] = 18;
            RGB[24*3 +2] = 38;
            RGB[25*3] = 45;
            RGB[25*3 +1] = 19;
            RGB[25*3 +2] = 40;
            RGB[26*3] = 45;
            RGB[26*3 +1] = 19;
            RGB[26*3 +2] = 44;
            RGB[27*3] = 45;
            RGB[27*3 +1] = 20;
            RGB[27*3 +2] = 49;
            RGB[28*3] = 45;
            RGB[28*3 +1] = 20;
            RGB[28*3 +2] = 54;
            RGB[29*3] = 45;
            RGB[29*3 +1] = 21;
            RGB[29*3 +2] = 59;
            RGB[30*3] = 44;
            RGB[30*3 +1] = 21;
            RGB[30*3 +2] = 63;
            RGB[31*3] = 44;
            RGB[31*3 +1] = 22;
            RGB[31*3 +2] = 67;
            RGB[32*3] = 44;
            RGB[32*3 +1] = 22;
            RGB[32*3 +2] = 70;
            RGB[33*3] = 44;
            RGB[33*3 +1] = 23;
            RGB[33*3 +2] = 74;
            RGB[34*3] = 44;
            RGB[34*3 +1] = 23;
            RGB[34*3 +2] = 76;
            RGB[35*3] = 44;
            RGB[35*3 +1] = 24;
            RGB[35*3 +2] = 80;
            RGB[36*3] = 44;
            RGB[36*3 +1] = 25;
            RGB[36*3 +2] = 83;
            RGB[37*3] = 44;
            RGB[37*3 +1] = 25;
            RGB[37*3 +2] = 86;
            RGB[38*3] = 43;
            RGB[38*3 +1] = 26;
            RGB[38*3 +2] = 89;
            RGB[39*3] = 44;
            RGB[39*3 +1] = 26;
            RGB[39*3 +2] = 92;
            RGB[40*3] = 44;
            RGB[40*3 +1] = 27;
            RGB[40*3 +2] = 95;
            RGB[41*3] = 43;
            RGB[41*3 +1] = 28;
            RGB[41*3 +2] = 98;
            RGB[42*3] = 43;
            RGB[42*3 +1] = 29;
            RGB[42*3 +2] = 100;
            RGB[43*3] = 43;
            RGB[43*3 +1] = 29;
            RGB[43*3 +2] = 103;
            RGB[44*3] = 43;
            RGB[44*3 +1] = 30;
            RGB[44*3 +2] = 106;
            RGB[45*3] = 43;
            RGB[45*3 +1] = 30;
            RGB[45*3 +2] = 108;
            RGB[46*3] = 43;
            RGB[46*3 +1] = 31;
            RGB[46*3 +2] = 110;
            RGB[47*3] = 43;
            RGB[47*3 +1] = 32;
            RGB[47*3 +2] = 113;
            RGB[48*3] = 43;
            RGB[48*3 +1] = 33;
            RGB[48*3 +2] = 115;
            RGB[49*3] = 43;
            RGB[49*3 +1] = 33;
            RGB[49*3 +2] = 118;
            RGB[50*3] = 42;
            RGB[50*3 +1] = 34;
            RGB[50*3 +2] = 120;
            RGB[51*3] = 39;
            RGB[51*3 +1] = 37;
            RGB[51*3 +2] = 120;
            RGB[52*3] = 37;
            RGB[52*3 +1] = 40;
            RGB[52*3 +2] = 119;
            RGB[53*3] = 35;
            RGB[53*3 +1] = 41;
            RGB[53*3 +2] = 119;
            RGB[54*3] = 34;
            RGB[54*3 +1] = 43;
            RGB[54*3 +2] = 119;
            RGB[55*3] = 32;
            RGB[55*3 +1] = 45;
            RGB[55*3 +2] = 120;
            RGB[56*3] = 31;
            RGB[56*3 +1] = 46;
            RGB[56*3 +2] = 120;
            RGB[57*3] = 30;
            RGB[57*3 +1] = 48;
            RGB[57*3 +2] = 121;
            RGB[58*3] = 29;
            RGB[58*3 +1] = 49;
            RGB[58*3 +2] = 121;
            RGB[59*3] = 29;
            RGB[59*3 +1] = 50;
            RGB[59*3 +2] = 122;
            RGB[60*3] = 28;
            RGB[60*3 +1] = 52;
            RGB[60*3 +2] = 122;
            RGB[61*3] = 28;
            RGB[61*3 +1] = 53;
            RGB[61*3 +2] = 124;
            RGB[62*3] = 27;
            RGB[62*3 +1] = 54;
            RGB[62*3 +2] = 125;
            RGB[63*3] = 27;
            RGB[63*3 +1] = 55;
            RGB[63*3 +2] = 126;
            RGB[64*3] = 26;
            RGB[64*3 +1] = 56;
            RGB[64*3 +2] = 127;
            RGB[65*3] = 27;
            RGB[65*3 +1] = 57;
            RGB[65*3 +2] = 128;
            RGB[66*3] = 27;
            RGB[66*3 +1] = 58;
            RGB[66*3 +2] = 128;
            RGB[67*3] = 26;
            RGB[67*3 +1] = 60;
            RGB[67*3 +2] = 129;
            RGB[68*3] = 26;
            RGB[68*3 +1] = 60;
            RGB[68*3 +2] = 131;
            RGB[69*3] = 26;
            RGB[69*3 +1] = 62;
            RGB[69*3 +2] = 132;
            RGB[70*3] = 26;
            RGB[70*3 +1] = 63;
            RGB[70*3 +2] = 133;
            RGB[71*3] = 27;
            RGB[71*3 +1] = 64;
            RGB[71*3 +2] = 134;
            RGB[72*3] = 26;
            RGB[72*3 +1] = 65;
            RGB[72*3 +2] = 134;
            RGB[73*3] = 26;
            RGB[73*3 +1] = 66;
            RGB[73*3 +2] = 136;
            RGB[74*3] = 26;
            RGB[74*3 +1] = 67;
            RGB[74*3 +2] = 137;
            RGB[75*3] = 26;
            RGB[75*3 +1] = 68;
            RGB[75*3 +2] = 138;
            RGB[76*3] = 24;
            RGB[76*3 +1] = 70;
            RGB[76*3 +2] = 136;
            RGB[77*3] = 22;
            RGB[77*3 +1] = 72;
            RGB[77*3 +2] = 133;
            RGB[78*3] = 20;
            RGB[78*3 +1] = 73;
            RGB[78*3 +2] = 133;
            RGB[79*3] = 19;
            RGB[79*3 +1] = 75;
            RGB[79*3 +2] = 131;
            RGB[80*3] = 17;
            RGB[80*3 +1] = 77;
            RGB[80*3 +2] = 129;
            RGB[81*3] = 15;
            RGB[81*3 +1] = 79;
            RGB[81*3 +2] = 127;
            RGB[82*3] = 14;
            RGB[82*3 +1] = 80;
            RGB[82*3 +2] = 126;
            RGB[83*3] = 13;
            RGB[83*3 +1] = 82;
            RGB[83*3 +2] = 124;
            RGB[84*3] = 9;
            RGB[84*3 +1] = 83;
            RGB[84*3 +2] = 122;
            RGB[85*3] = 9;
            RGB[85*3 +1] = 85;
            RGB[85*3 +2] = 122;
            RGB[86*3] = 8;
            RGB[86*3 +1] = 86;
            RGB[86*3 +2] = 120;
            RGB[87*3] = 7;
            RGB[87*3 +1] = 88;
            RGB[87*3 +2] = 118;
            RGB[88*3] = 6;
            RGB[88*3 +1] = 89;
            RGB[88*3 +2] = 118;
            RGB[89*3] = 5;
            RGB[89*3 +1] = 90;
            RGB[89*3 +2] = 116;
            RGB[90*3] = 3;
            RGB[90*3 +1] = 92;
            RGB[90*3 +2] = 115;
            RGB[91*3] = 3;
            RGB[91*3 +1] = 93;
            RGB[91*3 +2] = 114;
            RGB[92*3] = 3;
            RGB[92*3 +1] = 94;
            RGB[92*3 +2] = 113;
            RGB[93*3] = 2;
            RGB[93*3 +1] = 96;
            RGB[93*3 +2] = 112;
            RGB[94*3] = 0;
            RGB[94*3 +1] = 97;
            RGB[94*3 +2] = 110;
            RGB[95*3] = 0;
            RGB[95*3 +1] = 98;
            RGB[95*3 +2] = 110;
            RGB[96*3] = 0;
            RGB[96*3 +1] = 100;
            RGB[96*3 +2] = 108;
            RGB[97*3] = 0;
            RGB[97*3 +1] = 101;
            RGB[97*3 +2] = 107;
            RGB[98*3] = 0;
            RGB[98*3 +1] = 102;
            RGB[98*3 +2] = 107;
            RGB[99*3] = 0;
            RGB[99*3 +1] = 103;
            RGB[99*3 +2] = 106;
            RGB[100*3] = 0;
            RGB[100*3 +1] = 105;
            RGB[100*3 +2] = 104;
            RGB[101*3] = 0;
            RGB[101*3 +1] = 106;
            RGB[101*3 +2] = 103;
            RGB[102*3] = 0;
            RGB[102*3 +1] = 108;
            RGB[102*3 +2] = 101;
            RGB[103*3] = 0;
            RGB[103*3 +1] = 109;
            RGB[103*3 +2] = 99;
            RGB[104*3] = 0;
            RGB[104*3 +1] = 110;
            RGB[104*3 +2] = 99;
            RGB[105*3] = 0;
            RGB[105*3 +1] = 112;
            RGB[105*3 +2] = 98;
            RGB[106*3] = 0;
            RGB[106*3 +1] = 113;
            RGB[106*3 +2] = 96;
            RGB[107*3] = 0;
            RGB[107*3 +1] = 114;
            RGB[107*3 +2] = 94;
            RGB[108*3] = 0;
            RGB[108*3 +1] = 116;
            RGB[108*3 +2] = 93;
            RGB[109*3] = 0;
            RGB[109*3 +1] = 117;
            RGB[109*3 +2] = 91;
            RGB[110*3] = 0;
            RGB[110*3 +1] = 119;
            RGB[110*3 +2] = 90;
            RGB[111*3] = 0;
            RGB[111*3 +1] = 120;
            RGB[111*3 +2] = 89;
            RGB[112*3] = 0;
            RGB[112*3 +1] = 121;
            RGB[112*3 +2] = 88;
            RGB[113*3] = 0;
            RGB[113*3 +1] = 122;
            RGB[113*3 +2] = 85;
            RGB[114*3] = 0;
            RGB[114*3 +1] = 124;
            RGB[114*3 +2] = 84;
            RGB[115*3] = 0;
            RGB[115*3 +1] = 125;
            RGB[115*3 +2] = 82;
            RGB[116*3] = 0;
            RGB[116*3 +1] = 126;
            RGB[116*3 +2] = 81;
            RGB[117*3] = 0;
            RGB[117*3 +1] = 128;
            RGB[117*3 +2] = 79;
            RGB[118*3] = 0;
            RGB[118*3 +1] = 129;
            RGB[118*3 +2] = 78;
            RGB[119*3] = 0;
            RGB[119*3 +1] = 130;
            RGB[119*3 +2] = 76;
            RGB[120*3] = 0;
            RGB[120*3 +1] = 131;
            RGB[120*3 +2] = 74;
            RGB[121*3] = 0;
            RGB[121*3 +1] = 132;
            RGB[121*3 +2] = 72;
            RGB[122*3] = 0;
            RGB[122*3 +1] = 134;
            RGB[122*3 +2] = 71;
            RGB[123*3] = 0;
            RGB[123*3 +1] = 135;
            RGB[123*3 +2] = 69;
            RGB[124*3] = 0;
            RGB[124*3 +1] = 136;
            RGB[124*3 +2] = 67;
            RGB[125*3] = 0;
            RGB[125*3 +1] = 137;
            RGB[125*3 +2] = 65;
            RGB[126*3] = 0;
            RGB[126*3 +1] = 139;
            RGB[126*3 +2] = 64;
            RGB[127*3] = 0;
            RGB[127*3 +1] = 140;
            RGB[127*3 +2] = 63;
            RGB[128*3] = 0;
            RGB[128*3 +1] = 141;
            RGB[128*3 +2] = 62;
            RGB[129*3] = 0;
            RGB[129*3 +1] = 143;
            RGB[129*3 +2] = 60;
            RGB[130*3] = 0;
            RGB[130*3 +1] = 144;
            RGB[130*3 +2] = 59;
            RGB[131*3] = 0;
            RGB[131*3 +1] = 146;
            RGB[131*3 +2] = 59;
            RGB[132*3] = 0;
            RGB[132*3 +1] = 147;
            RGB[132*3 +2] = 58;
            RGB[133*3] = 0;
            RGB[133*3 +1] = 148;
            RGB[133*3 +2] = 57;
            RGB[134*3] = 0;
            RGB[134*3 +1] = 150;
            RGB[134*3 +2] = 55;
            RGB[135*3] = 0;
            RGB[135*3 +1] = 151;
            RGB[135*3 +2] = 54;
            RGB[136*3] = 0;
            RGB[136*3 +1] = 152;
            RGB[136*3 +2] = 52;
            RGB[137*3] = 0;
            RGB[137*3 +1] = 154;
            RGB[137*3 +2] = 52;
            RGB[138*3] = 0;
            RGB[138*3 +1] = 155;
            RGB[138*3 +2] = 50;
            RGB[139*3] = 0;
            RGB[139*3 +1] = 156;
            RGB[139*3 +2] = 49;
            RGB[140*3] = 0;
            RGB[140*3 +1] = 157;
            RGB[140*3 +2] = 48;
            RGB[141*3] = 0;
            RGB[141*3 +1] = 159;
            RGB[141*3 +2] = 46;
            RGB[142*3] = 0;
            RGB[142*3 +1] = 160;
            RGB[142*3 +2] = 44;
            RGB[143*3] = 0;
            RGB[143*3 +1] = 161;
            RGB[143*3 +2] = 43;
            RGB[144*3] = 0;
            RGB[144*3 +1] = 162;
            RGB[144*3 +2] = 41;
            RGB[145*3] = 0;
            RGB[145*3 +1] = 164;
            RGB[145*3 +2] = 40;
            RGB[146*3] = 0;
            RGB[146*3 +1] = 165;
            RGB[146*3 +2] = 38;
            RGB[147*3] = 0;
            RGB[147*3 +1] = 166;
            RGB[147*3 +2] = 35;
            RGB[148*3] = 0;
            RGB[148*3 +1] = 168;
            RGB[148*3 +2] = 33;
            RGB[149*3] = 0;
            RGB[149*3 +1] = 169;
            RGB[149*3 +2] = 30;
            RGB[150*3] = 0;
            RGB[150*3 +1] = 169;
            RGB[150*3 +2] = 22;
            RGB[151*3] = 0;
            RGB[151*3 +1] = 169;
            RGB[151*3 +2] = 7;
            RGB[152*3] = 8;
            RGB[152*3 +1] = 169;
            RGB[152*3 +2] = 0;
            RGB[153*3] = 19;
            RGB[153*3 +1] = 171;
            RGB[153*3 +2] = 0;
            RGB[154*3] = 28;
            RGB[154*3 +1] = 172;
            RGB[154*3 +2] = 0;
            RGB[155*3] = 35;
            RGB[155*3 +1] = 172;
            RGB[155*3 +2] = 0;
            RGB[156*3] = 39;
            RGB[156*3 +1] = 173;
            RGB[156*3 +2] = 0;
            RGB[157*3] = 44;
            RGB[157*3 +1] = 174;
            RGB[157*3 +2] = 0;
            RGB[158*3] = 47;
            RGB[158*3 +1] = 175;
            RGB[158*3 +2] = 0;
            RGB[159*3] = 51;
            RGB[159*3 +1] = 176;
            RGB[159*3 +2] = 0;
            RGB[160*3] = 54;
            RGB[160*3 +1] = 177;
            RGB[160*3 +2] = 0;
            RGB[161*3] = 58;
            RGB[161*3 +1] = 178;
            RGB[161*3 +2] = 0;
            RGB[162*3] = 60;
            RGB[162*3 +1] = 179;
            RGB[162*3 +2] = 0;
            RGB[163*3] = 63;
            RGB[163*3 +1] = 180;
            RGB[163*3 +2] = 0;
            RGB[164*3] = 65;
            RGB[164*3 +1] = 181;
            RGB[164*3 +2] = 0;
            RGB[165*3] = 68;
            RGB[165*3 +1] = 182;
            RGB[165*3 +2] = 0;
            RGB[166*3] = 70;
            RGB[166*3 +1] = 183;
            RGB[166*3 +2] = 0;
            RGB[167*3] = 73;
            RGB[167*3 +1] = 184;
            RGB[167*3 +2] = 0;
            RGB[168*3] = 75;
            RGB[168*3 +1] = 185;
            RGB[168*3 +2] = 0;
            RGB[169*3] = 77;
            RGB[169*3 +1] = 186;
            RGB[169*3 +2] = 0;
            RGB[170*3] = 78;
            RGB[170*3 +1] = 187;
            RGB[170*3 +2] = 0;
            RGB[171*3] = 81;
            RGB[171*3 +1] = 188;
            RGB[171*3 +2] = 0;
            RGB[172*3] = 83;
            RGB[172*3 +1] = 189;
            RGB[172*3 +2] = 0;
            RGB[173*3] = 85;
            RGB[173*3 +1] = 191;
            RGB[173*3 +2] = 0;
            RGB[174*3] = 86;
            RGB[174*3 +1] = 192;
            RGB[174*3 +2] = 0;
            RGB[175*3] = 91;
            RGB[175*3 +1] = 192;
            RGB[175*3 +2] = 0;
            RGB[176*3] = 99;
            RGB[176*3 +1] = 193;
            RGB[176*3 +2] = 0;
            RGB[177*3] = 106;
            RGB[177*3 +1] = 193;
            RGB[177*3 +2] = 0;
            RGB[178*3] = 113;
            RGB[178*3 +1] = 193;
            RGB[178*3 +2] = 0;
            RGB[179*3] = 119;
            RGB[179*3 +1] = 193;
            RGB[179*3 +2] = 0;
            RGB[180*3] = 124;
            RGB[180*3 +1] = 194;
            RGB[180*3 +2] = 0;
            RGB[181*3] = 129;
            RGB[181*3 +1] = 194;
            RGB[181*3 +2] = 0;
            RGB[182*3] = 135;
            RGB[182*3 +1] = 195;
            RGB[182*3 +2] = 0;
            RGB[183*3] = 140;
            RGB[183*3 +1] = 195;
            RGB[183*3 +2] = 0;
            RGB[184*3] = 145;
            RGB[184*3 +1] = 195;
            RGB[184*3 +2] = 0;
            RGB[185*3] = 148;
            RGB[185*3 +1] = 196;
            RGB[185*3 +2] = 0;
            RGB[186*3] = 153;
            RGB[186*3 +1] = 196;
            RGB[186*3 +2] = 0;
            RGB[187*3] = 157;
            RGB[187*3 +1] = 196;
            RGB[187*3 +2] = 0;
            RGB[188*3] = 161;
            RGB[188*3 +1] = 197;
            RGB[188*3 +2] = 0;
            RGB[189*3] = 165;
            RGB[189*3 +1] = 197;
            RGB[189*3 +2] = 0;
            RGB[190*3] = 169;
            RGB[190*3 +1] = 198;
            RGB[190*3 +2] = 0;
            RGB[191*3] = 172;
            RGB[191*3 +1] = 198;
            RGB[191*3 +2] = 0;
            RGB[192*3] = 176;
            RGB[192*3 +1] = 199;
            RGB[192*3 +2] = 0;
            RGB[193*3] = 180;
            RGB[193*3 +1] = 199;
            RGB[193*3 +2] = 0;
            RGB[194*3] = 184;
            RGB[194*3 +1] = 199;
            RGB[194*3 +2] = 0;
            RGB[195*3] = 186;
            RGB[195*3 +1] = 200;
            RGB[195*3 +2] = 0;
            RGB[196*3] = 190;
            RGB[196*3 +1] = 200;
            RGB[196*3 +2] = 0;
            RGB[197*3] = 193;
            RGB[197*3 +1] = 201;
            RGB[197*3 +2] = 0;
            RGB[198*3] = 196;
            RGB[198*3 +1] = 201;
            RGB[198*3 +2] = 0;
            RGB[199*3] = 200;
            RGB[199*3 +1] = 202;
            RGB[199*3 +2] = 0;
            RGB[200*3] = 201;
            RGB[200*3 +1] = 201;
            RGB[200*3 +2] = 23;
            RGB[201*3] = 203;
            RGB[201*3 +1] = 201;
            RGB[201*3 +2] = 50;
            RGB[202*3] = 206;
            RGB[202*3 +1] = 202;
            RGB[202*3 +2] = 65;
            RGB[203*3] = 207;
            RGB[203*3 +1] = 202;
            RGB[203*3 +2] = 76;
            RGB[204*3] = 209;
            RGB[204*3 +1] = 203;
            RGB[204*3 +2] = 86;
            RGB[205*3] = 212;
            RGB[205*3 +1] = 203;
            RGB[205*3 +2] = 94;
            RGB[206*3] = 213;
            RGB[206*3 +1] = 204;
            RGB[206*3 +2] = 102;
            RGB[207*3] = 215;
            RGB[207*3 +1] = 204;
            RGB[207*3 +2] = 109;
            RGB[208*3] = 218;
            RGB[208*3 +1] = 205;
            RGB[208*3 +2] = 116;
            RGB[209*3] = 219;
            RGB[209*3 +1] = 206;
            RGB[209*3 +2] = 120;
            RGB[210*3] = 221;
            RGB[210*3 +1] = 206;
            RGB[210*3 +2] = 127;
            RGB[211*3] = 223;
            RGB[211*3 +1] = 207;
            RGB[211*3 +2] = 132;
            RGB[212*3] = 226;
            RGB[212*3 +1] = 207;
            RGB[212*3 +2] = 137;
            RGB[213*3] = 227;
            RGB[213*3 +1] = 208;
            RGB[213*3 +2] = 142;
            RGB[214*3] = 229;
            RGB[214*3 +1] = 208;
            RGB[214*3 +2] = 147;
            RGB[215*3] = 231;
            RGB[215*3 +1] = 209;
            RGB[215*3 +2] = 151;
            RGB[216*3] = 232;
            RGB[216*3 +1] = 210;
            RGB[216*3 +2] = 155;
            RGB[217*3] = 235;
            RGB[217*3 +1] = 210;
            RGB[217*3 +2] = 159;
            RGB[218*3] = 237;
            RGB[218*3 +1] = 211;
            RGB[218*3 +2] = 163;
            RGB[219*3] = 238;
            RGB[219*3 +1] = 212;
            RGB[219*3 +2] = 167;
            RGB[220*3] = 240;
            RGB[220*3 +1] = 212;
            RGB[220*3 +2] = 172;
            RGB[221*3] = 243;
            RGB[221*3 +1] = 213;
            RGB[221*3 +2] = 175;
            RGB[222*3] = 243;
            RGB[222*3 +1] = 214;
            RGB[222*3 +2] = 179;
            RGB[223*3] = 245;
            RGB[223*3 +1] = 214;
            RGB[223*3 +2] = 183;
            RGB[224*3] = 248;
            RGB[224*3 +1] = 215;
            RGB[224*3 +2] = 186;
            RGB[225*3] = 248;
            RGB[225*3 +1] = 216;
            RGB[225*3 +2] = 189;
            RGB[226*3] = 248;
            RGB[226*3 +1] = 217;
            RGB[226*3 +2] = 192;
            RGB[227*3] = 247;
            RGB[227*3 +1] = 219;
            RGB[227*3 +2] = 195;
            RGB[228*3] = 247;
            RGB[228*3 +1] = 220;
            RGB[228*3 +2] = 197;
            RGB[229*3] = 247;
            RGB[229*3 +1] = 221;
            RGB[229*3 +2] = 200;
            RGB[230*3] = 248;
            RGB[230*3 +1] = 223;
            RGB[230*3 +2] = 204;
            RGB[231*3] = 247;
            RGB[231*3 +1] = 224;
            RGB[231*3 +2] = 206;
            RGB[232*3] = 247;
            RGB[232*3 +1] = 226;
            RGB[232*3 +2] = 209;
            RGB[233*3] = 247;
            RGB[233*3 +1] = 227;
            RGB[233*3 +2] = 211;
            RGB[234*3] = 247;
            RGB[234*3 +1] = 228;
            RGB[234*3 +2] = 213;
            RGB[235*3] = 247;
            RGB[235*3 +1] = 230;
            RGB[235*3 +2] = 216;
            RGB[236*3] = 247;
            RGB[236*3 +1] = 231;
            RGB[236*3 +2] = 218;
            RGB[237*3] = 247;
            RGB[237*3 +1] = 232;
            RGB[237*3 +2] = 220;
            RGB[238*3] = 248;
            RGB[238*3 +1] = 234;
            RGB[238*3 +2] = 224;
            RGB[239*3] = 247;
            RGB[239*3 +1] = 235;
            RGB[239*3 +2] = 225;
            RGB[240*3] = 247;
            RGB[240*3 +1] = 236;
            RGB[240*3 +2] = 228;
            RGB[241*3] = 247;
            RGB[241*3 +1] = 238;
            RGB[241*3 +2] = 231;
            RGB[242*3] = 247;
            RGB[242*3 +1] = 239;
            RGB[242*3 +2] = 232;
            RGB[243*3] = 248;
            RGB[243*3 +1] = 240;
            RGB[243*3 +2] = 235;
            RGB[244*3] = 248;
            RGB[244*3 +1] = 242;
            RGB[244*3 +2] = 237;
            RGB[245*3] = 247;
            RGB[245*3 +1] = 243;
            RGB[245*3 +2] = 239;
            RGB[246*3] = 247;
            RGB[246*3 +1] = 244;
            RGB[246*3 +2] = 241;
            RGB[247*3] = 248;
            RGB[247*3 +1] = 246;
            RGB[247*3 +2] = 244;
            RGB[248*3] = 248;
            RGB[248*3 +1] = 247;
            RGB[248*3 +2] = 246;
            RGB[249*3] = 248;
            RGB[249*3 +1] = 248;
            RGB[249*3 +2] = 248;
            RGB[250*3] = 249;
            RGB[250*3 +1] = 249;
            RGB[250*3 +2] = 249;
            RGB[251*3] = 250;
            RGB[251*3 +1] = 250;
            RGB[251*3 +2] = 250;
            RGB[252*3] = 252;
            RGB[252*3 +1] = 252;
            RGB[252*3 +2] = 252;
            RGB[253*3] = 253;
            RGB[253*3 +1] = 253;
            RGB[253*3 +2] = 253;
            RGB[254*3] = 254;
            RGB[254*3 +1] = 254;
            RGB[254*3 +2] = 254;
            RGB[255*3] = 255;
            RGB[255*3 +1] = 255;
            RGB[255*3 +2] = 255;
            break;
        case 4:
            RGB[0*3] = 255;
            RGB[0*3 +1] = 138;
            RGB[0*3 +2] = 255;
            RGB[1*3] = 255;
            RGB[1*3 +1] = 135;
            RGB[1*3 +2] = 255;
            RGB[2*3] = 255;
            RGB[2*3 +1] = 133;
            RGB[2*3 +2] = 255;
            RGB[3*3] = 255;
            RGB[3*3 +1] = 131;
            RGB[3*3 +2] = 255;
            RGB[4*3] = 255;
            RGB[4*3 +1] = 129;
            RGB[4*3 +2] = 255;
            RGB[5*3] = 255;
            RGB[5*3 +1] = 127;
            RGB[5*3 +2] = 255;
            RGB[6*3] = 255;
            RGB[6*3 +1] = 125;
            RGB[6*3 +2] = 255;
            RGB[7*3] = 255;
            RGB[7*3 +1] = 123;
            RGB[7*3 +2] = 255;
            RGB[8*3] = 255;
            RGB[8*3 +1] = 121;
            RGB[8*3 +2] = 255;
            RGB[9*3] = 255;
            RGB[9*3 +1] = 118;
            RGB[9*3 +2] = 255;
            RGB[10*3] = 255;
            RGB[10*3 +1] = 116;
            RGB[10*3 +2] = 255;
            RGB[11*3] = 255;
            RGB[11*3 +1] = 114;
            RGB[11*3 +2] = 255;
            RGB[12*3] = 255;
            RGB[12*3 +1] = 112;
            RGB[12*3 +2] = 255;
            RGB[13*3] = 255;
            RGB[13*3 +1] = 110;
            RGB[13*3 +2] = 255;
            RGB[14*3] = 253;
            RGB[14*3 +1] = 108;
            RGB[14*3 +2] = 255;
            RGB[15*3] = 251;
            RGB[15*3 +1] = 105;
            RGB[15*3 +2] = 255;
            RGB[16*3] = 249;
            RGB[16*3 +1] = 103;
            RGB[16*3 +2] = 255;
            RGB[17*3] = 245;
            RGB[17*3 +1] = 101;
            RGB[17*3 +2] = 255;
            RGB[18*3] = 239;
            RGB[18*3 +1] = 99;
            RGB[18*3 +2] = 255;
            RGB[19*3] = 232;
            RGB[19*3 +1] = 96;
            RGB[19*3 +2] = 255;
            RGB[20*3] = 225;
            RGB[20*3 +1] = 94;
            RGB[20*3 +2] = 255;
            RGB[21*3] = 219;
            RGB[21*3 +1] = 92;
            RGB[21*3 +2] = 255;
            RGB[22*3] = 212;
            RGB[22*3 +1] = 89;
            RGB[22*3 +2] = 255;
            RGB[23*3] = 206;
            RGB[23*3 +1] = 87;
            RGB[23*3 +2] = 255;
            RGB[24*3] = 200;
            RGB[24*3 +1] = 85;
            RGB[24*3 +2] = 255;
            RGB[25*3] = 194;
            RGB[25*3 +1] = 82;
            RGB[25*3 +2] = 255;
            RGB[26*3] = 188;
            RGB[26*3 +1] = 80;
            RGB[26*3 +2] = 255;
            RGB[27*3] = 181;
            RGB[27*3 +1] = 77;
            RGB[27*3 +2] = 255;
            RGB[28*3] = 175;
            RGB[28*3 +1] = 75;
            RGB[28*3 +2] = 255;
            RGB[29*3] = 169;
            RGB[29*3 +1] = 72;
            RGB[29*3 +2] = 255;
            RGB[30*3] = 163;
            RGB[30*3 +1] = 70;
            RGB[30*3 +2] = 255;
            RGB[31*3] = 157;
            RGB[31*3 +1] = 67;
            RGB[31*3 +2] = 255;
            RGB[32*3] = 150;
            RGB[32*3 +1] = 65;
            RGB[32*3 +2] = 255;
            RGB[33*3] = 144;
            RGB[33*3 +1] = 62;
            RGB[33*3 +2] = 255;
            RGB[34*3] = 138;
            RGB[34*3 +1] = 59;
            RGB[34*3 +2] = 255;
            RGB[35*3] = 132;
            RGB[35*3 +1] = 57;
            RGB[35*3 +2] = 255;
            RGB[36*3] = 125;
            RGB[36*3 +1] = 54;
            RGB[36*3 +2] = 255;
            RGB[37*3] = 119;
            RGB[37*3 +1] = 51;
            RGB[37*3 +2] = 255;
            RGB[38*3] = 112;
            RGB[38*3 +1] = 48;
            RGB[38*3 +2] = 255;
            RGB[39*3] = 106;
            RGB[39*3 +1] = 45;
            RGB[39*3 +2] = 255;
            RGB[40*3] = 99;
            RGB[40*3 +1] = 42;
            RGB[40*3 +2] = 255;
            RGB[41*3] = 93;
            RGB[41*3 +1] = 39;
            RGB[41*3 +2] = 255;
            RGB[42*3] = 86;
            RGB[42*3 +1] = 36;
            RGB[42*3 +2] = 255;
            RGB[43*3] = 79;
            RGB[43*3 +1] = 32;
            RGB[43*3 +2] = 255;
            RGB[44*3] = 72;
            RGB[44*3 +1] = 28;
            RGB[44*3 +2] = 255;
            RGB[45*3] = 63;
            RGB[45*3 +1] = 24;
            RGB[45*3 +2] = 255;
            RGB[46*3] = 53;
            RGB[46*3 +1] = 18;
            RGB[46*3 +2] = 255;
            RGB[47*3] = 43;
            RGB[47*3 +1] = 11;
            RGB[47*3 +2] = 255;
            RGB[48*3] = 33;
            RGB[48*3 +1] = 6;
            RGB[48*3 +2] = 255;
            RGB[49*3] = 24;
            RGB[49*3 +1] = 6;
            RGB[49*3 +2] = 255;
            RGB[50*3] = 15;
            RGB[50*3 +1] = 9;
            RGB[50*3 +2] = 254;
            RGB[51*3] = 7;
            RGB[51*3 +1] = 15;
            RGB[51*3 +2] = 253;
            RGB[52*3] = 1;
            RGB[52*3 +1] = 20;
            RGB[52*3 +2] = 252;
            RGB[53*3] = 1;
            RGB[53*3 +1] = 23;
            RGB[53*3 +2] = 252;
            RGB[54*3] = 7;
            RGB[54*3 +1] = 25;
            RGB[54*3 +2] = 253;
            RGB[55*3] = 17;
            RGB[55*3 +1] = 27;
            RGB[55*3 +2] = 254;
            RGB[56*3] = 27;
            RGB[56*3 +1] = 30;
            RGB[56*3 +2] = 255;
            RGB[57*3] = 33;
            RGB[57*3 +1] = 33;
            RGB[57*3 +2] = 255;
            RGB[58*3] = 39;
            RGB[58*3 +1] = 38;
            RGB[58*3 +2] = 255;
            RGB[59*3] = 43;
            RGB[59*3 +1] = 45;
            RGB[59*3 +2] = 255;
            RGB[60*3] = 47;
            RGB[60*3 +1] = 52;
            RGB[60*3 +2] = 255;
            RGB[61*3] = 50;
            RGB[61*3 +1] = 58;
            RGB[61*3 +2] = 255;
            RGB[62*3] = 52;
            RGB[62*3 +1] = 64;
            RGB[62*3 +2] = 255;
            RGB[63*3] = 54;
            RGB[63*3 +1] = 69;
            RGB[63*3 +2] = 255;
            RGB[64*3] = 55;
            RGB[64*3 +1] = 75;
            RGB[64*3 +2] = 255;
            RGB[65*3] = 56;
            RGB[65*3 +1] = 80;
            RGB[65*3 +2] = 255;
            RGB[66*3] = 56;
            RGB[66*3 +1] = 85;
            RGB[66*3 +2] = 255;
            RGB[67*3] = 57;
            RGB[67*3 +1] = 91;
            RGB[67*3 +2] = 255;
            RGB[68*3] = 57;
            RGB[68*3 +1] = 96;
            RGB[68*3 +2] = 255;
            RGB[69*3] = 57;
            RGB[69*3 +1] = 101;
            RGB[69*3 +2] = 255;
            RGB[70*3] = 57;
            RGB[70*3 +1] = 106;
            RGB[70*3 +2] = 255;
            RGB[71*3] = 57;
            RGB[71*3 +1] = 111;
            RGB[71*3 +2] = 255;
            RGB[72*3] = 56;
            RGB[72*3 +1] = 116;
            RGB[72*3 +2] = 255;
            RGB[73*3] = 55;
            RGB[73*3 +1] = 120;
            RGB[73*3 +2] = 255;
            RGB[74*3] = 55;
            RGB[74*3 +1] = 125;
            RGB[74*3 +2] = 255;
            RGB[75*3] = 54;
            RGB[75*3 +1] = 130;
            RGB[75*3 +2] = 255;
            RGB[76*3] = 53;
            RGB[76*3 +1] = 135;
            RGB[76*3 +2] = 255;
            RGB[77*3] = 52;
            RGB[77*3 +1] = 139;
            RGB[77*3 +2] = 255;
            RGB[78*3] = 51;
            RGB[78*3 +1] = 144;
            RGB[78*3 +2] = 255;
            RGB[79*3] = 50;
            RGB[79*3 +1] = 148;
            RGB[79*3 +2] = 255;
            RGB[80*3] = 48;
            RGB[80*3 +1] = 153;
            RGB[80*3 +2] = 255;
            RGB[81*3] = 47;
            RGB[81*3 +1] = 157;
            RGB[81*3 +2] = 255;
            RGB[82*3] = 46;
            RGB[82*3 +1] = 162;
            RGB[82*3 +2] = 255;
            RGB[83*3] = 44;
            RGB[83*3 +1] = 166;
            RGB[83*3 +2] = 255;
            RGB[84*3] = 43;
            RGB[84*3 +1] = 171;
            RGB[84*3 +2] = 255;
            RGB[85*3] = 41;
            RGB[85*3 +1] = 175;
            RGB[85*3 +2] = 255;
            RGB[86*3] = 40;
            RGB[86*3 +1] = 179;
            RGB[86*3 +2] = 255;
            RGB[87*3] = 38;
            RGB[87*3 +1] = 184;
            RGB[87*3 +2] = 255;
            RGB[88*3] = 36;
            RGB[88*3 +1] = 188;
            RGB[88*3 +2] = 255;
            RGB[89*3] = 34;
            RGB[89*3 +1] = 192;
            RGB[89*3 +2] = 255;
            RGB[90*3] = 32;
            RGB[90*3 +1] = 196;
            RGB[90*3 +2] = 255;
            RGB[91*3] = 30;
            RGB[91*3 +1] = 201;
            RGB[91*3 +2] = 255;
            RGB[92*3] = 28;
            RGB[92*3 +1] = 205;
            RGB[92*3 +2] = 255;
            RGB[93*3] = 25;
            RGB[93*3 +1] = 209;
            RGB[93*3 +2] = 255;
            RGB[94*3] = 22;
            RGB[94*3 +1] = 213;
            RGB[94*3 +2] = 255;
            RGB[95*3] = 19;
            RGB[95*3 +1] = 218;
            RGB[95*3 +2] = 255;
            RGB[96*3] = 15;
            RGB[96*3 +1] = 222;
            RGB[96*3 +2] = 255;
            RGB[97*3] = 11;
            RGB[97*3 +1] = 226;
            RGB[97*3 +2] = 255;
            RGB[98*3] = 8;
            RGB[98*3 +1] = 231;
            RGB[98*3 +2] = 255;
            RGB[99*3] = 4;
            RGB[99*3 +1] = 236;
            RGB[99*3 +2] = 255;
            RGB[100*3] = 1;
            RGB[100*3 +1] = 240;
            RGB[100*3 +2] = 254;
            RGB[101*3] = 0;
            RGB[101*3 +1] = 243;
            RGB[101*3 +2] = 254;
            RGB[102*3] = 0;
            RGB[102*3 +1] = 245;
            RGB[102*3 +2] = 252;
            RGB[103*3] = 0;
            RGB[103*3 +1] = 246;
            RGB[103*3 +2] = 249;
            RGB[104*3] = 0;
            RGB[104*3 +1] = 247;
            RGB[104*3 +2] = 245;
            RGB[105*3] = 0;
            RGB[105*3 +1] = 248;
            RGB[105*3 +2] = 239;
            RGB[106*3] = 0;
            RGB[106*3 +1] = 247;
            RGB[106*3 +2] = 234;
            RGB[107*3] = 0;
            RGB[107*3 +1] = 245;
            RGB[107*3 +2] = 227;
            RGB[108*3] = 0;
            RGB[108*3 +1] = 242;
            RGB[108*3 +2] = 220;
            RGB[109*3] = 0;
            RGB[109*3 +1] = 239;
            RGB[109*3 +2] = 212;
            RGB[110*3] = 0;
            RGB[110*3 +1] = 237;
            RGB[110*3 +2] = 205;
            RGB[111*3] = 0;
            RGB[111*3 +1] = 234;
            RGB[111*3 +2] = 199;
            RGB[112*3] = 0;
            RGB[112*3 +1] = 232;
            RGB[112*3 +2] = 192;
            RGB[113*3] = 0;
            RGB[113*3 +1] = 229;
            RGB[113*3 +2] = 185;
            RGB[114*3] = 0;
            RGB[114*3 +1] = 227;
            RGB[114*3 +2] = 179;
            RGB[115*3] = 0;
            RGB[115*3 +1] = 224;
            RGB[115*3 +2] = 172;
            RGB[116*3] = 0;
            RGB[116*3 +1] = 222;
            RGB[116*3 +2] = 166;
            RGB[117*3] = 0;
            RGB[117*3 +1] = 219;
            RGB[117*3 +2] = 160;
            RGB[118*3] = 0;
            RGB[118*3 +1] = 216;
            RGB[118*3 +2] = 154;
            RGB[119*3] = 0;
            RGB[119*3 +1] = 214;
            RGB[119*3 +2] = 148;
            RGB[120*3] = 0;
            RGB[120*3 +1] = 211;
            RGB[120*3 +2] = 142;
            RGB[121*3] = 0;
            RGB[121*3 +1] = 209;
            RGB[121*3 +2] = 136;
            RGB[122*3] = 0;
            RGB[122*3 +1] = 206;
            RGB[122*3 +2] = 130;
            RGB[123*3] = 0;
            RGB[123*3 +1] = 203;
            RGB[123*3 +2] = 125;
            RGB[124*3] = 0;
            RGB[124*3 +1] = 201;
            RGB[124*3 +2] = 120;
            RGB[125*3] = 0;
            RGB[125*3 +1] = 198;
            RGB[125*3 +2] = 115;
            RGB[126*3] = 0;
            RGB[126*3 +1] = 195;
            RGB[126*3 +2] = 110;
            RGB[127*3] = 0;
            RGB[127*3 +1] = 193;
            RGB[127*3 +2] = 105;
            RGB[128*3] = 0;
            RGB[128*3 +1] = 190;
            RGB[128*3 +2] = 101;
            RGB[129*3] = 0;
            RGB[129*3 +1] = 187;
            RGB[129*3 +2] = 97;
            RGB[130*3] = 0;
            RGB[130*3 +1] = 185;
            RGB[130*3 +2] = 93;
            RGB[131*3] = 0;
            RGB[131*3 +1] = 182;
            RGB[131*3 +2] = 90;
            RGB[132*3] = 0;
            RGB[132*3 +1] = 179;
            RGB[132*3 +2] = 87;
            RGB[133*3] = 0;
            RGB[133*3 +1] = 177;
            RGB[133*3 +2] = 84;
            RGB[134*3] = 0;
            RGB[134*3 +1] = 174;
            RGB[134*3 +2] = 81;
            RGB[135*3] = 0;
            RGB[135*3 +1] = 171;
            RGB[135*3 +2] = 79;
            RGB[136*3] = 0;
            RGB[136*3 +1] = 169;
            RGB[136*3 +2] = 77;
            RGB[137*3] = 0;
            RGB[137*3 +1] = 166;
            RGB[137*3 +2] = 76;
            RGB[138*3] = 0;
            RGB[138*3 +1] = 163;
            RGB[138*3 +2] = 74;
            RGB[139*3] = 0;
            RGB[139*3 +1] = 160;
            RGB[139*3 +2] = 73;
            RGB[140*3] = 0;
            RGB[140*3 +1] = 158;
            RGB[140*3 +2] = 72;
            RGB[141*3] = 0;
            RGB[141*3 +1] = 155;
            RGB[141*3 +2] = 71;
            RGB[142*3] = 0;
            RGB[142*3 +1] = 152;
            RGB[142*3 +2] = 70;
            RGB[143*3] = 0;
            RGB[143*3 +1] = 149;
            RGB[143*3 +2] = 68;
            RGB[144*3] = 0;
            RGB[144*3 +1] = 146;
            RGB[144*3 +2] = 66;
            RGB[145*3] = 0;
            RGB[145*3 +1] = 143;
            RGB[145*3 +2] = 65;
            RGB[146*3] = 0;
            RGB[146*3 +1] = 140;
            RGB[146*3 +2] = 62;
            RGB[147*3] = 0;
            RGB[147*3 +1] = 137;
            RGB[147*3 +2] = 60;
            RGB[148*3] = 0;
            RGB[148*3 +1] = 133;
            RGB[148*3 +2] = 58;
            RGB[149*3] = 0;
            RGB[149*3 +1] = 130;
            RGB[149*3 +2] = 55;
            RGB[150*3] = 0;
            RGB[150*3 +1] = 127;
            RGB[150*3 +2] = 52;
            RGB[151*3] = 0;
            RGB[151*3 +1] = 122;
            RGB[151*3 +2] = 47;
            RGB[152*3] = 0;
            RGB[152*3 +1] = 118;
            RGB[152*3 +2] = 43;
            RGB[153*3] = 0;
            RGB[153*3 +1] = 115;
            RGB[153*3 +2] = 39;
            RGB[154*3] = 0;
            RGB[154*3 +1] = 114;
            RGB[154*3 +2] = 38;
            RGB[155*3] = 0;
            RGB[155*3 +1] = 114;
            RGB[155*3 +2] = 39;
            RGB[156*3] = 0;
            RGB[156*3 +1] = 116;
            RGB[156*3 +2] = 40;
            RGB[157*3] = 0;
            RGB[157*3 +1] = 118;
            RGB[157*3 +2] = 43;
            RGB[158*3] = 0;
            RGB[158*3 +1] = 121;
            RGB[158*3 +2] = 45;
            RGB[159*3] = 0;
            RGB[159*3 +1] = 124;
            RGB[159*3 +2] = 47;
            RGB[160*3] = 0;
            RGB[160*3 +1] = 127;
            RGB[160*3 +2] = 50;
            RGB[161*3] = 0;
            RGB[161*3 +1] = 130;
            RGB[161*3 +2] = 52;
            RGB[162*3] = 0;
            RGB[162*3 +1] = 134;
            RGB[162*3 +2] = 54;
            RGB[163*3] = 0;
            RGB[163*3 +1] = 137;
            RGB[163*3 +2] = 56;
            RGB[164*3] = 0;
            RGB[164*3 +1] = 140;
            RGB[164*3 +2] = 57;
            RGB[165*3] = 0;
            RGB[165*3 +1] = 142;
            RGB[165*3 +2] = 58;
            RGB[166*3] = 0;
            RGB[166*3 +1] = 145;
            RGB[166*3 +2] = 58;
            RGB[167*3] = 0;
            RGB[167*3 +1] = 148;
            RGB[167*3 +2] = 58;
            RGB[168*3] = 0;
            RGB[168*3 +1] = 150;
            RGB[168*3 +2] = 58;
            RGB[169*3] = 0;
            RGB[169*3 +1] = 153;
            RGB[169*3 +2] = 57;
            RGB[170*3] = 0;
            RGB[170*3 +1] = 155;
            RGB[170*3 +2] = 56;
            RGB[171*3] = 0;
            RGB[171*3 +1] = 157;
            RGB[171*3 +2] = 55;
            RGB[172*3] = 0;
            RGB[172*3 +1] = 159;
            RGB[172*3 +2] = 52;
            RGB[173*3] = 0;
            RGB[173*3 +1] = 161;
            RGB[173*3 +2] = 49;
            RGB[174*3] = 0;
            RGB[174*3 +1] = 163;
            RGB[174*3 +2] = 45;
            RGB[175*3] = 1;
            RGB[175*3 +1] = 164;
            RGB[175*3 +2] = 36;
            RGB[176*3] = 3;
            RGB[176*3 +1] = 166;
            RGB[176*3 +2] = 21;
            RGB[177*3] = 5;
            RGB[177*3 +1] = 167;
            RGB[177*3 +2] = 7;
            RGB[178*3] = 9;
            RGB[178*3 +1] = 169;
            RGB[178*3 +2] = 0;
            RGB[179*3] = 17;
            RGB[179*3 +1] = 171;
            RGB[179*3 +2] = 0;
            RGB[180*3] = 31;
            RGB[180*3 +1] = 174;
            RGB[180*3 +2] = 0;
            RGB[181*3] = 47;
            RGB[181*3 +1] = 177;
            RGB[181*3 +2] = 0;
            RGB[182*3] = 61;
            RGB[182*3 +1] = 180;
            RGB[182*3 +2] = 0;
            RGB[183*3] = 70;
            RGB[183*3 +1] = 183;
            RGB[183*3 +2] = 0;
            RGB[184*3] = 78;
            RGB[184*3 +1] = 185;
            RGB[184*3 +2] = 0;
            RGB[185*3] = 86;
            RGB[185*3 +1] = 188;
            RGB[185*3 +2] = 0;
            RGB[186*3] = 93;
            RGB[186*3 +1] = 191;
            RGB[186*3 +2] = 0;
            RGB[187*3] = 101;
            RGB[187*3 +1] = 193;
            RGB[187*3 +2] = 0;
            RGB[188*3] = 108;
            RGB[188*3 +1] = 196;
            RGB[188*3 +2] = 0;
            RGB[189*3] = 115;
            RGB[189*3 +1] = 198;
            RGB[189*3 +2] = 0;
            RGB[190*3] = 122;
            RGB[190*3 +1] = 201;
            RGB[190*3 +2] = 0;
            RGB[191*3] = 129;
            RGB[191*3 +1] = 203;
            RGB[191*3 +2] = 0;
            RGB[192*3] = 136;
            RGB[192*3 +1] = 206;
            RGB[192*3 +2] = 0;
            RGB[193*3] = 143;
            RGB[193*3 +1] = 208;
            RGB[193*3 +2] = 0;
            RGB[194*3] = 150;
            RGB[194*3 +1] = 210;
            RGB[194*3 +2] = 0;
            RGB[195*3] = 157;
            RGB[195*3 +1] = 213;
            RGB[195*3 +2] = 0;
            RGB[196*3] = 164;
            RGB[196*3 +1] = 215;
            RGB[196*3 +2] = 0;
            RGB[197*3] = 170;
            RGB[197*3 +1] = 217;
            RGB[197*3 +2] = 0;
            RGB[198*3] = 177;
            RGB[198*3 +1] = 220;
            RGB[198*3 +2] = 0;
            RGB[199*3] = 184;
            RGB[199*3 +1] = 222;
            RGB[199*3 +2] = 0;
            RGB[200*3] = 191;
            RGB[200*3 +1] = 224;
            RGB[200*3 +2] = 0;
            RGB[201*3] = 198;
            RGB[201*3 +1] = 226;
            RGB[201*3 +2] = 0;
            RGB[202*3] = 205;
            RGB[202*3 +1] = 228;
            RGB[202*3 +2] = 0;
            RGB[203*3] = 212;
            RGB[203*3 +1] = 230;
            RGB[203*3 +2] = 0;
            RGB[204*3] = 221;
            RGB[204*3 +1] = 233;
            RGB[204*3 +2] = 0;
            RGB[205*3] = 230;
            RGB[205*3 +1] = 236;
            RGB[205*3 +2] = 0;
            RGB[206*3] = 235;
            RGB[206*3 +1] = 237;
            RGB[206*3 +2] = 0;
            RGB[207*3] = 236;
            RGB[207*3 +1] = 237;
            RGB[207*3 +2] = 0;
            RGB[208*3] = 236;
            RGB[208*3 +1] = 233;
            RGB[208*3 +2] = 0;
            RGB[209*3] = 236;
            RGB[209*3 +1] = 227;
            RGB[209*3 +2] = 0;
            RGB[210*3] = 236;
            RGB[210*3 +1] = 221;
            RGB[210*3 +2] = 0;
            RGB[211*3] = 237;
            RGB[211*3 +1] = 216;
            RGB[211*3 +2] = 0;
            RGB[212*3] = 237;
            RGB[212*3 +1] = 212;
            RGB[212*3 +2] = 0;
            RGB[213*3] = 237;
            RGB[213*3 +1] = 207;
            RGB[213*3 +2] = 0;
            RGB[214*3] = 237;
            RGB[214*3 +1] = 202;
            RGB[214*3 +2] = 0;
            RGB[215*3] = 237;
            RGB[215*3 +1] = 198;
            RGB[215*3 +2] = 0;
            RGB[216*3] = 237;
            RGB[216*3 +1] = 193;
            RGB[216*3 +2] = 0;
            RGB[217*3] = 237;
            RGB[217*3 +1] = 188;
            RGB[217*3 +2] = 0;
            RGB[218*3] = 237;
            RGB[218*3 +1] = 183;
            RGB[218*3 +2] = 0;
            RGB[219*3] = 237;
            RGB[219*3 +1] = 179;
            RGB[219*3 +2] = 0;
            RGB[220*3] = 237;
            RGB[220*3 +1] = 174;
            RGB[220*3 +2] = 0;
            RGB[221*3] = 237;
            RGB[221*3 +1] = 169;
            RGB[221*3 +2] = 0;
            RGB[222*3] = 237;
            RGB[222*3 +1] = 164;
            RGB[222*3 +2] = 0;
            RGB[223*3] = 237;
            RGB[223*3 +1] = 159;
            RGB[223*3 +2] = 0;
            RGB[224*3] = 237;
            RGB[224*3 +1] = 154;
            RGB[224*3 +2] = 0;
            RGB[225*3] = 237;
            RGB[225*3 +1] = 149;
            RGB[225*3 +2] = 0;
            RGB[226*3] = 237;
            RGB[226*3 +1] = 144;
            RGB[226*3 +2] = 0;
            RGB[227*3] = 236;
            RGB[227*3 +1] = 139;
            RGB[227*3 +2] = 0;
            RGB[228*3] = 236;
            RGB[228*3 +1] = 134;
            RGB[228*3 +2] = 0;
            RGB[229*3] = 235;
            RGB[229*3 +1] = 129;
            RGB[229*3 +2] = 0;
            RGB[230*3] = 235;
            RGB[230*3 +1] = 124;
            RGB[230*3 +2] = 0;
            RGB[231*3] = 234;
            RGB[231*3 +1] = 119;
            RGB[231*3 +2] = 0;
            RGB[232*3] = 234;
            RGB[232*3 +1] = 113;
            RGB[232*3 +2] = 0;
            RGB[233*3] = 233;
            RGB[233*3 +1] = 108;
            RGB[233*3 +2] = 0;
            RGB[234*3] = 232;
            RGB[234*3 +1] = 103;
            RGB[234*3 +2] = 0;
            RGB[235*3] = 231;
            RGB[235*3 +1] = 97;
            RGB[235*3 +2] = 0;
            RGB[236*3] = 230;
            RGB[236*3 +1] = 91;
            RGB[236*3 +2] = 0;
            RGB[237*3] = 229;
            RGB[237*3 +1] = 85;
            RGB[237*3 +2] = 0;
            RGB[238*3] = 228;
            RGB[238*3 +1] = 80;
            RGB[238*3 +2] = 0;
            RGB[239*3] = 227;
            RGB[239*3 +1] = 73;
            RGB[239*3 +2] = 0;
            RGB[240*3] = 225;
            RGB[240*3 +1] = 67;
            RGB[240*3 +2] = 0;
            RGB[241*3] = 224;
            RGB[241*3 +1] = 61;
            RGB[241*3 +2] = 0;
            RGB[242*3] = 222;
            RGB[242*3 +1] = 53;
            RGB[242*3 +2] = 0;
            RGB[243*3] = 220;
            RGB[243*3 +1] = 45;
            RGB[243*3 +2] = 0;
            RGB[244*3] = 218;
            RGB[244*3 +1] = 32;
            RGB[244*3 +2] = 0;
            RGB[245*3] = 216;
            RGB[245*3 +1] = 17;
            RGB[245*3 +2] = 1;
            RGB[246*3] = 214;
            RGB[246*3 +1] = 5;
            RGB[246*3 +2] = 3;
            RGB[247*3] = 212;
            RGB[247*3 +1] = 0;
            RGB[247*3 +2] = 4;
            RGB[248*3] = 210;
            RGB[248*3 +1] = 0;
            RGB[248*3 +2] = 10;
            RGB[249*3] = 208;
            RGB[249*3 +1] = 0;
            RGB[249*3 +2] = 21;
            RGB[250*3] = 206;
            RGB[250*3 +1] = 0;
            RGB[250*3 +2] = 31;
            RGB[251*3] = 204;
            RGB[251*3 +1] = 0;
            RGB[251*3 +2] = 37;
            RGB[252*3] = 202;
            RGB[252*3 +1] = 0;
            RGB[252*3 +2] = 40;
            RGB[253*3] = 200;
            RGB[253*3 +1] = 0;
            RGB[253*3 +2] = 43;
            RGB[254*3] = 197;
            RGB[254*3 +1] = 0;
            RGB[254*3 +2] = 44;
            RGB[255*3] = 194;
            RGB[255*3 +1] = 0;
            RGB[255*3 +2] = 45;
            break;
        default:
            for (int i=0; i<256; i++)
            {
                RGB[i*3] = (unsigned char) i;
                RGB[i*3 + 1] = (unsigned char) i;
                RGB[i*3 + 2] = (unsigned char) i;
            }
            break;
    }
}

void Image2dframe::OnNextClicked(wxCommandEvent& event)
{
   cmpnumber = cmpnumber + dcmp; 
   if(cmpnumber > maxcmp) cmpnumber = maxcmp;
   wxCommandEvent parevent(SelectCmp, GetId());
   parevent.SetEventObject(this);
   parevent.SetInt(cmpnumber);
   // Send event to App
   ProcessWindowEvent(parevent);
}

void Image2dframe::OnCMPinterval(wxCommandEvent& event)
{
    long val;

    val  = wxGetNumberFromUser(	_("Changes the increment of the next and previous tools"), _("Enter a value"), _("Change gather interval"), this->getDcmp(), 1, this->getMaxcmp(), this, wxDefaultPosition);
    if(val != -1){
        this->setDcmp(val);
    }
}

void Image2dframe::OnPreviousClicked(wxCommandEvent& event)
{
   cmpnumber = cmpnumber - dcmp; 
   if(cmpnumber < 0) cmpnumber = 0;
   wxCommandEvent parevent(SelectCmp, GetId());
   parevent.SetEventObject(this);
   parevent.SetInt(cmpnumber);
   // Send event to App
   ProcessWindowEvent(parevent);
}


void Image2dframe::OnSavepicks(wxCommandEvent& event)
{
    if(nlayers > 0){
        wxCommandEvent parevent(SaveEvent, GetId());
        parevent.SetEventObject(this);
        parevent.SetInt(0);
        // Send event to App
        ProcessWindowEvent(parevent);
    }
}

void Image2dframe::setStatusbar(int cmp, float x, float y, float val)
{
    char label[48];
    snprintf(label, 48, "CMP: %d, X: %.2f, Z: %.2f, Val: %f", cmp, x, y, val);
    StatusBar1->SetStatusText(_(label));
}

void Image2dframe::createPicks()
{
   Picks *newpicks ; 
   newpicks = new Picks(MAXPOS, this->getMaxcmp()+1, this->getN1(), this->getD1(), this->getO1(), PICK_HORIZONTAL);
   picks.push_back(newpicks);
   nlayers++;
}

void Image2dframe::createPicks(int direction)
{
    Picks *newpicks ; 
    if(direction == PICK_HORIZONTAL){
        newpicks = new Picks(MAXPOS, this->getMaxcmp()+1, this->getN1(), this->getD1(), this->getO1(), PICK_HORIZONTAL);
    }else{
        newpicks = new Picks(MAXPOS, this->getMaxcmp()+1, this->getN2(), this->getD2(), this->getO2(), PICK_VERTICAL);
    }
   picks.push_back(newpicks);
   nlayers++;
}

void Image2dframe::Plotcrosshair(wxDC &dc, int w, int h)
{
        float v0,v1, t0, t1;
        v0=zoom->Getx0();
        v1=zoom->Getx1();
        t0=zoom->Gety0();
        t1=zoom->Gety1();
        float ax;
        float ay;
        ay=(t1 - t0)/(h-1);
        ax=(v1 - v0)/(w-1);

        pos[0].x=(int) ((crosshair_pt[0] - v0)/ax);
        pos[0].y=(int) ((crosshair_pt[1] - t0)/ay);

        wxPen myCrossPen(*wxWHITE,2,wxDOT_DASH);
        dc.SetPen(myCrossPen);
        dc.CrossHair(pos[0]);
}


void Image2dframe::Plotpicks(wxDC &dc, int w, int h)
{
        int current_cmp=this->getCmpnumber();
        int *np;
        int nump;
        np=picks[layer]->Getnpicks();
        rockseis::Point2D<float> *points;
        points=picks[layer]->Getpicks();
        int maxpicks=picks[layer]->Getmaxpicks();
        float v0,v1, t0, t1;
        v0=zoom->Getx0();
        v1=zoom->Getx1();
        t0=zoom->Gety0();
        t1=zoom->Gety1();
        float ax;
        float ay;
        ay=(t1 - t0)/(h-1);
        ax=(v1 - v0)/(w-1);

        int i;
        //Find previous and next picks
        int previous,next;
        previous=current_cmp;
        nump=0;
        while(previous>0){
            previous--;
            nump=np[previous];
            if(nump) break;
        }
        if(nump){
            for(i=0; i<nump; i++)
            {
                pos[i].x=(int) ((points[previous*maxpicks + i].x - v0)/ax);
                pos[i].y=(int) ((points[previous*maxpicks + i].y - t0)/ay);
            }
                // Create a red pen 2 pixels wide drawing a dotted line
                wxPen myPen(*wxRED,2,wxDOT_DASH);
                // Tell dc to start using this pen to draw.
                dc.SetPen( myPen );
                dc.DrawLines(nump, pos);
        }

        next=current_cmp;
        nump=0;
        while(next<this->getMaxcmp()){
            next++;
            nump=np[next];
            if(nump) break;
        }
        if(nump){
            for(i=0; i<nump; i++)
            {
                pos[i].x=(int) ((points[next*maxpicks + i].x - v0)/ax);
                pos[i].y=(int) ((points[next*maxpicks + i].y - t0)/ay);
            }
                // Create a red pen 2 pixels wide drawing a dotted line
                wxPen myPen(*wxGREEN,2,wxDOT_DASH);
                // Tell dc to start using this pen to draw.
                dc.SetPen( myPen );
                dc.DrawLines(nump, pos);
        }

        nump=np[current_cmp];
        for(i=0; i<nump; i++)
        {
            pos[i].x=(int) ((points[current_cmp*maxpicks + i].x - v0)/ax);
            pos[i].y=(int) ((points[current_cmp*maxpicks + i].y - t0)/ay);
        }
        if(nump){
        // Create a green pen 3 pixels wide drawing a solid line
            wxPen myWhitePen(*wxWHITE,2,wxSOLID);
        // Tell dc to start using this pen to draw.
            dc.SetPen( myWhitePen );
            for(i=0; i<nump; i++){
                // Draw a Point
                dc.DrawCircle(pos[i].x, pos[i].y, 3);
            }
            dc.DrawLines(nump, pos);
        }
}

void Image2dframe::setZoom(float o1, float max1, float o2, float max2, int ix0, int nx, int iy0, int ny)
{
            zoom->Setx0(o1);
            zoom->Setix0(ix0);
            zoom->Setnx(nx);
            zoom->Setx1(max1);
            zoom->Sety0(o2);
            zoom->Setiy0(iy0);
            zoom->Setny(ny);
            zoom->Sety1(max2);
}

void Image2dframe::OnPick(wxCommandEvent& event)
{
    if(ToolBarItem4->IsToggled()){
        ToolBar1->ToggleTool(idToolzoom, false);
    }
}

void Image2dframe::OnZoom(wxCommandEvent& event)
{
    if(ToolBarItem3->IsToggled()){
        ToolBar1->ToggleTool(idToolpick, false);
    }
}
