#include "image2dframe.h"

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
const long Image2dframe::ID_MENUITEM1 = wxNewId();
const long Image2dframe::ID_MENUITEM2 = wxNewId();
const long Image2dframe::ID_MENUITEM3 = wxNewId();
const long Image2dframe::ID_MENUITEM4 = wxNewId();
const long Image2dframe::ID_MENUITEM5 = wxNewId();
const long Image2dframe::ID_MENUITEM6 = wxNewId();
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
    menubarset = false;

    cmpnumber = 0;
    dcmp = 1;
    maxcmp = 1;

    nlayers=0;
    layer=0;

    displaycrosshair = false;
    getcrosshair = false;

	//(*Initialize(Image2dframe)
	wxFlexGridSizer* FlexGridSizer1;

	Create(parent, id, _("2D image window"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE, _T("id"));
	SetClientSize(wxDefaultSize);
	Move(wxDefaultPosition);
	FlexGridSizer1 = new wxFlexGridSizer(2, 3, 0, 0);
	this->SetSizer(FlexGridSizer1);
	FlexGridSizer1->AddGrowableCol(1);
	FlexGridSizer1->AddGrowableRow(1);
	LeftCorner = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1->Add(LeftCorner, 1, wxALL|wxEXPAND, 1);
	Xaxis = new wxPanel(this, ID_PANEL2, wxDefaultPosition, wxSize(600,50), wxTAB_TRAVERSAL, _T("ID_PANEL2"));
	FlexGridSizer1->Add(Xaxis, 1, wxALL|wxEXPAND, 1);
	RightCorner = new wxPanel(this, ID_PANEL3, wxDefaultPosition, wxSize(50,50), wxTAB_TRAVERSAL, _T("ID_PANEL3"));
	FlexGridSizer1->Add(RightCorner, 1, wxALL|wxEXPAND, 1);
	Zaxis = new wxPanel(this, ID_PANEL4, wxDefaultPosition, wxSize(50,600), wxTAB_TRAVERSAL, _T("ID_PANEL4"));
	FlexGridSizer1->Add(Zaxis, 1, wxALL|wxEXPAND, 1);
	Imagewindow = new wxPanel(this, ID_IMAGEWINDOW1, wxDefaultPosition, wxSize(600,600), wxWANTS_CHARS, _T("ID_IMAGEWINDOW1"));
	FlexGridSizer1->Add(Imagewindow, 1, wxALL|wxEXPAND, 1);
    StatusBar1 = new wxStatusBar(this, ID_STATUSBAR1, wxSIMPLE_BORDER, _T("ID_STATUSBAR1"));
    int __wxStatusBarWidths_1[1] = { -10 };
    int __wxStatusBarStyles_1[1] = { wxSB_NORMAL };
    StatusBar1->SetFieldsCount(1,__wxStatusBarWidths_1);
    StatusBar1->SetStatusStyles(1,__wxStatusBarStyles_1);
    SetStatusBar(StatusBar1);

    FlexGridSizer1->Fit(this);
	FlexGridSizer1->SetSizeHints(this);

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
    getRgb(color, &RGB[0]);
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

    wxClientDC dc( Imagewindow );
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
    wxClientDC dc( Imagewindow );
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
            zoom->Sety0(o2);
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

	    wxClientDC dc( Imagewindow );
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
	wxClientDC dc( Imagewindow );
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
        getRgb(color, &RGB[0]);
        this->LoadImage(zoom->Getix0(), zoom->Getnx(), zoom->Getiy0(), zoom->Getny());
        Refresh();
    }

    if(event.GetKeyCode() == wxKeyCode('R')){
        color++; 
        color = color%NCOLORS;
        getRgb(color, &RGB[0]);
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
	    Connect(idToolLoad,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&Image2dframe::OnLoadpicks);
    }
}

void Image2dframe::createMenubar()
{
    if(!this->getMenubarset()){
        MenuBar = new wxMenuBar();
        Menu1 = new wxMenu();
        MenuItem1 = new wxMenuItem(Menu1, ID_MENUITEM1, _("Load picks"), wxEmptyString, wxITEM_NORMAL);
        Menu1->Append(MenuItem1);
        MenuItem2 = new wxMenuItem(Menu1, ID_MENUITEM2, _("Save picks"), wxEmptyString, wxITEM_NORMAL);
        Menu1->Append(MenuItem2);
        MenuItem3 = new wxMenuItem(Menu1, ID_MENUITEM3, _("Clear picks"), _("Clear all picks on the dataset"), wxITEM_NORMAL);
        Menu1->Append(MenuItem3);
        MenuItem4 = new wxMenuItem(Menu1, ID_MENUITEM4, _("Close"), wxEmptyString, wxITEM_NORMAL);
        Menu1->Append(MenuItem4);
        MenuBar->Append(Menu1, _("File"));

        Menu2 = new wxMenu();
        MenuItem5 = new wxMenuItem(Menu2, ID_MENUITEM5, _("Mute left/above"), _("Mute events to the left of vertical picks or above horizontal picks"), wxITEM_NORMAL);
        Menu2->Append(MenuItem5);
        MenuItem6 = new wxMenuItem(Menu2, ID_MENUITEM6, _("Mute right/below"), _("Mute events to the right of vertical picks or below horizontal picks"), wxITEM_NORMAL);
        Menu2->Append(MenuItem6);
        MenuBar->Append(Menu2, _("Mute"));
        SetMenuBar(MenuBar);
        this->setMenubarset(true);


	    Connect(ID_MENUITEM1,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnLoadpicks);
	    Connect(ID_MENUITEM2,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnSavepicks);
	    Connect(ID_MENUITEM3,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnClearpicks);
	    Connect(ID_MENUITEM4,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnClose);
	    Connect(ID_MENUITEM5,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnMute_la);
	    Connect(ID_MENUITEM6,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&Image2dframe::OnMute_rb);
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

void Image2dframe::OnLoadpicks(wxCommandEvent& event)
{
    if(nlayers > 0){
        wxCommandEvent parevent(LoadEvent, GetId());
        parevent.SetEventObject(this);
        parevent.SetInt(0);
        // Send event to App
        ProcessWindowEvent(parevent);
    }
}

void Image2dframe::OnClearpicks(wxCommandEvent& event)
{

    if(nlayers > 0){
        wxMessageDialog dialog( this, wxT("This will clear all picks in the data.\nDo you really want to proceed?"),
                wxT("Clear all picks"),
                wxNO_DEFAULT | wxYES_NO);

        switch( dialog.ShowModal() )
        {
            case wxID_YES:
                picks[layer]->Clearpicks();
                Refresh();
                break;
            case wxID_NO:
                break;
            default:
                break;
        }
    }
}

void Image2dframe::OnMute_la(wxCommandEvent& event)
{
    if(nlayers > 0){
        wxCommandEvent parevent(MuteEvent, GetId());
        parevent.SetEventObject(this);
        parevent.SetInt(0);
        // Send event to App
        ProcessWindowEvent(parevent);
    }
}

void Image2dframe::OnMute_rb(wxCommandEvent& event)
{
    if(nlayers > 0){
        wxCommandEvent parevent(MuteEvent, GetId());
        parevent.SetEventObject(this);
        parevent.SetInt(1);
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

void Image2dframe::createPicks(int direction, long n, float d, float o)
{
    Picks *newpicks ; 
    if(direction == PICK_HORIZONTAL){
        newpicks = new Picks(MAXPOS, this->getMaxcmp()+1, n, d, o, PICK_HORIZONTAL);
    }else{
        newpicks = new Picks(MAXPOS, this->getMaxcmp()+1, n, d, o, PICK_VERTICAL);
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
        if(nump>1){
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
        if(nump>1){
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
            if(nump>1){
               dc.DrawLines(nump, pos);
            }
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
