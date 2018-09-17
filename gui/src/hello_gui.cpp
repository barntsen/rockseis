#include "mainframe.h"

// Class declaration of the app
class MyApp: public wxApp
{
public:
	virtual bool OnInit();
};

wxBEGIN_EVENT_TABLE(MainFrame, wxFrame)
	EVT_MENU(ID_Hello,   MainFrame::OnHello)
	EVT_MENU(wxID_EXIT,  MainFrame::OnExit)
	EVT_MENU(wxID_ABOUT, MainFrame::OnAbout)
wxEND_EVENT_TABLE()


// Main implementation of app
wxIMPLEMENT_APP(MyApp);


// Class implementation of the app
bool MyApp::OnInit()
{
	wxString hello = "Hello World";
	MainFrame *frame1 = new MainFrame( hello, wxPoint(50, 50), wxSize(600, 600) );
	frame1->Show( true );
    
	return true;
}

