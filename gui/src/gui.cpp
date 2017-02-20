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
	MainFrame *frame = new MainFrame( "Hello World", wxPoint(50, 50), wxSize(600, 600) );
	frame->Show( true );
	return true;
}

