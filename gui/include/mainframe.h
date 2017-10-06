// Include statements
#include <wx/wxprec.h>
#include <wx/app.h>
#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif


class MainFrame: public wxFrame {
public:
	MainFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
private:
	void OnHello(wxCommandEvent& event);
	void OnExit(wxCommandEvent& event);
	void OnAbout(wxCommandEvent& event);
	wxDECLARE_EVENT_TABLE();
};


enum
{
    ID_Hello = 1
};

