#include "image2dframe.h"
#include "file.h"
#include "utils.h"
#include "geometry.h"
#include <wx/filedlg.h>
#include <wx/msgdlg.h>
#include <args.hxx>

// Class declaration of the app
class MyApp: public wxApp
{
public:
	virtual bool OnInit();
	void readCIP(wxCommandEvent& event);
	void Save(wxCommandEvent& event);
	void Load(wxCommandEvent& event);
	void Cross(wxCommandEvent& event);
	int FilterEvent(wxEvent& event);


private:
    void OnZoframeDestroy(wxWindowDestroyEvent& event);
    void OnCipframeDestroy(wxWindowDestroyEvent& event);
    std::string infile, picksfile, mutedfile;
    size_t n1,n3,n4,n6,ntot;
    float d1,d3,d4,d6;
    float o1,o3,o4,o6;
    float *cipdata=NULL;
    float *zodata=NULL;
    double *ddata=NULL;
    int dsize;
    int ix;
    Image2dframe *cipframe;
    Image2dframe *zoframe;
};

// Main implementation of app
wxIMPLEMENT_APP(MyApp);

wxDEFINE_EVENT(SelectCmp, wxCommandEvent);
wxDEFINE_EVENT(SaveEvent, wxCommandEvent);
wxDEFINE_EVENT(LoadEvent, wxCommandEvent);
wxDEFINE_EVENT(Crosshair, wxCommandEvent);

// Class implementation of the app
bool MyApp::OnInit()
{
    args::ArgumentParser parser("Program to pick a mute on IMAGE type files.", "");
    parser.LongPrefix("-");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> parinput(parser, "FILENAME", "Input file", {"i"});
    args::ValueFlag<std::string> parpicksfile(parser, "FILENAME", "Output picks file", {"picks"});
    args::ValueFlag<std::string> paroutputfile(parser, "FILENAME", "Output muted file", {"mutedfile"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cerr << parser;
        return false;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return false;
    }
    if (parinput){
         infile = args::get(parinput);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing input file -i");
    }

    if (parpicksfile){
         picksfile = args::get(parpicksfile);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing picks file -picks");
    }

    if (paroutputfile){
        mutedfile = args::get(paroutputfile);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing picks file -picks");
    }



    std::shared_ptr<rockseis::File> in (new rockseis::File());

    bool status;
    status = in->input(infile);
    if(status == FILE_ERR){
        rockseis::rs_error("rsrss2rsf: Error opening file for input.");
        std::cerr << status << std::endl;
    }

    // Get parameters
    n1 = in->getN(1);
    n3 = in->getN(3);
    n4 = in->getN(4);
    n6 = in->getN(6);

    d1 = (float) in->getD(1);
    d3 = (float) in->getD(3);
    d4 = (float) in->getD(4);
    d6 = (float) in->getD(6);

    o1 = (float) in->getO(1);
    o3 = (float) in->getO(3);
    o4 = (float) in->getO(4);
    o6 = (float) in->getO(6);

    if(n3 == 0 || n4 == 0){
        rockseis::rs_error("Image is not the correct format, or is not 2D");
    }
    ntot = in->getNtot();
    if(n1*n3*n4 != ntot) rockseis::rs_error("Image data is not 2 dimensional");

    cipdata = (float *) calloc(n3*n4, sizeof(float));
    if(cipdata == NULL) rockseis::rs_error("Failed to allocate memory for the cip data.");

    zodata = (float *) calloc(n1*n3, sizeof(float));
    if(zodata == NULL) rockseis::rs_error("Failed to allocate memory for the zo image.");
    dsize = in->getData_format();

    int Nheader = in->getNheader();
    if(Nheader > 0) rockseis::rs_error("Image should not contain coordinate headers.");

    rockseis::Index I2D(n1,n3,n4,n6); 

    //Get Zo data
    in->seekg(I2D(0,0,(n4-1)/2, 0)*dsize + in->getStartofdata());
    switch(dsize){
        case 4:
            in->read(&zodata[0], n1*n3);
            break;
        case 8:
            ddata = (double *) calloc(n1*n3, sizeof(double));
            if(ddata == NULL) rockseis::rs_error("Failed to allocate memory for the image.");
            in->read(&zodata[0], n1*n3);
            for (size_t i=0; i < n1*n3; i++){
                zodata[i] = (float) ddata[i];
            }
            free(ddata);
            break;
        default:
            rockseis::rs_error("Unsupported esize.");
            break;
    }

    //Get Cip data
    switch(dsize){
        case 4:
            // Seek to data points
            for (size_t i=0; i < n3; i++){
                for (size_t j=0; j < n4; j++){
                    in->seekg(I2D(ix,i,j, 0)*dsize + in->getStartofdata());
                    in->read(&cipdata[i*n4 + j], 1);
                }
            }
            break;
        case 8:
            ddata = (double *) calloc(n3*n4, sizeof(double));
            if(ddata == NULL) rockseis::rs_error("Failed to allocate memory for the cip data.");
            for (size_t i=0; i < n3; i++){
                for (size_t j=0; j < n4; j++){
                    in->seekg(I2D(ix,i,j, 0)*dsize + in->getStartofdata());
                    in->read(&cipdata[i*n4 + j], 1);
                }
            }
            for (size_t i=0; i < n1*n4; i++){
                cipdata[i] = (float) ddata[i];
            }
            free(ddata);
            break;
        default:
            rockseis::rs_error("Unsupported esize.");
            break;
    }
    in->close();

    //(*AppInitialize
    bool wxsOK = true;
    wxInitAllImageHandlers();
    if ( wxsOK ){

	    cipframe = new Image2dframe(n4, d4, o4, n3, d3, o3, cipdata, 0);
	    cipframe->createToolbar();
	    cipframe->SetLabel (wxT("Picking window"));
	    cipframe->setMaxcmp(n1-1);
        cipframe->createPicks(PICK_VERTICAL, n3, d3, o3);
        cipframe->setGetcrosshair(true);
        cipframe->Connect(wxEVT_DESTROY, wxWindowDestroyEventHandler(MyApp::OnCipframeDestroy),NULL, this);
	    cipframe->Show( true );
	    Bind(SelectCmp, &MyApp::readCIP, this, cipframe->GetId());
	    Bind(SaveEvent, &MyApp::Save, this, cipframe->GetId());
	    Bind(LoadEvent, &MyApp::Load, this, cipframe->GetId());
	    Bind(Crosshair, &MyApp::Cross, this, cipframe->GetId());

	    zoframe = new Image2dframe(n1, d1, o1, n3, d3, o3, zodata, 0);
	    zoframe->SetLabel (wxT("Zero offset image window"));
        zoframe->setDisplaycrosshair(true);
        zoframe->Connect(wxEVT_DESTROY, wxWindowDestroyEventHandler(MyApp::OnZoframeDestroy),NULL, this);
	    zoframe->Show( true );
    }

    return wxsOK;
}

void MyApp::readCIP(wxCommandEvent& event)
{
    int ix = event.GetInt();
    if(ix < 0) ix = 0;
    if(ix > n1-1) ix = n1-1;

    rockseis::Index I2D(n1,n3,n4,n6); 
    
    std::shared_ptr<rockseis::File> in (new rockseis::File());

    bool status;
    status = in->input(infile);
    if(status == FILE_ERR){
        rockseis::rs_error("rsrss2rsf: Error opening file for input.");
        std::cerr << status << std::endl;
    }

    //Get Cip data
    switch(dsize){
        case 4:
            // Seek to data points
            for (size_t i=0; i < n3; i++){
                for (size_t j=0; j < n4; j++){
                    in->seekg(I2D(ix,i,j, 0)*dsize + in->getStartofdata());
                    in->read(&cipdata[i*n4 + j], 1);
                }
            }
            break;
        case 8:
            ddata = (double *) calloc(n3*n4, sizeof(double));
            if(ddata == NULL) rockseis::rs_error("Failed to allocate memory for the cip data.");
            for (size_t i=0; i < n3; i++){
                for (size_t j=0; j < n4; j++){
                    in->seekg(I2D(ix,i,j, 0)*dsize + in->getStartofdata());
                    in->read(&cipdata[i*n4 + j], 1);
                }
            }
            for (size_t i=0; i < n1*n4; i++){
                cipdata[i] = (float) ddata[i];
            }
            free(ddata);
            break;
        default:
            rockseis::rs_error("Unsupported esize.");
            break;
    }
    cipframe->ComputeClip();
    cipframe->LoadImage(0, n4, 0, n3);
    cipframe->setStatusbar(ix, 0., 0., 0.); 
    zoframe->setCmpnumber(ix);
    cipframe->Refresh();
    zoframe->Refresh();
    in->close();
}

int MyApp::FilterEvent(wxEvent& event)
{
	if (event.GetEventType() == wxEVT_KEY_UP)

	{
		event.Skip(true);
	}

	return -1;
}

void MyApp::Save(wxCommandEvent& event)
{
    std::vector<Picks*> *picks;
    picks = cipframe->getPicks();
    int *npicks; 
    npicks = (*picks)[0]->Getnpicks();
    bool haspicks=false;
    for(int i=0; i< n1; i++)
    {
        if(npicks[i] > 0){
            haspicks=true;
        }
    }
    if(haspicks){
        if((*picks)[0]->Savepicks(picksfile) != PICKS_OK) {
            wxMessageBox(_("Something went wrong when saving the picks file."), _("Error saving picks"), wxICON_INFORMATION);
            return;
        }
    }
}

void MyApp::Load(wxCommandEvent& event)
{
    std::string filename;
    std::vector<Picks*> *picks;
    picks = cipframe->getPicks();
    int *npicks; 
    npicks = (*picks)[0]->Getnpicks();

    bool haspicks=false;
    for(int i=0; i< n1; i++)
    {
        if(npicks[i] > 0){
            haspicks=true;
        }
    }

    if (haspicks && (*picks)[0]->Getchanged())
    {
        if (wxMessageBox(_("Current picks have not been saved! Proceed?"), _("Please confirm"),
                         wxICON_QUESTION | wxYES_NO, cipframe) == wxNO )
            return;
        //else: proceed asking to the user the new file to open
    }

    wxFileDialog openFileDialog(cipframe, _("Open Picks file"), "", "", "Pick files (*.rss)|*.rss", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;     // the user changed idea...

    // proceed loading the file chosen by the user;
    filename = openFileDialog.GetPath();
    if((*picks)[0]->Loadpicks(filename) != PICKS_OK) {
        wxMessageBox(_("Something went wrong when reading from the picks file."), _("Error loading picks"), wxICON_INFORMATION);
	(*picks)[0]->Clearpicks();
        return;
    }
    cipframe->Refresh();
}

void MyApp::Cross(wxCommandEvent& event)
{
    float *pt = (float *) event.GetClientData();
    float *ptzo = zoframe->getCrosshair_pt();
    ptzo[0] = pt[0]*d1 + o1;
    ptzo[1] = pt[1];
    zoframe->Refresh();
}

void MyApp::OnZoframeDestroy(wxWindowDestroyEvent& event)
{
    cipframe->Close();
    exit(0);
}

void MyApp::OnCipframeDestroy(wxWindowDestroyEvent& event)
{
    zoframe->Close();
    exit(0);
}
