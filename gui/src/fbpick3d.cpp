#include "image2dframe.h"
#include "file.h"
#include "sort.h"
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
	void readGath(wxCommandEvent& event);
	void Save(wxCommandEvent& event);
	void Load(wxCommandEvent& event);
	int FilterEvent(wxEvent& event);
    void transposeData();

private:
    Image2dframe *frame;
    std::shared_ptr<rockseis::Sort<float>> Sort;
    std::shared_ptr<rockseis::Data3D<float>> gather3d;
    std::string infile;
    std::string picksfile;
    std::string tdatafile;
    int sort_type;
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
    args::ArgumentParser parser("Program to pick first break on DATA type files.", "");
    parser.LongPrefix("-");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> parinput(parser, "FILENAME", "Input file", {"i"});
    args::ValueFlag<std::string> parpicksfile(parser, "FILENAME", "Output picks file", {"picks"});
    args::ValueFlag<std::string> partdatafile(parser, "FILENAME", "Output traveltime data file", {"tdata"});
    args::ValueFlag<int> parsort(parser, "Integer", "Sort order: 0 - Shot gathers (default), 1 - Receiver gathers", {"sort"});
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

    if (partdatafile){
         tdatafile = args::get(partdatafile);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing output traveltime data file -tdata");
    }

    if (parsort){
         sort_type = args::get(parsort);
    }else{
         sort_type = 0;
    }


    Sort = std::make_shared<rockseis::Sort<float>>();
    Sort->setDatafile(infile);
    switch(sort_type){
        case 0:
            Sort->createShotmap(infile); 
            break;
        case 1:
            Sort->createReceivermap(infile); 
            break;
        default:
            Sort->createShotmap(infile); 
            break;
    }


    gather3d = Sort->get3DGather(0);
    this->transposeData();

    // Get parameters
	wxString pick = "Gather";
    frame = new Image2dframe(gather3d->getNtrace(), 1.0, 0, gather3d->getNt(), gather3d->getDt(), gather3d->getOt(), gather3d->getData(), 0);
    frame->createToolbar();
    frame->SetLabel (wxT("Picking window"));
    frame->setMaxcmp(Sort->getNensemb()-1);
    frame->createPicks();
    frame->Show( true );

	Bind(SelectCmp, &MyApp::readGath, this, frame->GetId());
	Bind(SaveEvent, &MyApp::Save, this, frame->GetId());
	Bind(LoadEvent, &MyApp::Load, this, frame->GetId());
	return true;
}

void MyApp::readGath(wxCommandEvent& event)
{
    int ix = event.GetInt();
    if(ix < 0) ix = 0;
    if(ix > Sort->getNensemb()-1) ix = Sort->getNensemb()-1;
    gather3d.reset();
    gather3d = Sort->get3DGather(ix);
    if (gather3d == NULL){
        rockseis::rs_error("Error getting gather.");
        exit(1);
    }
    this->transposeData();

    frame->setN1(gather3d->getNtrace());
    frame->setN2(gather3d->getNt());
    frame->setImagedata(gather3d->getData());
    frame->ComputeClip();
    frame->setZoom(frame->getO1(), (frame->getN1()-1)*frame->getD1() + frame->getO1(), frame->getO2(), (frame->getN2()-1)*frame->getD2() + frame->getO2(), 0, frame->getN1(), 0, frame->getN2());
    frame->LoadImage(0, gather3d->getNtrace(), 0, gather3d->getNt());
    frame->setStatusbar(ix, 0., 0., 0.); 
    frame->Refresh();
}

void MyApp::transposeData()
{
    int nt = gather3d->getNt();
    int ntr = gather3d->getNtrace();
    rockseis::Index Idata(nt,ntr);
    rockseis::Index Itransp(ntr,nt);
    float *fdata = gather3d->getData();
    float *ftrans = (float *) calloc(nt*ntr, sizeof(float));
    for(int i=0; i < nt; i++){
        for (int j=0; j<ntr; j++){
            ftrans[Itransp(j,i)] = fdata[Idata(i,j)];
        }
    }
    for(int i=0; i < nt*ntr; i++){
            fdata[i] = ftrans[i];
    }
    free(ftrans);
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
    std::shared_ptr<rockseis::Data3D<float>> wrkgather;
    std::shared_ptr<rockseis::Data3D<float>> outgather;
    rockseis::Point3D<float> *scoords3d;
    std::vector<Picks*> *picks;
    picks = frame->getPicks();
    int *npicks; 
    int maxpicks=(*picks)[0]->Getmaxpicks();
    npicks = (*picks)[0]->Getnpicks();
    rockseis::Point2D<float> *points;
    points=(*picks)[0]->Getpicks();
    rockseis::Index Idata(2*maxpicks,Sort->getNensemb());
    rockseis::Index Ipicks(maxpicks,Sort->getNensemb());
    float *vrms;
    float *data;
    bool haspicks=false;
    for(int i=0; i< Sort->getNensemb(); i++)
    {
        if(npicks[i] > 0){
            haspicks=true;
        }
    }
    if(haspicks){
	    if((*picks)[0]->Savepicks(picksfile) != PICKS_OK){
		    wxMessageBox(_("Something went wrong when saving the picks file."), _("Error saving picks"), wxICON_INFORMATION);
		    return;
	    }
	    /*
        outgather = std::make_shared<rockseis::Data3D<float>>(Sort->getNensemb(), 2*(*picks)[0]->Getmaxpicks(), 1.0, 0.0);
        outgather->setFile(picksfile);
        scoords3d = (outgather->getGeom())->getScoords();
        data = outgather->getData();
        // Save only picks in a special format
        for(int i=0; i< Sort->getNensemb(); i++)
        {
            // Save number of picks in sx
            scoords3d[i].x = npicks[i];
            if(npicks[i]){
                // Save coordinates in trace
                for(int j=0; j<npicks[i]; j++){
                    data[Idata(2*j,i)] = points[Ipicks(j,i)].x;
                    data[Idata(2*j+1,i)] = points[Ipicks(j,i)].y;
                }
            }
        }
        outgather->write();
        outgather->close();
	*/

        // Save data in another file
        wrkgather = Sort->get3DGather(0);
        bool open = false;
        for(int i=0; i< Sort->getNensemb(); i++)
        {
            if(npicks[i] > 0){
                //Get gather3d
                wrkgather.reset();
                wrkgather = Sort->get3DGather(i);
                outgather.reset();
                outgather = std::make_shared<rockseis::Data3D<float>>(wrkgather->getNtrace(), 1, 1.0, 0.0);
                outgather->copyCoords(wrkgather);

                data = outgather->getData();
                //Compute interpolated picks
                (*picks)[0]->Interp(i);
                vrms  = (*picks)[0]->Getvrms();
                for(int j=0; j<outgather->getNtrace(); j++){
                    data[j] = vrms[j];
                }
                outgather->setFile(tdatafile);
                if(!open) {
                    outgather->open("o");
                    open = true;
                }else{
                    outgather->open("a");
                }
                outgather->writeTraces();
                outgather->close();
            }
        }
    }
    (*picks)[0]->Setchanged(false);
}

void MyApp::Load(wxCommandEvent& event)
{
    std::string filename;
    std::vector<Picks*> *picks;
    picks = frame->getPicks();
    int *npicks; 
    npicks = (*picks)[0]->Getnpicks();

    bool haspicks=false;
    for(int i=0; i< Sort->getNensemb(); i++)
    {
        if(npicks[i] > 0){
            haspicks=true;
        }
    }

    if (haspicks && (*picks)[0]->Getchanged())
    {
        if (wxMessageBox(_("Current picks have not been saved! Proceed?"), _("Please confirm"),
                         wxICON_QUESTION | wxYES_NO, frame) == wxNO )
            return;
        //else: proceed asking to the user the new file to open
    }

    wxFileDialog openFileDialog(frame, _("Open Picks file"), "", "", "Pick files (*.rss)|*.rss", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog.ShowModal() == wxID_CANCEL)
        return;     // the user changed idea...

    // proceed loading the file chosen by the user;
    filename = openFileDialog.GetPath();
    if((*picks)[0]->Loadpicks(filename) != PICKS_OK) {
        wxMessageBox(_("Something went wrong when reading from the picks file."), _("Error loading picks"), wxICON_INFORMATION);
        return;
    }
    frame->Refresh();
}

