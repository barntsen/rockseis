#include "image2dframe.h"
#include "file.h"
#include "sort.h"
#include "utils.h"
#include "geometry.h"
#include <args.hxx>

// Class declaration of the app
class MyApp: public wxApp
{
public:
	virtual bool OnInit();
	void readGath(wxCommandEvent& event);
	void Save(wxCommandEvent& event);
	int FilterEvent(wxEvent& event);
    void transposeData();

private:
    Image2dframe *frame;
    std::shared_ptr<rockseis::Sort<float>> Sort;
    std::shared_ptr<rockseis::Data2D<float>> gather;
    std::string infile;
    std::string picksfile;

};

// Main implementation of app
wxIMPLEMENT_APP(MyApp);

wxDEFINE_EVENT(SelectCmp, wxCommandEvent);
wxDEFINE_EVENT(SaveEvent, wxCommandEvent);
wxDEFINE_EVENT(LoadEvent, wxCommandEvent);

// Class implementation of the app
bool MyApp::OnInit()
{
    args::ArgumentParser parser("Program to pick first break on DATA type files.", "");
    parser.LongPrefix("-");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> input(parser, "FILENAME", "Input file", {"i"});
    args::ValueFlag<std::string> output(parser, "FILENAME", "Output picks file", {"o"});
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
    if (input){
         infile = args::get(input);
         std::cerr << "-i  " << infile << std::endl;
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing input file -i");
    }
    if (output){
         picksfile = args::get(output);
         std::cerr << "-o  " << picksfile << std::endl;
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing input file -o");
    }


    Sort = std::make_shared<rockseis::Sort<float>>();
    Sort->setDatafile(infile);
    Sort->createShotmap(infile); 

    gather = Sort->get2DGather(0);
    this->transposeData();

    // Get parameters
	wxString pick = "Gather";
    frame = new Image2dframe(gather->getNtrace(), 1, 0, gather->getNt(), gather->getDt(), gather->getOt(), gather->getData(), 0);
    frame->createToolbar();
    frame->SetLabel (wxT("Picking window"));
    frame->setMaxcmp(Sort->getNensemb()-1);
    frame->createPicks();
    frame->Show( true );

	Bind(SelectCmp, &MyApp::readGath, this, frame->GetId());
	Bind(SaveEvent, &MyApp::Save, this, frame->GetId());
	return true;
}

void MyApp::readGath(wxCommandEvent& event)
{
    int ix = event.GetInt();
    if(ix < 0) ix = 0;
    if(ix > Sort->getNensemb()-1) ix = Sort->getNensemb()-1;
    gather.reset();
    gather = Sort->get2DGather(ix);
    if (gather == NULL){
        rockseis::rs_error("Error getting gather.");
        exit(1);
    }
    this->transposeData();

    frame->setN1(gather->getNtrace());
    frame->setN2(gather->getNt());
    frame->setImagedata(gather->getData());
    frame->ComputeClip();
    frame->setZoom(frame->getO1(), (frame->getN1()-1)*frame->getD1() + frame->getO1(), frame->getO2(), (frame->getN2()-1)*frame->getD2() + frame->getO2(), 0, frame->getN1(), 0, frame->getN2());
    frame->LoadImage(0, gather->getNtrace(), 0, gather->getNt());
    frame->setStatusbar(ix, 0., 0., 0.); 
    frame->Refresh();
}

void MyApp::transposeData()
{
    int nt = gather->getNt();
    int ntr = gather->getNtrace();
    rockseis::Index Idata(nt,ntr);
    rockseis::Index Itransp(ntr,nt);
    float *fdata = gather->getData();
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
    std::shared_ptr<rockseis::Data2D<float>> wrkgather;
    std::shared_ptr<rockseis::Data2D<float>> outgather;
    std::vector<Picks*> *picks;
    picks = frame->getPicks();
    int *npicks; 
    npicks = (*picks)[0]->Getnpicks();
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
        wrkgather = Sort->get2DGather(0);
        outgather = std::make_shared<rockseis::Data2D<float>>(wrkgather->getNtrace(), 1, 1.0, 0.0);
        outgather->setFile(picksfile);
        outgather->open("o");
        outgather->close();
        for(int i=0; i< Sort->getNensemb(); i++)
        {
            if(npicks[i] > 0){
                std::cerr << "Starting to save picks" << std::endl;
                //Get gather
                wrkgather.reset();
                wrkgather = Sort->get2DGather(i);
                outgather.reset();
                outgather = std::make_shared<rockseis::Data2D<float>>(wrkgather->getNtrace(), 1, 1.0, 0.0);
                data = outgather->getData();
                //Compute interpolated picks
                (*picks)[0]->Interp(i);
                vrms  = (*picks)[0]->Getvrms();
                for(int j=0; j<gather->getNtrace(); j++){
                    data[j] = vrms[j];
                }
                outgather->setFile(picksfile);
                outgather->open("a");
                outgather->writeTraces();
            }
        }
        wrkgather->close();
    }
}
