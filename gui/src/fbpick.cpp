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
	int FilterEvent(wxEvent& event);
    void transposeData();

private:
    Image2dframe *frame;
    std::shared_ptr<rockseis::Sort<float>> Sort;
    std::shared_ptr<rockseis::Data2D<float>> gather;

};

// Main implementation of app
wxIMPLEMENT_APP(MyApp);

wxDEFINE_EVENT(SelectCmp, wxCommandEvent);

// Class implementation of the app
bool MyApp::OnInit()
{
    args::ArgumentParser parser("Program to pick first break on DATA type files.", "");
    parser.LongPrefix("-");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> input(parser, "FILENAME", "Input file", {"if"});
    //args::ValueFlag<std::string> output(parser, "FILENAME", "Output file", {"of"});
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
    std::string infile;
    if (input){
         infile = args::get(input);
         std::cerr << "-if = " << infile << std::endl;
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing input file -if");
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
