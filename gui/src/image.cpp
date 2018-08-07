#include "image2dframe.h"
#include "file.h"
#include "utils.h"
#include "geometry.h"
#include <args.hxx>

// Class declaration of the app
class MyApp: public wxApp
{
public:
	virtual bool OnInit();
	int FilterEvent(wxEvent& event);
};

// Main implementation of app
wxIMPLEMENT_APP(MyApp);

wxDEFINE_EVENT(SelectCmp, wxCommandEvent);

// Class implementation of the app
bool MyApp::OnInit()
{
    std::shared_ptr<rockseis::File> in (new rockseis::File());
    bool status;
    status = in->input();
    if(status == FILE_ERR){
        rockseis::rs_error("rsrss2rsf: Error opening file for input.");
        std::cerr << status << std::endl;
    }

    size_t n1,n2,ntot;
    float d1,d2;
    float o1,o2;

    // Get parameters
    n1 = in->getN(1);
    n2 = in->getN(3);

    d1 = (float) in->getD(1);
    d2 = (float) in->getD(3);

    o1 = (float) in->getO(1);
    o2 = (float) in->getO(3);

    if(n2 == 0){
        n2 = in->getN(2);
        d2 = (float) in->getD(2);
        o2 = (float) in->getO(2);
    }

    ntot = in->getNtot();

    if(n1*n2 != ntot) rockseis::rs_error("Image data is not 2 dimensional");

    float *fdata=NULL;
    double *ddata=NULL;
    fdata = (float *) calloc(ntot, sizeof(float));
    if(fdata == NULL) rockseis::rs_error("Failed to allocate memory for the image.");

    int dsize, hsize;
    dsize = in->getData_format();

    int Nheader = in->getNheader();
    hsize = in->getHeader_format();

    char *hdata = NULL;
    if(Nheader){
        hdata = (char *) calloc(Nheader*hsize, sizeof(char));
    }
    switch(dsize){
        case 4:
            for (size_t j=0; j < n2; j++){
                in->read(&hdata[0], Nheader*hsize);
                in->read(&fdata[j*n1], n1);
            }
            break;
        case 8:
            ddata = (double *) calloc(ntot, sizeof(double));
            if(ddata == NULL) rockseis::rs_error("Failed to allocate memory for the image.");
            for (size_t j=0; j < n2; j++){
                if(Nheader){
                    in->read(&hdata[0], Nheader*hsize);
                }
                in->read(&ddata[j*n1], n1);
                for (size_t i=0; i < n1; i++){
                    fdata[j*n1 + i] = (float) ddata[j*n1 + i];
                }
            }
            free(ddata);
            break;
        default:
            rockseis::rs_error("Unsupported esize.");
            break;
    }

	wxString pick = "2D image";
	Image2dframe *frame = new Image2dframe(n1, d1, o1, n2, d2, o2, fdata, 0);
	frame->Show( true );

	return true;
}

int MyApp::FilterEvent(wxEvent& event)
{
	if (event.GetEventType() == wxEVT_KEY_UP)

	{
		event.Skip(true);
	}

	return -1;
}
