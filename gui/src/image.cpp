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
wxDEFINE_EVENT(SaveEvent, wxCommandEvent);
wxDEFINE_EVENT(LoadEvent, wxCommandEvent);
wxDEFINE_EVENT(Crosshair, wxCommandEvent);
wxDEFINE_EVENT(MuteEvent, wxCommandEvent);

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

    size_t n[2],ntot;
    float d[2];
    float o[2];

    int ndims=0;
    // Get parameters
    for(int i=0; i<MAXDIMS; i++){
        if(in->getN(i+1) > 1){
            if(ndims < 2){
                n[ndims]  = in->getN(i+1);
                d[ndims]  = (float) in->getD(i+1);
                o[ndims]  = (float) in->getO(i+1);
            }
            ndims ++;
        }
    }
    if(ndims > 2) rockseis::rs_error("Image data is not 2 dimensional");
    ntot = n[0]*n[1];


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
            for (size_t j=0; j < n[1]; j++){
                in->read(&hdata[0], Nheader*hsize);
                in->read(&fdata[j*n[0]], n[0]);
            }
            break;
        case 8:
            ddata = (double *) calloc(ntot, sizeof(double));
            if(ddata == NULL) rockseis::rs_error("Failed to allocate memory for the image.");
            for (size_t j=0; j < n[1]; j++){
                if(Nheader){
                    in->read(&hdata[0], Nheader*hsize);
                }
                in->read(&ddata[j*n[0]], n[0]);
                for (size_t i=0; i < n[0]; i++){
                    fdata[j*n[0] + i] = (float) ddata[j*n[0] + i];
                }
            }
            free(ddata);
            break;
        default:
            rockseis::rs_error("Unsupported esize.");
            break;
    }

	wxString pick = "2D image";
	Image2dframe *frame = new Image2dframe(n[0], d[0], o[0], n[1], d[1], o[1], fdata, 0);
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
