#include "file.h"
#include "data.h"
#include "fft.h"
#include "utils.h"
#include "inparse.h"
#include <args.hxx>

using namespace rockseis;

int main(int argc, char* argv[])
{
   std::shared_ptr<rockseis::File> in (new rockseis::File());
   std::shared_ptr<rockseis::File> out (new rockseis::File());
   bool status;

    args::ArgumentParser parser("Program to bandpass filter data rss files.", "");
    parser.LongPrefix("");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<float> parflo(parser, "Double", "Minimum frequency", {"flo"});
    args::ValueFlag<float> parfhi(parser, "Double", "Maximum frequency", {"fhi"});
    args::ValueFlag<float> parf0(parser, "Double", "Frequency anomaly", {"f0"});
    args::ValueFlag<bool> parinv(parser, "Bool", "Inverse", {"inv"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cerr << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 0;
    }

    // Read input file
    status = in->input();
    if(status == FILE_ERR){
       rockseis::rs_error("Error reading from input file.");
    }

    float flo, fhi, f0;
    bool inv;

    if (parflo) { 
        flo = (float) args::get(parflo);
        if(flo < 0) rs_error("Negative flo");
        flo *= in->getD(1);
    }else{ 
        flo = 0.0;
    }
    if (parfhi) { 
        fhi = (float) args::get(parfhi);
        fhi *= in->getD(1);
    }else{
        fhi=0.5;
    }
    if(fhi/in->getD(1) < flo/in->getD(1)) rs_error ("Need flo < fhi");

    if (parinv) { 
        inv = (bool) args::get(parinv);
    }else{ 
        inv = false;
    }

   if (parf0) { 
        f0 = (float) args::get(parf0);
    }else{
        f0=10.0;
    }

    rockseis::rs_datatype type = static_cast<rockseis::rs_datatype>(in->getType());
    int esize_hdr = in->getHeader_format();
    int esize_data = in->getData_format();
    if(esize_hdr > 4 || esize_data > 4) rs_error("Double precision data not supported at the moment.");

    if((type == rockseis::DATA2D || type == rockseis::DATA3D))
    {

        std::shared_ptr<rockseis::Data2D<float>> Indata2d;
        std::shared_ptr<rockseis::Data3D<float>> Indata3d;
        std::shared_ptr<rockseis::Data2D<float>> Outdata2d;
        std::shared_ptr<rockseis::Data3D<float>> Outdata3d;
        std::shared_ptr<rockseis::Fft<float>> fft1d (new rockseis::Fft<float>(in->getN(1)));
        std::shared_ptr<rockseis::Fft<float>> ifft1d (new rockseis::Fft<float>(in->getN(1)));
        float *tracein;
        float *traceout;
        float *stran;
        float *amp;
        size_t ntr;

        size_t nflo, nfhi, nf;
        nflo = (size_t) ( flo*(fft1d->getNfft()) + 0.5 );
        nfhi = (size_t) ( fhi*(fft1d->getNfft()) + 0.5 );
        nf = nfhi-nflo + 1;
        float df = 1./(fft1d->getNfft()*in->getD(1));
        amp = (float *) calloc(in->getN(1)*nf, sizeof(float));
        stran = (float *) calloc(2*in->getN(1)*nf, sizeof(float));
        /*
        std::cerr << "main nf: " << nf << std::endl; 
        std::shared_ptr<File> Fout = std::make_shared<File>(); 
        Fout->output("stdout");
        Fout->setN(1,in->getN(1));
        Fout->setD(1,in->getD(1));
        Fout->setO(1,in->getO(1));
        Fout->setN(2,nf);
        Fout->setD(2,df);
        Fout->setO(2,flo/in->getD(1));
        Fout->setN(3,in->getN(2));
        Fout->setD(3,1);
        Fout->setNheader(0);
        Fout->setData_format(sizeof(float));
        Fout->setHeader_format(sizeof(float));
        Fout->writeHeader();
        Fout->seekp(Fout->getStartofdata());
        */


        // Create amplitude spectrum of Ricker window 
        float *ricker = (float *) calloc(nf, sizeof(float));
        float freq = 0.;
        for(int i=0; i<nf; i++){
            freq = flo/in->getD(1) + i*df;
            ricker[i] = ((2*freq*freq)/((sqrtf(PI))*(f0*f0*f0)))*expf(-(freq*freq/(f0*f0)));
        }

        float maxval = 0.0;
        for(int i=0; i<nf; i++){
            if(ricker[i] > maxval)
                maxval = ricker[i];
        }

        if(maxval > 0){
            for(int i=0; i<nf; i++){
                ricker[i] /= maxval;
            }
        }

        switch(type)
        {
            case DATA2D:
                Indata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Indata2d->setFdata(in);

                Outdata2d->setFile("stdout");
                status = Outdata2d->open("o");
                if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

                tracein = Indata2d->getData();
                traceout = Outdata2d->getData();
                ntr = in->getN(2);
                for(size_t i=0; i< ntr; i++)
                {
                    if(Indata2d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                    // Run S-transform
                    Indata2d->St1D(fft1d, ifft1d, tracein, nflo, nfhi, stran);
                    // Output amplitude spectrum
                    for(size_t j=0; j<in->getN(1)*nf; j++){
                       amp[j] = sqrtf(stran[2*j]*stran[2*j] + stran[2*j+1]*stran[2*j+1]);
                    }
                    for(size_t j=0; j<in->getN(1); j++){
                        traceout[j] = 0.0;
                        for (size_t k=0; k < nf; k++){
                            traceout[j] += amp[k*in->getN(1) + j]*ricker[k];
                            traceout[j] -= amp[k*in->getN(1) + j]*(1-ricker[k]);
                        }
                    }

                    for(size_t j=0; j<in->getN(1); j++){
                        if(traceout[j] < 0) traceout[j] = 0.0;
                    }
                    Outdata2d->copyCoords(Indata2d);
                    if(Outdata2d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                }
                Outdata2d->close();
                break;
            case DATA3D:
                Indata3d = std::make_shared<rockseis::Data3D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata3d = std::make_shared<rockseis::Data3D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Indata3d->setFdata(in);
                Outdata3d->copyCoords(Indata3d);

                Outdata3d->setFile("stdout");
                status = Outdata3d->open("o");
                if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

                tracein = Indata3d->getData();
                traceout = Outdata3d->getData();
                ntr = in->getN(2);
                for(size_t i=0; i< ntr; i++)
                {
                    if(Indata3d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                    // Run S-transform
                    Indata3d->St1D(fft1d, ifft1d, tracein, nflo, nfhi, stran);
                    // Output amplitude spectrum
                    for(size_t j=0; j<in->getN(1)*nf; j++){
                       amp[j] = sqrtf(stran[2*j]*stran[2*j] + stran[2*j+1]*stran[2*j+1]);
                    }
                    for(size_t j=0; j<in->getN(1); j++){
                        traceout[j] = 0.0;
                        for (size_t k=0; k < nf; k++){
                            traceout[j] += amp[k*in->getN(1) + j]*ricker[k];
                            traceout[j] -= amp[k*in->getN(1) + j]*(1-ricker[k]);
                        }
                    }

                    for(size_t j=0; j<in->getN(1); j++){
                        if(traceout[j] < 0) traceout[j] = 0.0;
                    }

                    Outdata3d->copyCoords(Indata3d);
                    if(Outdata3d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                }
                Outdata3d->close();
                break;
            default:
                break;
        }

        free(stran);
        free(amp);
        exit (0);
    }

    std::cerr << "Type not supported." << std::endl; 
    exit(1);

}
