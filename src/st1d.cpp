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

    float flo, fhi;

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
        std::shared_ptr<rockseis::Fft<float>> fft1d (new rockseis::Fft<T>(in->getN(1)));
        float *tracein;
        float *traceout;
        size_t ntr;

        size_t nflo, nfhi, nf;
        nflo = (size_t) ( flo*(fft1d->getNfft()/2 + 1) + 0.5 )
        nfhi = (size_t) ( fhi*(fft1d->getNfft()/2 + 1) + 0.5 )
        nf = nfhi-nflo + 1;
        switch(type)
        {
            case DATA2D:
                Indata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata2d = std::make_shared<rockseis::Data2D<float>>(nf, in->getN(1), in->getD(1), in->getO(1));
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
                    for(long i0=0; i0<ipad; i0++){
                        traceout[i0] = 0.0;
                    }
                    for(long i0=0; i0<in->getN(1); i0++){
                        traceout[i0] = tracein[i0];
                    }
                    Outdata2d->apply_filter(&freqs[0]);
                    Outdata2d->copyCoords(Indata2d);
                    if(Outdata2d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                }
                Outdata2d->close();
                break;
            case DATA3D:
                Indata3d = std::make_shared<rockseis::Data3D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata3d = std::make_shared<rockseis::Data3D<float>>(1, (ipad + in->getN(1)), in->getD(1), in->getO(1));
                Indata3d->setFdata(in);

                Outdata3d->setFile("stdout");
                status = Outdata3d->open("o");
                if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
                tracein = Indata3d->getData();
                traceout = Outdata3d->getData();
                ntr = in->getN(2);
                for(size_t i=0; i< ntr; i++)
                {
                    if(Indata3d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                    for(long i0=0; i0<ipad; i0++){
                        traceout[i0] = 0.0;
                    }
                    for(long i0=0; i0<in->getN(1); i0++){
                        traceout[i0+ipad] = tracein[i0];
                    }
                    
                    Outdata3d->apply_filter(&freqs[0]);
                    Outdata3d->copyCoords(Indata3d);
                    if(Outdata3d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                }
                Outdata3d->close();

                break;
            default:
                break;
        }

        exit (0);
    }

    std::cerr << "Type not supported." << std::endl; 
    exit(1);

}
