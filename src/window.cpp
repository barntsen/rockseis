#include "file.h"
#include "data.h"
#include "utils.h"
#include "inparse.h"
#include <args.hxx>
#define Iin(i0,i1,i2,i3,i4,i5,i6,i7,i8) ((i8)*n_in[7]*n_in[6]*n_in[5]*n_in[4]*n_in[3]*n_in[2]*n_in[1]*n_in[0] + (i7)*n_in[6]*n_in[5]*n_in[4]*n_in[3]*n_in[2]*n_in[1]*n_in[0] + (i6)*n_in[5]*n_in[4]*n_in[3]*n_in[2]*n_in[1]*n_in[0] + (i5)*n_in[4]*n_in[3]*n_in[2]*n_in[1]*n_in[0] + (i4)*n_in[3]*n_in[2]*n_in[1]*n_in[0]  + (i3)*n_in[2]*n_in[1]*n_in[0] + (i2)*n_in[1]*n_in[0] + (i1)*n_in[0] + (i0))
#define Iout(i0,i1,i2,i3,i4,i5,i6,i7,i8) ((i8)*n[7]*n[6]*n[5]*n[4]*n[3]*n[2]*n[1]*n[0] + (i7)*n[6]*n[5]*n[4]*n[3]*n[2]*n[1]*n[0] + (i6)*n[5]*n[4]*n[3]*n[2]*n[1]*n[0] + (i5)*n[4]*n[3]*n[2]*n[1]*n[0] + (i4)*n[3]*n[2]*n[1]*n[0]  + (i3)*n[2]*n[1]*n[0] + (i2)*n[1]*n[0] + (i1)*n[0] + (i0))

using namespace rockseis;

int main(int argc, char* argv[])
{
	std::shared_ptr<rockseis::File> in (new rockseis::File());
	std::shared_ptr<rockseis::File> out (new rockseis::File());
    bool status;
    
    //Input varialbles
    long n[MAXDIMS], f[MAXDIMS], j[MAXDIMS];
    double min[MAXDIMS], max[MAXDIMS];
    
    // Cordinates
    std::string key;
    double minkey, maxkey;

    /* Parameters */
    size_t n_in[MAXDIMS];
    double d_in[MAXDIMS];
    double o_in[MAXDIMS];
    Point3D<float> s;
    Point3D<float> g;
    Point3D<float> offset;
    Point3D<float> aoffset;
    float foffset;

    args::ArgumentParser parser("Program to window rss data.", "Integer values f,n take precedence over double values in min, max.");
    parser.LongPrefix("");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<bool> parseto(parser, "Bool", "Do not set O's", {"notseto"});
    args::ValueFlag<double> parminkey(parser, "Double", "Mininum value in key", {"min"});
    args::ValueFlag<double> parmaxkey(parser, "Double", "Maximum value in key", {"max"});
    args::ValueFlag<std::string> parkey(parser, "String", "Key (s,g,foffset,aoffset) ex: sx, offsetx, foffset", {"key"});
    args::ValueFlag<double> *parmin[MAXDIMS];
    args::ValueFlag<double> *parmax[MAXDIMS];
    args::ValueFlag<long> *parn[MAXDIMS];
    args::ValueFlag<long> *parf[MAXDIMS];
    args::ValueFlag<long> *parj[MAXDIMS];
   
    std::string help_text;
    std::string key_text;
    for(int i=0; i<MAXDIMS; i++){
        help_text = "Minimum value in dim " + std::to_string(i+1);
        key_text = "min" + std::to_string(i+1);
        parmin[i] = new args::ValueFlag<double>(parser, "Double", help_text, {key_text});

        help_text = "Maximum value in dim " + std::to_string(i+1);
        key_text = "max" + std::to_string(i+1);
        parmax[i] = new args::ValueFlag<double>(parser, "Double", help_text, {key_text});

        help_text = "Number of elements to pass in dim " + std::to_string(i+1);
        key_text = "n" + std::to_string(i+1);
        parn[i] = new args::ValueFlag<long>(parser, "Integer", help_text, {key_text});

        help_text = "First element to pass in dim " + std::to_string(i+1);
        key_text = "f" + std::to_string(i+1);
        parf[i] = new args::ValueFlag<long>(parser, "Integer", help_text, {key_text});

        help_text = "How many elements to jump along dim " + std::to_string(i+1);
        key_text = "j" + std::to_string(i+1);
        parj[i] = new args::ValueFlag<long>(parser, "Integer", help_text, {key_text});
    }
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

    for(int i=0; i < MAXDIMS; i++)
    {
        n_in[i] = in->getN(i+1);
        d_in[i] = in->getD(i+1);
        o_in[i] = in->getO(i+1);
    }

    rockseis::rs_datatype type = static_cast<rockseis::rs_datatype>(in->getType());
    int esize_hdr = in->getHeader_format();
    int esize_data = in->getData_format();
    if(esize_hdr > 4 || esize_data > 4) rs_error("Double precision data not supported at the moment.");


    for(int i=0; i<MAXDIMS; i++){
        if(n_in[i]){
            // Get j's
            if (*parj[i]){
                j[i] = args::get(*parj[i]);
            }else{
                j[i] = 1;
            }

            // Get f's
            if (*parf[i]){
                f[i] = args::get(*parf[i]);
            }else{
                if(*parmin[i]){
                    min[i] = args::get(*parmin[i]);
                    f[i] = (long) (0.5 + (min[i] - o_in[i])/d_in[i]);
                }else{
                    f[i] = 0;
                }
            }
            if(f[i] < 0){
                rs_error("f[", std::to_string(i), "] is negative.");
            }

            //New values for d[i] and o[i]
            o_in[i] += f[i]*d_in[i];
            d_in[i] *= j[i];

            // Get n's
            if (*parn[i]){
                n[i] = args::get(*parn[i]);
            }else{
                if(*parmax[i]){
                    max[i] = args::get(*parmax[i]);
                    n[i] = (long) (1.5 + (max[i] - o_in[i])/d_in[i]);
                }else{
                    n[i] = (long) (1.5 + (n_in[i] - 1 - f[i])/j[i]);
                }
            }
            if(parseto){
                o_in[i] -= f[i]*d_in[i]/j[i];
            }
        }else{
            n[i] = 0;
            d_in[i] = 0.;
            o_in[i] = 0.;
            f[i] = 0;
            j[i] = 1;
        }
        if(n[i] > n_in[i]) rs_error("n[", std::to_string(i), "] is larger than input size.");

    }

    float *tracein = (float *) calloc(n_in[0], sizeof(float));
    float *traceout = (float *) calloc(n[0], sizeof(float));
    if((type == rockseis::GENERIC || type == rockseis::REGULAR || type == rockseis::SNAPSHOT))
    {
        out->output();
        for(long i=0; i<MAXDIMS; i++){
            if(n[i]){
                out->setN(i+1,n[i]);
                out->setD(i+1,d_in[i]);
                out->setO(i+1,o_in[i]);
            }
        }
        out->setType(type);
        out->setData_format(esize_data);
        out->setHeader_format(esize_hdr);
        out->writeHeader();

        for(int i=0; i<MAXDIMS; i++){
            if(n_in[i] == 0){
                n_in[i] = 1;
            }
            if(n[i] == 0){
                n[i] = 1;
            }
        }
        for(long i8=f[8]; i8<n_in[8]; i8+=j[8]){
        for(long i7=f[7]; i7<n_in[7]; i7+=j[7]){
        for(long i6=f[6]; i6<n_in[6]; i6+=j[6]){
        for(long i5=f[5]; i5<n_in[5]; i5+=j[5]){
        for(long i4=f[4]; i4<n_in[4]; i4+=j[4]){
        for(long i3=f[3]; i3<n_in[3]; i3+=j[3]){
        for(long i2=f[2]; i2<n_in[2]; i2+=j[2]){
        for(long i1=f[1]; i1<n_in[1]; i1+=j[1]){
            if((i1 -f[1])/j[1] > n[1]-1) break;
            if((i2 -f[2])/j[2] > n[2]-1) break;
            if((i3 -f[3])/j[3] > n[3]-1) break;
            if((i4 -f[4])/j[4] > n[4]-1) break;
            if((i5 -f[5])/j[5] > n[5]-1) break;
            if((i6 -f[6])/j[6] > n[6]-1) break;
            if((i7 -f[7])/j[7] > n[7]-1) break;
            if((i8 -f[8])/j[8] > n[8]-1) break;
            in->read(tracein, n_in[0], Iin(0,i1,i2,i3,i4,i5,i6,i7,i8)*esize_data);
            for(long i0=0; i0<n[0]; i0++){
                traceout[i0] = tracein[i0*j[0] + f[0]];
            }
            out->write(traceout, n[0], Iout(0,(i1-f[1])/j[1],(i2-f[2])/j[2],(i3-f[3])/j[3],(i4-f[4])/j[4],(i5-f[5])/j[5],(i6-f[6])/j[6],(i7-f[7])/j[7],(i8-f[8])/j[8])*esize_data);
        }
        }
        }
        }
        }
        }
        }
        }
        out->close();
        free(tracein);
        free(traceout);
        exit (0);
    }

    if((type == rockseis::DATA2D || type == rockseis::DATA3D))
    {

        std::shared_ptr<rockseis::Data2D<float>> Indata2d;
        std::shared_ptr<rockseis::Data3D<float>> Indata3d;
        std::shared_ptr<rockseis::Data2D<float>> Outdata2d;
        std::shared_ptr<rockseis::Data3D<float>> Outdata3d;
        if (parkey){
             key = args::get(parkey);
        }
        if (parminkey){
             minkey = args::get(parminkey);
        }
        if (parmaxkey){
             maxkey = args::get(parmaxkey);
        }
        if(parminkey && parmaxkey){
            if(minkey > maxkey){
                rs_error("Check Input: min > max.");
            }
        }
        rockseis::Point2D<float> *scoords2d;
        rockseis::Point2D<float> *gcoords2d;
        rockseis::Point3D<float> *scoords3d;
        rockseis::Point3D<float> *gcoords3d;
        float *tracein;
        float *traceout;

        size_t ntr;
        size_t count_traces = 0;
        switch(type)
        {
            case DATA2D:
                Indata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata2d = std::make_shared<rockseis::Data2D<float>>(1, n[0], d_in[0], o_in[0]);
                Indata2d->setFdata(in);
                scoords2d = (Indata2d->getGeom())->getScoords();
                gcoords2d = (Indata2d->getGeom())->getGcoords();

                Outdata2d->setFile("stdout");
                status = Outdata2d->open("o");
                if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
                tracein = Indata2d->getData();
                traceout = Outdata2d->getData();
                ntr = in->getN(2);
                for(size_t i=0; i< ntr; i++)
                {
                    if(Indata2d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                    if(i<f[1]) continue;
                    if((i-f[1])/j[1] > n[1]-1) break;
                    if((i-f[1]) % j[1]) continue;
                    s.x = scoords2d[0].x;
                    g.x = gcoords2d[0].x;
                    s.y = scoords2d[0].y;
                    g.y = gcoords2d[0].y;
                    offset.x = gcoords2d[0].x - scoords2d[0].x;
                    offset.y = gcoords2d[0].y - scoords2d[0].y;
                    aoffset.x = fabsf(gcoords2d[0].x - scoords2d[0].x);
                    aoffset.y = fabsf(gcoords2d[0].y - scoords2d[0].y);
                    foffset = sqrtf((aoffset.x*aoffset.x) + (aoffset.y*aoffset.y)); 
                    if(strncmp(key.c_str(), "sx", 2) == 0) {
                        if(parminkey && s.x < minkey) continue; 
                        if(parmaxkey && s.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "gx", 2) == 0) {
                        if(parminkey && g.x < minkey) continue; 
                        if(parmaxkey && g.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "sz", 2) == 0) {
                        if(parminkey && s.y < minkey) continue; 
                        if(parmaxkey && s.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "gz", 2) == 0) {
                        if(parminkey && g.y < minkey) continue; 
                        if(parmaxkey && g.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "offsetx", 7) == 0) {
                        if(parminkey && offset.x < minkey) continue; 
                        if(parmaxkey && offset.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "offsetz", 7) == 0) {
                        if(parminkey && offset.y < minkey) continue; 
                        if(parmaxkey && offset.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "aoffsetx", 8) == 0) {
                        if(parminkey && aoffset.x < minkey) continue; 
                        if(parmaxkey && aoffset.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "aoffsetz", 8) == 0) {
                        if(parminkey && aoffset.y < minkey) continue; 
                        if(parmaxkey && aoffset.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "foffset", 6) == 0) {
                        if(parminkey && foffset < minkey) continue; 
                        if(parmaxkey && foffset > maxkey) continue; 
                    }

                    for(long i0=0; i0<n[0]; i0++){
                        traceout[i0] = tracein[i0*j[0] + f[0]];
                    }
                    Outdata2d->copyCoords(Indata2d);
                    if(Outdata2d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                    count_traces++;
                }
                Outdata2d->close();
                if(count_traces == 0){
                    rs_error("Output has zero traces!");
                }
                break;
            case DATA3D:
                Indata3d = std::make_shared<rockseis::Data3D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
                Outdata3d = std::make_shared<rockseis::Data3D<float>>(1, n[0], d_in[0], o_in[0]);
                Indata3d->setFdata(in);
                scoords3d = (Indata3d->getGeom())->getScoords();
                gcoords3d = (Indata3d->getGeom())->getGcoords();

                Outdata3d->setFile("stdout");
                status = Outdata3d->open("o");
                if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

                tracein = Indata3d->getData();
                traceout = Outdata3d->getData();
                ntr = in->getN(2);
                for(size_t i=0; i< ntr; i++)
                {
                    if(Indata3d->readTraces() == FILE_ERR) rs_error("Error reading from input file");
                    if(i<f[1]) continue;
                    if((i-f[1])/j[1] > n[1]-1) break;
                    if((i-f[1]) % j[1]) continue;
                    s.x = scoords3d[0].x;
                    g.x = gcoords3d[0].x;
                    s.y = scoords3d[0].y;
                    g.y = gcoords3d[0].y;
                    s.z = scoords3d[0].z;
                    g.z = gcoords3d[0].z;
                    offset.x = gcoords3d[0].x - scoords3d[0].x;
                    offset.y = gcoords3d[0].y - scoords3d[0].y;
                    offset.z = gcoords3d[0].z - scoords3d[0].z;

                    aoffset.x = fabsf(gcoords3d[0].x - scoords3d[0].x);
                    aoffset.y = fabsf(gcoords3d[0].y - scoords3d[0].y);
                    aoffset.z = fabsf(gcoords3d[0].z - scoords3d[0].z);

                    foffset = sqrtf((aoffset.x*aoffset.x) + (aoffset.y*aoffset.y) + (aoffset.z*aoffset.z)); 

                    if(strncmp(key.c_str(), "sx", 2) == 0) {
                        if(parminkey && s.x < minkey) continue; 
                        if(parmaxkey && s.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "gx", 2) == 0) {
                        if(parminkey && g.x < minkey) continue; 
                        if(parmaxkey && g.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "sy", 2) == 0) {
                        if(parminkey && s.y < minkey) continue; 
                        if(parmaxkey && s.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "gy", 2) == 0) {
                        if(parminkey && g.y < minkey) continue; 
                        if(parmaxkey && g.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "sz", 2) == 0) {
                        if(parminkey && s.z < minkey) continue; 
                        if(parmaxkey && s.z > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "gz", 2) == 0) {
                        if(parminkey && g.z < minkey) continue; 
                        if(parmaxkey && g.z > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "offsetx", 7) == 0) {
                        if(parminkey && offset.x < minkey) continue; 
                        if(parmaxkey && offset.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "offsety", 7) == 0) {
                        if(parminkey && offset.y < minkey) continue; 
                        if(parmaxkey && offset.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "offsetz", 7) == 0) {
                        if(parminkey && offset.z < minkey) continue; 
                        if(parmaxkey && offset.z > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "aoffsetx", 8) == 0) {
                        if(parminkey && aoffset.x < minkey) continue; 
                        if(parmaxkey && aoffset.x > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "aoffsety", 8) == 0) {
                        if(parminkey && aoffset.y < minkey) continue; 
                        if(parmaxkey && aoffset.y > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "aoffsetz", 8) == 0) {
                        if(parminkey && aoffset.z < minkey) continue; 
                        if(parmaxkey && aoffset.z > maxkey) continue; 
                    }
                    if(strncmp(key.c_str(), "foffset", 6) == 0) {
                        if(parminkey && foffset < minkey) continue; 
                        if(parmaxkey && foffset > maxkey) continue; 
                    }

                    for(long i0=0; i0<n[0]; i0++){
                        traceout[i0] = tracein[i0*j[0] + f[0]];
                    }
                    Outdata3d->copyCoords(Indata3d);
                    if(Outdata3d->writeTraces() == FILE_ERR) rs_error("Error writting to output file");
                    count_traces++;
                }
                Outdata3d->close();
                if(count_traces == 0){
                    rs_error("Output has zero traces!");
                }
                break;
            default:
                break;
        }

        exit (0);
    }

    std::cerr << "Type not supported." << std::endl; 
    exit(1);

}
