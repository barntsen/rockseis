#include "file.h"
#include "data.h"
#include "utils.h"
#include "inparse.h"
#include <args.hxx>

using namespace rockseis;

int main(int argc, char* argv[])
{
   std::shared_ptr<rockseis::File> in (new rockseis::File());
   std::shared_ptr<rockseis::File> ttimesfile (new rockseis::File());
   std::shared_ptr<rockseis::File> out (new rockseis::File());
   bool status;

   double v0;
   double delay;
   int ntapermax;
   int type;
   std::string ttimesfilename;

   args::ArgumentParser parser("Program to mute data rss files.", "");
   parser.LongPrefix("");
   parser.LongSeparator("=");
   args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
   args::ValueFlag<int> partype(parser, "Int", "type=0 mute below, type=1 mute above", {"type"});
   args::ValueFlag<double> parv0(parser, "Double", "Constant velocity used to compute direct wave traveltime (in m/s)", {"v0"});
   args::ValueFlag<double> pardelay(parser, "Double", "Time delay to add to the traveltime (in seconds)", {"delay"});
   args::ValueFlag<int> parntaper(parser, "Int", "Number of samples to use as a taper", {"ntaper"});
   args::ValueFlag<std::string> parttimes(parser, "FILENAME", "Optional file of same form as input data containing traveltimes (overrules v0)", {"ttimes"});
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

   if (parv0) { 
      v0 = (double) args::get(parv0);
   }else{ 
      if(!parttimes)
        rs_error("Must provide v0"); 
   }
   if (pardelay) { 
      delay = (double) args::get(pardelay);
   }else{
      delay = 0.0;
   }
   if (partype) { 
      type = (int) args::get(partype);
      if(type > 1 || type < 0 ) rs_error("type must be 0 or 1");
   }else{
      type = 0;
   }

   if(parntaper){
      ntapermax = (int) args::get(parntaper);
      if(ntapermax < 0) rs_error("ntaper must be positive");
   }else{
      ntapermax = 0;
   }

   if(parttimes){
      ttimesfilename = args::get(parttimes);
   }

   // Read input file
   status = in->input();
   if(status == FILE_ERR){
      rockseis::rs_error("Error reading from input file.");
   }

   if(parttimes){
      // Read traveltime file
      status = ttimesfile->input(ttimesfilename);
      if(status == FILE_ERR){
         rockseis::rs_error("Error reading from traveltime input file.");
      }
   }

   rockseis::rs_datatype rstype = static_cast<rockseis::rs_datatype>(in->getType());
   int esize_hdr = in->getHeader_format();
   int esize_data = in->getData_format();
   if(esize_hdr > 4 || esize_data > 4) rs_error("Double precision data not supported at the moment.");

   if((rstype == rockseis::DATA2D || rstype == rockseis::DATA3D))
   {

      std::shared_ptr<rockseis::Data2D<float>> Indata2d;
      std::shared_ptr<rockseis::Data3D<float>> Indata3d;
      std::shared_ptr<rockseis::Data2D<float>> Tdata2d;
      std::shared_ptr<rockseis::Data3D<float>> Tdata3d;
      std::shared_ptr<rockseis::Data2D<float>> Outdata2d;
      std::shared_ptr<rockseis::Data3D<float>> Outdata3d;

      rockseis::Point2D<float> *scoords2D = NULL;
      rockseis::Point2D<float> *gcoords2D = NULL;
      rockseis::Point3D<float> *scoords3D = NULL;
      rockseis::Point3D<float> *gcoords3D = NULL;
      float *tdata;
      double time;
      float *tracein;
      float *traceout;
      size_t ntr;
      int nd;
      int i0;
      int start, end, ntaper, itime;
      double taper;
      switch(rstype)
      {
         case DATA2D:
            Indata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
            Outdata2d = std::make_shared<rockseis::Data2D<float>>(1, in->getN(1), in->getD(1), in->getO(1));
            Indata2d->setFdata(in);
            if(parttimes){
               Tdata2d = std::make_shared<rockseis::Data2D<float>>(1, 1, 1.0, 0.0);
               Tdata2d->setFdata(ttimesfile);
               tdata = Tdata2d->getData();
            }

            Outdata2d->setFile("stdout");
            status = Outdata2d->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
            tracein = Indata2d->getData();
            traceout = Outdata2d->getData();
            ntr = in->getN(2);
            for(size_t i=0; i< ntr; i++)
            {
               if(Indata2d->readTraces() == FILE_ERR) rs_error("Error reading from input file");

               // Computing traveltimes
               if(parttimes){
                  if(Tdata2d->readTraces() == FILE_ERR) rs_error("Error reading from traveltimes input file");
                  time = tdata[0] + delay;
               }else{
                  scoords2D = (Indata2d->getGeom())->getScoords();
                  gcoords2D = (Indata2d->getGeom())->getGcoords();
                  time = std::sqrt((double) (SQ(scoords2D->x - gcoords2D->x) + SQ(scoords2D->y - gcoords2D->y)))/v0 + delay;
               }
               itime = (int) std::round((double) (time/Indata2d->getDt()));
               nd = (int) in->getN(1);
               for(int i0=0; i0<nd; i0++){
                  traceout[i0] = tracein[i0];
               }
               if(type){
                  start = itime - ntapermax;
                  end = itime;
                  if(start < 0) start = 0;
                  if(start > nd-1) start = nd-1;
                  ntaper = end - start;

                  for(i0=0; i0 < start; i0++){
                     traceout[i0]=0.0;
                  }
                  if(ntaper > 0){
                     for(i0=0; i0<ntaper; i0++){
                        taper = 0.5*(1.0-std::cos(PI*((double) i0)/(ntaper-1)));
                        taper=taper*taper;
                        if((i0+start)<nd){
                           traceout[i0+start] *= taper;
                        }
                     }
                  }
               }else{
                  start = itime;
                  end = itime + ntapermax;
                  if(start < 0) start = 0;
                  if(start > nd-1) start = nd-1;
                  if(end > nd-1) end = nd-1;
                  ntaper = end - start;

                  if(ntaper > 0){
                     for(i0=0; i0<ntaper; i0++){
                        taper = 0.5*(1.0-cosf(PI*((float) i0)/(ntaper-1)));
                        taper=taper*taper;
                        if((end - i0)>=0){
                           traceout[end-i0] *= taper;
                        }
                     }
                  }

                  for(i0=end+1; i0 < nd ; i0++){
                     traceout[i0]=0.0;
                  }
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
            if(parttimes){
               Tdata3d = std::make_shared<rockseis::Data3D<float>>(1, 1, 1.0, 0.0);
               Tdata3d->setFdata(ttimesfile);
               tdata = Tdata3d->getData();
            }

            Outdata3d->setFile("stdout");
            status = Outdata3d->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");
            tracein = Indata3d->getData();
            traceout = Outdata3d->getData();
            ntr = in->getN(2);
            for(size_t i=0; i< ntr; i++)
            {
               if(Indata3d->readTraces() == FILE_ERR) rs_error("Error reading from input file");

               // Computing traveltimes
               if(parttimes){
                  if(Tdata3d->readTraces() == FILE_ERR) rs_error("Error reading from traveltimes input file");
                  time = tdata[0] + delay;
               }else{
                  scoords3D = (Indata3d->getGeom())->getScoords();
                  gcoords3D = (Indata3d->getGeom())->getGcoords();
                  time = std::sqrt((double) (SQ(scoords3D->x - gcoords3D->x) + SQ(scoords3D->y - gcoords3D->y) + SQ(scoords3D->z - gcoords3D->z)))/v0 + delay;
               }
               itime = (int) std::round((double) (time/Indata3d->getDt()));
               nd = (int) in->getN(1);
               for(int i0=0; i0<nd; i0++){
                  traceout[i0] = tracein[i0];
               }
               if(type){
                  start = itime - ntapermax;
                  end = itime;
                  if(start < 0) start = 0;
                  if(start > nd-1) start = nd-1;
                  ntaper = end - start;

                  for(i0=0; i0 < start; i0++){
                     traceout[i0]=0.0;
                  }
                  if(ntaper > 0){
                     for(i0=0; i0<ntaper; i0++){
                        taper = 0.5*(1.0-std::cos(PI*((double) i0)/(ntaper-1)));
                        taper=taper*taper;
                        if((i0+start)<nd){
                           traceout[i0+start] *= taper;
                        }
                     }
                  }
               }else{
                  start = itime;
                  end = itime + ntapermax;
                  if(start < 0) start = 0;
                  if(start > nd-1) start = nd-1;
                  if(end > nd-1) end = nd-1;
                  ntaper = end - start;

                  if(ntaper > 0){
                     for(i0=0; i0<ntaper; i0++){
                        taper = 0.5*(1.0-cosf(PI*((float) i0)/(ntaper-1)));
                        taper=taper*taper;
                        if((end - i0)>=0){
                           traceout[end-i0] *= taper;
                        }
                     }
                  }

                  for(i0=end+1; i0 < nd ; i0++){
                     traceout[i0]=0.0;
                  }
               }

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
