#include <iostream>
#include <stdio.h>
#include <memory>
#include <fstream>
#include <math.h>
#include <file.h>
#include "data.h"
#include "sort.h"
#include "utils.h"
#include "geometry.h"
#include <args.hxx>

int main(int argc, char* argv[])
{

    /* Variables */
    std::string Cproj;
    std::string Surveyname;
    FILE *fp;

    args::ArgumentParser parser("Program to convert RSS data file to SPS survey files.", "");
    parser.LongPrefix("");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> parproj(parser, "ex. UTM32N", "Projection", {"proj"});
    args::ValueFlag<std::string> parname(parser, "ex. Line1_Hafrsfjord", "Survey name", {"name"});
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cerr << parser;
        return 1;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (parproj){
         Cproj = args::get(parproj);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing projection.");
    }

    if (parname){
         Surveyname = args::get(parname);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing survey name.");
    }

    //Creating SPS survey files

    char sfilename[256];
    char rfilename[256];
    char xfilename[256];
    snprintf(sfilename, 256, "%s_SPS_Sources.s", Surveyname.c_str());
    snprintf(rfilename, 256, "%s_SPS_Receivers.r", Surveyname.c_str());
    snprintf(xfilename, 256, "%s_SPS_Templates.x", Surveyname.c_str());

    // Writting Source file
    fp = fopen(sfilename, "w");
    fprintf(fp,"H00 SPS format version num.     SPS 2.1;\n");
    fprintf(fp,"H01 Survey Description          rss2sps;\n");
    fprintf(fp,"H02 Date of Survey              date;predicted ;\n");
    fprintf(fp,"H03 Client                      University of Stavanger;\n");
    fprintf(fp,"H04 Geophysical Contractor      University of Stavanger;\n");
    fprintf(fp,"H05 Positioning Contractor      Wiktor Weibull;\n");
    fprintf(fp,"H06 Post Proc. Contractor       Wiktor Weibull;\n");
    fprintf(fp,"H07 Field Computer System(s)    RT System 2;\n");
    fprintf(fp,"H08 Coordination location       N/A;\n");
    fprintf(fp,"H09 Offset from coord. location N/A;\n");
    fprintf(fp,"H11 Spare                       N/A;\n");
    fprintf(fp,"H12 Geodetic datum,-spheroid    WGS84 , WGS84 ;\n");
    fprintf(fp,"H13 Spare                       N/A;\n");
    fprintf(fp,"H14 Geodetic datum parameters   0,0,0;\n");
    fprintf(fp,"H26 H12 is WGS datum transformation parameters ;\n");
    fprintf(fp,"H15 Spare                       N/A;\n");
    fprintf(fp,"H16 Spare                       N/A;\n");
    fprintf(fp,"H17 Vertical datum description  MSL - mean sea level;\n");
    fprintf(fp,"H18 Projection type             utm;\n");
    fprintf(fp,"H19 Projection zone             32 ;\n");
    fprintf(fp,"H20 Description of grid units   m ;\n");
    fprintf(fp,"H260Coord Sys Friendly Name     WGS 84 / UTM zone 32N(epsg:32632)(M);\n");
    fprintf(fp,"H26 Receiver Type               R1,Three Component ;\n");
    fprintf(fp,"H600Type, Model, Polarity       R1,Three Component,OTHER,Seg ;\n");
    fprintf(fp,"H601Damp coeff, natural freq.   R1,70,10 ;\n");
    fprintf(fp,"H603Unit spacing X, Y           R1,10,10 ;\n");
    fprintf(fp,"H26 I/O ITN of Instrument Type  1 ;\n");
    fprintf(fp,"H26 Proc Flow,Acq. Module Name  1,Processing Flow,Acquisition ;\n");
    fprintf(fp,"H400Type, Model, Polarity       1,N/A,N/A,Normal ;\n");
    fprintf(fp,"H401Crew Name,Comment           1,This Crew,Climate or User Observation ;\n");
    fprintf(fp,"H402Sample int.,Record Len.     1,2,1 ;\n");
    fprintf(fp,"H403Number of Channels          1,N/A ;\n");
    fprintf(fp,"H404Tape type, format, density  1,N/A,SEG-Y WSI,N/A ;\n");
    fprintf(fp,"H405Filter_alias Hz,dB pnt,slope1,400,-3,315 ;\n");
    fprintf(fp,"H406Filter_notch Hz,-3dB points 1,N/A,N/A ;\n");
    fprintf(fp,"H407Filter_low Hz,dB pnt,slope  1,3,-3,6 ;\n");
    fprintf(fp,"H408Time delay FTB-SOD app Y/N  1,0 ;\n");
    fprintf(fp,"H409Multi component recording   1,N/A ;\n");
    fprintf(fp,"H410Aux. channel 1 contents     1,NONE ;\n");
    fprintf(fp,"H411Aux. channel 2 contents     1,NONE ;\n");
    fprintf(fp,"H412Aux. channel 3 contents     1,NONE ;\n");
    fprintf(fp,"H413Aux. channel 4 contents     1,NONE ;\n");
    fprintf(fp,"H26 ............................................................................;\n");
    fprintf(fp,"H26                                                                             ;\n");
    fprintf(fp,"H26      Sp     Sp        Static          Water                Surface  Dayof   ;\n");
    fprintf(fp,"H26      Line   Station   Correction      Depth               Elevation Year    ;\n");
    fprintf(fp,"H26       |         |BLANK   |               |                        |  |   Sec;\n");
    fprintf(fp,"H26       |         | |Sp    | Pt   Uphole   |                        |  |     |;\n");
    fprintf(fp,"H26       |         | |Index |Depth  Time    |           Northing     |  | Min |;\n");
    fprintf(fp,"H26       |         | ||     |   |     |     |                  |     |  |   | |;\n");
    fprintf(fp,"H26       |         | ||Sp   |   |Seis |     |                  |     |  |Hr | |;\n");
    fprintf(fp,"H26       |         | ||Code |   |Datum|     |  Easting         |     |  | | | |;\n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |;\n");
    fprintf(fp,"H26      1         2         3         4         5         6         7         8;\n");
    fprintf(fp,"H26 5678901234567890123456789012345678901234567890123456789012345678901234567890;\n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |;\n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |;\n");
    fclose(fp);

    std::string filename = "stdin";
	std::shared_ptr<rockseis::Data3D<float>> OneShot;
	std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
    Sort->createShotmap(filename); 
    int ns = Sort->getNensemb();
    std::cerr << "Number of shot gathers: " << ns << std::endl;
    rockseis::Point3D<float> *scoords;
    fp = fopen(sfilename, "a");
    float Spstation = 1001.0;
    while((OneShot = Sort->get3DGather()) != nullptr ){
        scoords = (OneShot->getGeom())->getScoords();
        fprintf(fp, "S%10.1f%10.1f %2dS1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%6d\n", 1001.0,Spstation,1,0,0.0,0,0,0.0,scoords[0].x,scoords[0].y,scoords[0].z,0,0);
        Spstation += 1.0;
        OneShot.reset();
    }
    fclose(fp);


    // Writting Receiver file
    fp = fopen(rfilename, "w");
    fprintf(fp,"H00 SPS format version num.     SPS 2.1;\n");
    fprintf(fp,"H01 Survey Description          rss2sps;\n");
    fprintf(fp,"H02 Date of Survey              date;predicted ;\n");
    fprintf(fp,"H03 Client                      University of Stavanger;\n");
    fprintf(fp,"H04 Geophysical Contractor      University of Stavanger;\n");
    fprintf(fp,"H05 Positioning Contractor      Wiktor Weibull;\n");
    fprintf(fp,"H06 Post Proc. Contractor       Wiktor Weibull;\n");
    fprintf(fp,"H07 Field Computer System(s)    RT System 2;\n");
    fprintf(fp,"H08 Coordination location       N/A;\n");
    fprintf(fp,"H09 Offset from coord. location N/A;\n");
    fprintf(fp,"H11 Spare                       N/A;\n");
    fprintf(fp,"H12 Geodetic datum,-spheroid    WGS84 , WGS84 ;\n");
    fprintf(fp,"H13 Spare                       N/A;\n");
    fprintf(fp,"H14 Geodetic datum parameters   0,0,0;\n");
    fprintf(fp,"H26 H12 is WGS datum transformation parameters ;\n");
    fprintf(fp,"H15 Spare                       N/A;\n");
    fprintf(fp,"H16 Spare                       N/A;\n");
    fprintf(fp,"H17 Vertical datum description  MSL - mean sea level;\n");
    fprintf(fp,"H18 Projection type             utm;\n");
    fprintf(fp,"H19 Projection zone             32 ;\n");
    fprintf(fp,"H20 Description of grid units   m ;\n");
    fprintf(fp,"H260Coord Sys Friendly Name     WGS 84 / UTM zone 32N(epsg:32632)(M);\n");
    fprintf(fp,"H26 Receiver Type               R1,Three Component ;\n");
    fprintf(fp,"H600Type, Model, Polarity       R1,Three Component,OTHER,Seg ;\n");
    fprintf(fp,"H601Damp coeff, natural freq.   R1,70,10 ;\n");
    fprintf(fp,"H603Unit spacing X, Y           R1,10,10 ;\n");
    fprintf(fp,"H26 I/O ITN of Instrument Type  1 ;\n");
    fprintf(fp,"H26 Proc Flow,Acq. Module Name  1,Processing Flow,Acquisition ;\n");
    fprintf(fp,"H400Type, Model, Polarity       1,N/A,N/A,Normal ;\n");
    fprintf(fp,"H401Crew Name,Comment           1,This Crew,Climate or User Observation ;\n");
    fprintf(fp,"H402Sample int.,Record Len.     1,2,1 ;\n");
    fprintf(fp,"H403Number of Channels          1,N/A ;\n");
    fprintf(fp,"H404Tape type, format, density  1,N/A,SEG-Y WSI,N/A ;\n");
    fprintf(fp,"H405Filter_alias Hz,dB pnt,slope1,400,-3,315 ;\n");
    fprintf(fp,"H406Filter_notch Hz,-3dB points 1,N/A,N/A ;\n");
    fprintf(fp,"H407Filter_low Hz,dB pnt,slope  1,3,-3,6 ;\n");
    fprintf(fp,"H408Time delay FTB-SOD app Y/N  1,0 ;\n");
    fprintf(fp,"H409Multi component recording   1,N/A ;\n");
    fprintf(fp,"H410Aux. channel 1 contents     1,NONE ;\n");
    fprintf(fp,"H411Aux. channel 2 contents     1,NONE ;\n");
    fprintf(fp,"H412Aux. channel 3 contents     1,NONE ;\n");
    fprintf(fp,"H413Aux. channel 4 contents     1,NONE ;\n");
    fprintf(fp,"H26 ............................................................................; \n");
    fprintf(fp,"H26                                                                             ; \n");
    fprintf(fp,"H26       RCV     RCV     Static          Water                Surface  Dayof   ; \n");
    fprintf(fp,"H26       Line   Station  Correction      Depth               Elevation Year    ; \n");
    fprintf(fp,"H26       |         |BLANK   |               |                        |  |   Sec; \n");
    fprintf(fp,"H26       |         | |RCV   | Pt   Uphole   |                        |  |     |; \n");
    fprintf(fp,"H26       |         | |Index |Depth  Time    |           Northing     |  | Min |; \n");
    fprintf(fp,"H26       |         | ||     |   |     |     |                  |     |  |   | |; \n");
    fprintf(fp,"H26       |         | ||RCV  |   |Seis |     |                  |     |  |Hr | |; \n");
    fprintf(fp,"H26       |         | ||Code |   |Datum|     |  Easting         |     |  | | | |; \n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |; \n");
    fprintf(fp,"H26      1         2         3         4         5         6         7         8; \n");
    fprintf(fp,"H26 5678901234567890123456789012345678901234567890123456789012345678901234567890; \n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |; \n");
    fprintf(fp,"H26       |         | || |   |   |   | |     |        |         |     |  | | | |; \n");
    fclose(fp);

    Sort->createReceivermap(filename); 
    int nr = Sort->getNensemb();
    std::cerr << "Number of receiver gathers: " << nr << std::endl;
    rockseis::Point3D<float> *gcoords;
    Spstation = 1001.0;
    fp = fopen(rfilename, "a");
    while((OneShot = Sort->get3DGather()) != nullptr ){
        gcoords = (OneShot->getGeom())->getGcoords();
        fprintf(fp, "R%10.1f%10.1f %2dR1%4d%4.1f%4d%2d%6.1f%9.1f%10.1f%6.1f%3d%6d\n", 1001.0,Spstation,1,0,0.0,0,0,0.0,gcoords[0].x,gcoords[0].y,gcoords[0].z,0,0);
        Spstation += 1.0;
        OneShot.reset();
    }
    fclose(fp);

    // Writting Template file
    fp = fopen(xfilename, "w");
    fprintf(fp,"H00 SPS format version num.     SPS 2.1;\n");
    fprintf(fp,"H01 Survey Description          rss2sps;\n");
    fprintf(fp,"H02 Date of Survey              date;predicted ;\n");
    fprintf(fp,"H03 Client                      University of Stavanger;\n");
    fprintf(fp,"H04 Geophysical Contractor      University of Stavanger;\n");
    fprintf(fp,"H05 Positioning Contractor      Wiktor Weibull;\n");
    fprintf(fp,"H06 Post Proc. Contractor       Wiktor Weibull;\n");
    fprintf(fp,"H07 Field Computer System(s)    RT System 2;\n");
    fprintf(fp,"H08 Coordination location       N/A;\n");
    fprintf(fp,"H09 Offset from coord. location N/A;\n");
    fprintf(fp,"H11 Spare                       N/A;\n");
    fprintf(fp,"H12 Geodetic datum,-spheroid    WGS84 , WGS84 ;\n");
    fprintf(fp,"H13 Spare                       N/A;\n");
    fprintf(fp,"H14 Geodetic datum parameters   0,0,0;\n");
    fprintf(fp,"H26 H12 is WGS datum transformation parameters ;\n");
    fprintf(fp,"H15 Spare                       N/A;\n");
    fprintf(fp,"H16 Spare                       N/A;\n");
    fprintf(fp,"H17 Vertical datum description  MSL - mean sea level;\n");
    fprintf(fp,"H18 Projection type             utm;\n");
    fprintf(fp,"H19 Projection zone             32 ;\n");
    fprintf(fp,"H20 Description of grid units   m ;\n");
    fprintf(fp,"H260Coord Sys Friendly Name     WGS 84 / UTM zone 32N(epsg:32632)(M);\n");
    fprintf(fp,"H26 Receiver Type               R1,Three Component ;\n");
    fprintf(fp,"H600Type, Model, Polarity       R1,Three Component,OTHER,Seg ;\n");
    fprintf(fp,"H601Damp coeff, natural freq.   R1,70,10 ;\n");
    fprintf(fp,"H603Unit spacing X, Y           R1,10,10 ;\n");
    fprintf(fp,"H26 I/O ITN of Instrument Type  1 ;\n");
    fprintf(fp,"H26 Proc Flow,Acq. Module Name  1,Processing Flow,Acquisition ;\n");
    fprintf(fp,"H400Type, Model, Polarity       1,N/A,N/A,Normal ;\n");
    fprintf(fp,"H401Crew Name,Comment           1,This Crew,Climate or User Observation ;\n");
    fprintf(fp,"H402Sample int.,Record Len.     1,2,1 ;\n");
    fprintf(fp,"H403Number of Channels          1,N/A ;\n");
    fprintf(fp,"H404Tape type, format, density  1,N/A,SEG-Y WSI,N/A ;\n");
    fprintf(fp,"H405Filter_alias Hz,dB pnt,slope1,400,-3,315 ;\n");
    fprintf(fp,"H406Filter_notch Hz,-3dB points 1,N/A,N/A ;\n");
    fprintf(fp,"H407Filter_low Hz,dB pnt,slope  1,3,-3,6 ;\n");
    fprintf(fp,"H408Time delay FTB-SOD app Y/N  1,0 ;\n");
    fprintf(fp,"H409Multi component recording   1,N/A ;\n");
    fprintf(fp,"H410Aux. channel 1 contents     1,NONE ;\n");
    fprintf(fp,"H411Aux. channel 2 contents     1,NONE ;\n");
    fprintf(fp,"H412Aux. channel 3 contents     1,NONE ;\n");
    fprintf(fp,"H413Aux. channel 4 contents     1,NONE ;\n");
    fprintf(fp,"H26 ............................................................................;\n");
    fprintf(fp,"H26                                                                             ;\n");
    fprintf(fp,"H26 Tape Record           Sp      Sp     1st  Last       Rcvr    1st        Rcvr;\n");
    fprintf(fp,"H26 Num     Num           Line  Station  Chan Chan       Line    Rcvr       Indx;\n");
    fprintf(fp,"H26   |       |Rec        |         |     |    |          |         |          |;\n");
    fprintf(fp,"H26   |       |Increment  |         |Sp   |    |Chan      |         |      Last|;\n");
    fprintf(fp,"H26   |       ||          |         |Indx |    |Increment |         |      Rcvr|;\n");
    fprintf(fp,"H26   |       ||Itn       |         ||    |    ||         |         |         ||;\n");
    fprintf(fp,"H26   |       |||         |         ||    |    ||         |         |         ||;\n");
    fprintf(fp,"H26      1         2         3         4         5         6         7         8;\n");
    fprintf(fp,"H26 5678901234567890123456789012345678901234567890123456789012345678901234567890;\n");
    fprintf(fp,"H26   |       |||         |         ||    |    ||         |         |         ||;\n");
    fprintf(fp,"H26   |       |||         |         ||    |    ||         |         |         ||;\n");

    Spstation = 1001.0;
    int i;
    for(i = 0; i < ns; i++){
        fprintf(fp, "X%14d11%10.1f%10.1f%d%5d%5d%d%10.1f%10.1f%10.1f%d\n", 1,1000.1,Spstation,1,1,nr*3,3, 1001.0, 1001.0, 1001.0+(nr-1), 1);
        Spstation += 1.0;
    }
    fclose(fp);



	return 0;
}
