#include <valarray>
#include <rsf.hh>
#include "file.h"
#include "data.h"
#include "utils.h"
#include <args.hxx>

#define MAXDIM 8

#define SF_SEGY_FORMAT  24
#define SF_SEGY_NS      20
#define SF_SEGY_DT      16
/*^*/

const int ENDI = 1;
#define IS_BIGENDIAN() ( (*(char*)&ENDI) == 0 )

enum {
    SF_EBCBYTES=3200,	/* Bytes in the card image EBCDIC block */
    SF_BNYBYTES=400,	/* Bytes in the binary coded block	*/
    SF_HDRBYTES=240,	/* Bytes in the tape trace header	*/
    SF_NKEYS=91,	/* Number of mandated header fields	*/
    SF_BHKEYS=27,	/* Number of mandated binary fields	*/
    SF_MAXKEYS=256      /* Maximum number of keys               */
};

typedef unsigned char byte;

static bool little_endian = false;

static byte EBCtoASC[256] = {
    0x00,0x01,0x02,0x03,0xCF,0x09,0xD3,0x7F,
    0xD4,0xD5,0xC3,0x0B,0x0C,0x0D,0x0E,0x0F,
    0x10,0x11,0x12,0x13,0xC7,0xB4,0x08,0xC9,
    0x18,0x19,0xCC,0xCD,0x83,0x1D,0xD2,0x1F,
    0x81,0x82,0x1C,0x84,0x86,0x0A,0x17,0x1B,
    0x89,0x91,0x92,0x95,0xA2,0x05,0x06,0x07,
    0xE0,0xEE,0x16,0xE5,0xD0,0x1E,0xEA,0x04,
    0x8A,0xF6,0xC6,0xC2,0x14,0x15,0xC1,0x1A,
    0x20,0xA6,0xE1,0x80,0xEB,0x90,0x9F,0xE2,
    0xAB,0x8B,0x9B,0x2E,0x3C,0x28,0x2B,0x7C,
    0x26,0xA9,0xAA,0x9C,0xDB,0xA5,0x99,0xE3,
    0xA8,0x9E,0x21,0x24,0x2A,0x29,0x3B,0x5E,
    0x2D,0x2F,0xDF,0xDC,0x9A,0xDD,0xDE,0x98,
    0x9D,0xAC,0xBA,0x2C,0x25,0x5F,0x3E,0x3F,
    0xD7,0x88,0x94,0xB0,0xB1,0xB2,0xFC,0xD6,
    0xFB,0x60,0x3A,0x23,0x40,0x27,0x3D,0x22,
    0xF8,0x61,0x62,0x63,0x64,0x65,0x66,0x67,
    0x68,0x69,0x96,0xA4,0xF3,0xAF,0xAE,0xC5,
    0x8C,0x6A,0x6B,0x6C,0x6D,0x6E,0x6F,0x70,
    0x71,0x72,0x97,0x87,0xCE,0x93,0xF1,0xFE,
    0xC8,0x7E,0x73,0x74,0x75,0x76,0x77,0x78,
    0x79,0x7A,0xEF,0xC0,0xDA,0x5B,0xF2,0xF9,
    0xB5,0xB6,0xFD,0xB7,0xB8,0xB9,0xE6,0xBB,
    0xBC,0xBD,0x8D,0xD9,0xBF,0x5D,0xD8,0xC4,
    0x7B,0x41,0x42,0x43,0x44,0x45,0x46,0x47,
    0x48,0x49,0xCB,0xCA,0xBE,0xE8,0xEC,0xED,
    0x7D,0x4A,0x4B,0x4C,0x4D,0x4E,0x4F,0x50,
    0x51,0x52,0xA1,0xAD,0xF5,0xF4,0xA3,0x8F,
    0x5C,0xE7,0x53,0x54,0x55,0x56,0x57,0x58,
    0x59,0x5A,0xA0,0x85,0x8E,0xE9,0xE4,0xD1,
    0x30,0x31,0x32,0x33,0x34,0x35,0x36,0x37,
    0x38,0x39,0xB3,0xF7,0xF0,0xFA,0xA7,0xFF
};

static byte ASCtoEBC[256] = {
    0x00,0x01,0x02,0x03,0x37,0x2D,0x2E,0x2F,
    0x16,0x05,0x15,0x0B,0x0C,0x0D,0x0E,0x0F,
    0x10,0x11,0x12,0x13,0x3C,0x15,0x32,0x26,
    0x18,0x19,0x3F,0x27,0x1C,0x1D,0x1E,0x1F,
    0x40,0x5A,0x7F,0x7B,0x5B,0x6C,0x50,0x7D,
    0x4D,0x5D,0x5C,0x4E,0x6B,0x60,0x4B,0x61,
    0xF0,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,
    0xF8,0xF9,0x7A,0x5E,0x4C,0x7E,0x6E,0x6F,
    0x7C,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,
    0xC8,0xC9,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,
    0xD7,0xD8,0xD9,0xE2,0xE3,0xE4,0xE5,0xE6,
    0xE7,0xE8,0xE9,0xAD,0xE0,0xBD,0x5F,0x6D,
    0x79,0x81,0x82,0x83,0x84,0x85,0x86,0x87,
    0x88,0x89,0x91,0x92,0x93,0x94,0x95,0x96,
    0x97,0x98,0x99,0xA2,0xA3,0xA4,0xA5,0xA6,
    0xA7,0xA8,0xA9,0xC0,0x4F,0xD0,0xA1,0x07,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,
    0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0x3F,0xFF
};

typedef struct Segy {
    const char *name;
    unsigned int size;
    const char *desc;
} segy;


/* Big-endian to Little-endian conversion and back */
static int convert2(const char* buf);
static int convert4(const char* buf);
static float fconvert4(const char* buf);
static void insert2(int y, char* buf);
static void swapb(byte *x, byte *y);

/* IBM to IEEE float conversion and back */
static float ibm2float (const char* num);
//static void float2ibm (float y, char* num);

static void swapb(byte *x, byte *y) 
    /* swap two bytes */
{
    byte tmp; 

    tmp = *x; 
    *x = *y; 
    *y = tmp;
}

static int convert2(const char* buf)
    /* convert buf to 2-byte int */
{
    union {
        byte b[2];
        short s;
    } x;

    memcpy(x.b,buf,2);

    if (little_endian) swapb(x.b,x.b+1);

    return (int) x.s;
}

static void insert2(int y, char* buf)
    /* convert 2-byte int to buf */
{
    union {
        byte b[2];
        short s;
    } x;

    x.s = (short) y;

    if (little_endian) swapb(x.b,x.b+1);

    memcpy(buf,x.b,2);
}

static int convert4(const char* buf)
    /* convert buf to 4-byte int */
{
    union {
        byte b[4];
        int s;
    } x;

    memcpy(x.b,buf,4);

    if (little_endian) {
        swapb(x.b,x.b+3);
        swapb(x.b+1,x.b+2);
    }

    return x.s;
}

static float fconvert4(const char* buf)
    /* convert buf to 4-byte float */
{
    union {
        byte b[4];
        float s;
    } x;

    memcpy(x.b,buf,4);

    if (little_endian) {
        swapb(x.b,x.b+3);
        swapb(x.b+1,x.b+2);
    }

    return x.s;
}

void ebc2asc (int narr, char* arr)
    /*< Convert char array arrr[narr]: EBC to ASCII >*/
{
    int i;
    unsigned char j;

    for (i=0; i < narr; i++) {
        j = (unsigned char) arr[i];
        arr[i] = (char) EBCtoASC[j];
    }
}

void asc2ebc (int narr, char* arr)
    /*< Convert char array arrr[narr]: ASCII to EBC >*/
{
    int i;
    unsigned char j;

    for (i=0; i < narr; i++) {
        j = (unsigned char) arr[i];
        arr[i] = (char) ASCtoEBC[j];
    }
}

int segyformat (const char* bhead)
    /*< extracts SEGY format from binary header >*/
{
    return convert2(bhead+SF_SEGY_FORMAT);
}

void set_segyformat (char* bhead, int format)
    /*< set SEGY format in binary header >*/
{
    insert2(format,bhead+SF_SEGY_FORMAT);
}

int segyns (const char* bhead)
    /*< extracts ns (number of samples) from binary header >*/
{
    return convert2(bhead+SF_SEGY_NS);
}

void set_segyns(char* bhead, int ns)
    /*< set ns (number of samples) in binary header >*/
{
    insert2(ns,bhead+SF_SEGY_NS);
}

float segydt (const char* bhead)
    /*< extracts dt (sampling) from binary header >*/
{
    return (float) (convert2(bhead+SF_SEGY_DT)/1000000.);
}

void set_segydt(char* bhead, float dt)
    /*< set dt (sampling) in binary header >*/    
{
    float scale;

    /* input in seconds or miliseconds? */
    scale = (dt < 1.0)? 1000000.:1000.;

    insert2((int) (scale*dt),bhead+SF_SEGY_DT);
}

static float ibm2float (const char* num)
    /* floating point conversion from IBM format */
{
    unsigned int x, s, f;
    const unsigned int fMAXIEEE = 0x7F7FFFFF;
    int e;         
    float y;

    x = convert4 (num);

    /* check for special case of zero */
    if ((x & 0x7fffffff) == 0) return 0.0; 

    /* fetch the sign, exponent (removing excess 64), and fraction */   
    s =   x & 0x80000000;                                               
    e = ((x & 0x7f000000) >> 24) - 64;                                   
    f =   x & 0x00ffffff;                                                

    /* convert scale factor from base-16 to base-2 */        
    if (e >= 0) {
        e <<= 2;  
    } else { 
        e = -((-e) << 2); 
    }

    /* convert exponent for 24 bit fraction to 23 bit fraction */           
    e -= 1;                                                               

    /* normalize the fraction */                                            
    if (0 != f) {
        while ((f & 0x00800000) == 0) {         
            f <<= 1;
            e -= 1;
        }       
    }                                                               

    /* drop the '1' preceeding the binary point */                       
    f &= 0x007fffff;                                                         

    /* convert exponent to excess 127 and store the number */
    if ((e += 127) >= 255) {
        s |= fMAXIEEE;
    } else if (e > 0) {
        s |= (e << 23) | f; 	    
    }    

    memcpy (&y,&s,4);
    return y;
}

bool endian (void)
    /*< Set endianness >*/
{
    if(IS_BIGENDIAN()){
        little_endian = false;
    }else{
        little_endian = true;
    }
    return little_endian;
}

static const segy standard_segy_key[] = {
    {"tracl",  4,  "trace sequence number within line 0"},
    {"tracr",  4,  "trace sequence number within reel 4"},
    {"fldr",   4,  "field record number 8"},
    {"tracf",  4,  "trace number within field record 12"},
    {"ep",     4,  "energy source point number 16"},
    {"cdp",    4,  "CDP ensemble number 20"},
    {"cdpt",   4,  "trace number within CDP ensemble 24"},
    {"trid",   2,  "trace identification code:"
        "\n\t1 = seismic data"
            "\n\t2 = dead"
            "\n\t3 = dummy"
            "\n\t4 = time break"
            "\n\t5 = uphole"
            "\n\t6 = sweep"
            "\n\t7 = timing"
            "\n\t8 = water break"
            "\n\t9---, N = optional use (N = 32,767) 28"},
    {"nvs",    2,  "number of vertically summed traces (see vscode in bhed structure) 30"},
    {"nhs",    2,  "number of horizontally summed traces (see vscode in bhed structure) 32"},
    {"duse",   2,  "data use:"
        "\n\t1 = production"
            "\n\t2 = test 34"},
    {"offset", 4,  "distance from source point to receiver group (negative if opposite to direction in which the line was shot) 36"},
    {"gelev",  4,  "receiver group elevation from sea level (above sea level is positive) 40"},
    {"selev",  4,  "source elevation from sea level (above sea level is positive) 44"},
    {"sdepth", 4,  "source depth (positive) 48"},
    {"gdel",   4,  "datum elevation at receiver group 52"},
    {"sdel",   4,  "datum elevation at source 56"},
    {"swdep",  4,  "water depth at source 60"},
    {"gwdep",  4,  "water depth at receiver group 64"},
    {"scalel", 2,  "scale factor for previous 7 entries with value plus or minus 10 to the power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) 68"},
    {"scalco", 2,  "scale factor for next 4 entries with value plus or minus 10 to the power 0, 1, 2, 3, or 4 (if positive, multiply, if negative divide) 70"},
    {"sx",     4,  "X source coordinate 72"},
    {"sy",     4,  "Y source coordinate 76"},
    {"gx",     4,  "X group coordinate 80"},
    {"gy",     4,  "Y source coordinate 84"},
    {"counit", 2,  "coordinate units code:"
        "\n\tfor previoius four entries"
            "\n\t1 = length (meters or feet)"
            "\n\t2 = seconds of arc (in this case, the"
            "\n\tX values are unsigned intitude and the Y values"
            "\n\tare latitude, a positive value designates"
            "\n\tthe number of seconds east of Greenwich or north of the equator 88"},
    {"wevel",   2,  "weathering velocity 90"},
    {"swevel",  2,  "subweathering velocity 92"},
    {"sut",     2,  "uphole time at source 94"},
    {"gut",     2,  "uphole time at receiver group 96"},
    {"sstat",   2,  "source static correction 98"},
    {"gstat",   2,  "group static correction 100"},
    {"tstat",   2,  "total static applied 102"},
    {"laga",    2,  "lag time A, time in ms between end of 240-byte trace identification header and time "
        "break, positive if time break occurs after end of header, time break is defined as "
            "the initiation pulse which maybe recorded on an auxiliary trace or as otherwise "
            "specified by the recording system 104"},
    {"lagb",    2,  "lag time B, time in ms between the time break and the initiation time of the energy source, "
        "may be positive or negative 106"},
    {"delrt",   2,  "delay recording time, time in ms between initiation time of energy source and time "
        "when recording of data samples begins (for deep water work if recording does not start at zero time) 108"},
    {"muts",    2,  "mute time--start 110"},
    {"mute",    2,  "mute time--end 112"},
    {"ns",      2,  "number of samples in this trace 114"},
    {"dt",      2,  "sample interval, in micro-seconds 116"},
    {"gain",    2,  "gain type of field instruments code:"
        "\n\t1 = fixed"
            "\n\t2 = binary"
            "\n\t3 = floating point"
            "\n\t4 ---- N = optional use 118"},
    {"igc",    2,   "instrument gain constant 120"},
    {"igi",    2,   "instrument early or initial gain 122"},
    {"corr",   2,   "correlated:"
        "\n\t1 = no"
            "\n\t2 = yes 124"},    
    {"sfs",    2,   "sweep frequency at start 126"},
    {"sfe",    2,   "sweep frequency at end 128"},
    {"slen",   2,   "sweep length in ms 130"},
    {"styp",   2,   "sweep type code:"
        "\n\t1 = linear"
            "\n\t2 = cos-squared"
            "\n\t3 = other 132"},   
    {"stas",   2,   "sweep trace length at start in ms 134"},
    {"stae",   2,   "sweep trace length at end in ms 136"},
    {"tatyp",  2,   "taper type: 1=linear, 2=cos^2, 3=other 138"},
    {"afilf",  2,   "alias filter frequency if used 140"},
    {"afils",  2,   "alias filter slope 142"},
    {"nofilf", 2,   "notch filter frequency if used 144"},
    {"nofils", 2,   "notch filter slope 146"},
    {"lcf",    2,   "low cut frequency if used 148"},
    {"hcf",    2,   "high cut frequncy if used 150"},
    {"lcs",    2,   "low cut slope 152"},
    {"hcs",    2,   "high cut slope 154"},
    {"year",   2,   "year data recorded 156"},
    {"day",    2,   "day of year 158"},
    {"hour",   2,   "hour of day (24 hour clock) 160"},
    {"minute", 2,   "minute of hour 162"},
    {"sec",    2,   "second of minute 164"},
    {"timbas", 2,   "time basis code:"
        "\n\t1 = local"
            "\n\t2 = GMT"
            "\n\t3 = other 166"},   
    {"trwf",   2,   "trace weighting factor, defined as 1/2^N volts for the least sigificant bit 168"},
    {"grnors", 2,   "geophone group number of roll switch position one 170"},
    {"grnofr", 2,   "geophone group number of trace one within original field record 172"},
    {"grnlof", 2,   "geophone group number of last trace within original field record 174"},
    {"gaps",   2,   "gap size (total number of groups dropped) 176"},
    {"otrav",  2,   "overtravel taper code:"
        "\n\t1 = down (or behind)"
            "\n\t2 = up (or ahead) 71/178"},
    {"cdpx",   4,   "X coordinate of CDP 180"},
    {"cdpy",   4,   "Y coordinate of CDP 184"},
    {"iline",  4,   "in-line number 188"},
    {"xline",  4,   "cross-line number 192"},
    {"shnum",  4,   "shotpoint number 196"},
    {"shsca",  2,   "shotpoint scalar 200"},
    {"tval",   2,   "trace value meas. 202"},
    {"tconst4",4,   "transduction const 204"},
    {"tconst2",2,   "transduction const 208"},
    {"tunits", 2,   "transduction units 210"},
    {"device", 2,   "device identifier 212"},
    {"tscalar",2,   "time scalar 214"},
    {"stype",  2,   "source type 216"},
    {"sendir", 4,   "source energy dir. 218"},
    {"unknown",2,   "unknown 222"},
    {"smeas4", 4,   "source measurement 224"},
    {"smeas2", 2,   "source measurement 228"},
    {"smeasu", 2,   "source measurement unit 230"},
    {"unass1", 4,    "unassigned 232"},
    {"unass2", 4,    "unassigned 236"},
};

int segykey (const char* key, char *trace) 
    /*< Extract a SEGY key value >*/
{
    int i;

    int count = 0;
    bool found = false;
    for (i=0; i < SF_NKEYS; i++) {
        if (0==strcmp(key,standard_segy_key[i].name)){
            found = true;
            break;
        }else{
            count += standard_segy_key[i].size;
        }
    }
    if(!found){
        std::cerr << "no such key " << key << std::endl;
        exit(1);
    }
    switch (standard_segy_key[i].size ){
        case 2:
            return convert2(&trace[0] + count);
            break;
        case 4:
            return convert4(&trace[0] + count);
            break;
    }
    return 0;
}

void segy2trace(const char* buf, float* trace, int ns, int format)
/*< Extract a floating-point trace[nt] from buffer buf.
---
format: 1: IBM, 2: int4, 3: int2, 5: IEEE
>*/
{
    int i, nb;

    nb = (3==format)? 2:4;

    for (i=0; i < ns; i++, buf += nb) {
        switch (format) {
            case 1: trace[i] = ibm2float (buf);       break; /* IBM float */
            case 2: trace[i] = (float) convert4(buf); break; /* int4 */
            case 3: trace[i] = (float) convert2(buf); break; /* int2 */
            case 5: trace[i] = fconvert4(buf);        break; /* IEEE float */
            default: 
                    std::cerr <<"Unknown format " << format; 
                    exit(1);
                    break;
        }
    }
}

int main(int argc, char* argv[])
{
    /* Variables */
    char ahead[SF_EBCBYTES], bhead[SF_BNYBYTES];
    int format, ns; 
    size_t ntr;
    char *trace;
    float dt=0.0, t0;
    FILE *file;
    off_t pos, start, nsegy=0;
    bool status = false;
    std::string Datatype;

    args::ArgumentParser parser("Program to convert RSF file to Data RSS file.", "");
    parser.LongPrefix("");
    parser.LongSeparator("=");
    args::HelpFlag help(parser, "help", "Display this help menu", {"h", "help"});
    args::ValueFlag<std::string> partype(parser, "DATA2D or DATA3D", "Dimension", {"type"});
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
    if (partype){
         Datatype = args::get(partype);
    }else{
        std::cerr << parser;
        rockseis::rs_error("Missing type.");
    }

    int type = -1; 
    if(!Datatype.compare("DATA2D"))
    {
        type = 2;
    }
    if(!Datatype.compare("DATA3D"))
    {
        type = 3;
    }
    if(type == -1)
    {
        type = 3;
    }

    rockseis::rs_datatype dtype = static_cast<rockseis::rs_datatype>(type);

    //Using stdin as input 
    file = stdin;

    //Determine endian
    endian();

    // Getting filesize
    ntr=0; // For now
    if (0==ntr) {
        fseeko(file,0,SEEK_END);
        pos = ftello(file); /* pos is the filesize in bytes */
        fseeko(file,0,SEEK_SET);
    }else{
        pos=0;
    }

    if (SF_EBCBYTES != fread(ahead, 1, SF_EBCBYTES, file)){ 
        std::cerr <<  "Error reading ascii/ebcdic header." << std::endl;
        exit(1);
    }

    if (SF_BNYBYTES != fread(bhead, 1, SF_BNYBYTES, file)){
        std::cerr <<  "Error reading binary header." << std::endl;
        exit(1);
    }

    //Getting parameters from binary header
    format = segyformat (bhead);

    ns = segyns (bhead);
    if (0>=ns){
        std::cerr << "Number of samples is not set in binary header" << std::endl;
        exit(1);
    }
    dt = segydt (bhead);

    nsegy = SF_HDRBYTES + ((3 == format)? ns*2: ns*4);    

    if (0==ntr) ntr = (pos - SF_EBCBYTES - SF_BNYBYTES)/nsegy;

    start = ftello(file);

    /* read first trace header */
    trace = sf_charalloc (SF_HDRBYTES);
    if (SF_HDRBYTES != fread(trace, 1, SF_HDRBYTES, file)){
        std::cerr<< "Error reading first trace header" << std::endl;
        exit(1);
    }
    fseeko(file,start,SEEK_SET); /* rewind */
    int scalel, scalco;
    int sx, sy, sz, gx, gy, gz;
    scalel = segykey("scalel", trace);
    scalco = segykey("scalco", trace);
    float fscalco, fscalel;
    bool sign; 
    sign = (scalco < 0);
    if (sign) {
        fscalco = -1./scalco;
    }else{
        fscalco = (float) scalco;
    }

    sign = (scalel < 0);
    if (sign) {
        fscalel = -1./scalel;
    }else{
        fscalel = (float) scalel;
    }

    t0 = segykey("delrt", trace)/1000.;
    free(trace);

    std::shared_ptr<rockseis::Data2D<float>> Bdata2D = NULL;
    std::shared_ptr<rockseis::Data3D<float>> Bdata3D = NULL;

    rockseis::Point2D<float> *scoords2D = NULL;
    rockseis::Point2D<float> *gcoords2D = NULL;
    rockseis::Point3D<float> *scoords3D = NULL;
    rockseis::Point3D<float> *gcoords3D = NULL;
    float *ftrace;

    switch(dtype){
        case rockseis::DATA2D:
            Bdata2D = std::make_shared<rockseis::Data2D<float>>(1, ns, dt, t0);
            scoords2D = (Bdata2D->getGeom())->getScoords();
            gcoords2D = (Bdata2D->getGeom())->getGcoords();
            //Using stdout as output 
            Bdata2D->setFile("stdout");
            // Open data for output
            status = Bdata2D->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

            nsegy = SF_HDRBYTES + ns*4;
            trace = sf_charalloc (nsegy);
            ftrace = Bdata2D->getData();
            for (size_t i=0; i < ntr; i++)
            {
                if(nsegy != fread(trace, 1, nsegy, file)){
                    std::cerr << "Error reading data in segy file.\n";
                    exit(1);
                }
                // Get coordinates from header
                sx = segykey("sx", trace);
                sz = segykey("selev", trace);
                gx = segykey("gx", trace);
                gz = segykey("gelev", trace);
                scoords2D[0].x = fscalco * sx;
                scoords2D[0].y = fscalel * sz;
                gcoords2D[0].x = fscalco * gx;
                gcoords2D[0].y = fscalel * gz;
                segy2trace(trace + SF_HDRBYTES, ftrace, ns, format);
                Bdata2D->writeTraces();
            }	
            Bdata2D->close();

            break;
        case rockseis::DATA3D:
            Bdata3D = std::make_shared<rockseis::Data3D<float>>(1, ns, dt, t0);
            scoords3D = (Bdata3D->getGeom())->getScoords();
            gcoords3D = (Bdata3D->getGeom())->getGcoords();
            //Using stdout as output 
            Bdata3D->setFile("stdout");
            // Open data for output
            status = Bdata3D->open("o");
            if(status == FILE_ERR) rockseis::rs_error("Error opening file for writting");

            nsegy = SF_HDRBYTES + ns*4;
            trace = sf_charalloc (nsegy);
            ftrace = Bdata3D->getData();
            for (size_t i=0; i < ntr; i++)
            {
                if(nsegy != fread(trace, 1, nsegy, file)){
                    std::cerr << "Error reading data in segy file.\n";
                    exit(1);
                }
                // Get coordinates from header
                sx = segykey("sx", trace);
                sy = segykey("sy", trace);
                sz = segykey("selev", trace);
                gx = segykey("gx", trace);
                gy = segykey("gy", trace);
                gz = segykey("gelev", trace);
                scoords3D[0].x = fscalco * sx;
                scoords3D[0].y = fscalco * sy;
                scoords3D[0].z = fscalel * sz;
                gcoords3D[0].x = fscalco * gx;
                gcoords3D[0].y = fscalco * gy;
                gcoords3D[0].z = fscalel * gz;
                segy2trace(trace + SF_HDRBYTES, ftrace, ns, format);
                Bdata3D->writeTraces();
            }	
            Bdata3D->close();
            break;
        default:
            rockseis::rs_error("Invalid data type.");
            break;
    }

    free (trace);
    exit (0);
}

