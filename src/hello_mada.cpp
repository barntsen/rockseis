#include <valarray>

#include <rsf.hh>

int main(int argc, char* argv[])
{
    sf_init(argc,argv);

    madagascar::iRSF par(0);
    madagascar::oRSF out;
    int n1 = 10;
    int n2 = 5;
    
    std::valarray<int> trace(n1);
    for (int i=0; i < n1; i++){
        trace[i] = 0;
    }
    
    trace[5]=1; 

    float d1=1.0;
    float o1=0.0;
    float d2=1.0;
    float o2=0.0;
    out.type(SF_INT);
    out.put("n1",n1);
    out.put("d1",d1);
    out.put("o1",o1);
    out.put("n2",n2);
    out.put("d2",d2);
    out.put("o2",o2);

    for (int i=0; i < 5; i++) {
	out << trace;
    }

    exit (0);
}

