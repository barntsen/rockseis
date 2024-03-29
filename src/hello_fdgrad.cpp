#include <iostream>
#include <mpi.h>
#include "inversion.h"
#include "inparse.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<InversionElastic2D<float>> inv;

/* Global functions */
void evaluate(double *f, double *x)
{
    inv->writeLog("##### Starting new evaluation #####");

    inv->writeLog("Saving linesearch models");
    // Save linesearch model
    inv->saveLinesearch(x);
    inv->writeLog("Linesearch models saved");

    // Start new gradient evaluation
    inv->writeLog("Starting gradient computation");
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv->runGrad();
    inv->writeLog("Gradient computation finished");

        //Read misfit
    inv->writeLog("Reading misfits");
    inv->readMisfit(f);
    inv->writeLog("##### Evaluation finished #####");

    // Writing progress information to log file
    char buffer[512],c_iter[32],c_misfit[32];
    snprintf(c_iter,32,"Linesearch\t");
	snprintf(c_misfit,32,"%15.10e     ",*f);
    time_t tempo = time(NULL);

	// Creating string for file print
	strcpy(buffer,c_iter);
	strcat(buffer,c_misfit);
	strcat(buffer,ctime(&tempo));
    buffer[strlen(buffer)-1] ='\0';
	inv->writeProgress(buffer);
}

int main(int argc, char** argv) {

    // Initializing MPI
    MPImodeling mpi(&argc,&argv);
    int task; 

    // Initialize Inversion class
    inv = std::make_shared<rockseis::InversionElastic2D<float>>(&mpi);
	/* General input parameters */
    bool status;
	int lpml;
	bool fs;
    bool incore = false;
    bool dataweightx;
    bool dataweightz;
    bool mute;
	int order;
	int snapinc;
	int nsnaps = 0;
	int _snapmethod;
    int source_type;
	int misfit_type;
    float apertx;
    float dtx=-1;
    float dtz=-1;
    float kvp = 1.0, kvs = 1.0, krho=1.0, ksource=1.0;
    float vpregalpha, vsregalpha, rhoregalpha;
    int max_linesearch, max_iterations;
    bool update_vp, update_vs, update_rho, update_source;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Vsfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Vsgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightxfile;
    std::string Dataweightzfile;
    std::string Misfitfile;
    std::string Snapfile;
    std::string Uxrecordfile;
    std::string Uxmodelledfile;
    std::string Uxresidualfile;
    std::string Uzrecordfile;
    std::string Uzmodelledfile;
    std::string Uzresidualfile;
    std::string Modmutefile;


    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file ", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("order", &order) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vs", &Vsfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("source_type", &source_type) == INPARSE_ERR) status = true;
    rockseis::rs_paramtype paramtype = PAR_GRID;
    if(Inpar->getPar("Uxrecordfile", &Uxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uzrecordfile", &Uzrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Snapfile", &Snapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapmethod", &_snapmethod) == INPARSE_ERR) status = true;
    rockseis::rs_snapmethod snapmethod = static_cast<rockseis::rs_snapmethod>(_snapmethod);
    switch(snapmethod){
        case rockseis::FULL:
            break;
        case rockseis::OPTIMAL:
            if(Inpar->getPar("nsnaps", &nsnaps) == INPARSE_ERR) status = true;
            if(Inpar->getPar("incore", &incore) == INPARSE_ERR) status = true;
            break;
        default:
            rockseis::rs_error("Invalid option of snapshot saving (snapmethod)."); 
    }
    if(Inpar->getPar("misfit_type", &misfit_type) == INPARSE_ERR) status = true;
    rockseis::rs_fwimisfit fwimisfit = static_cast<rockseis::rs_fwimisfit>(misfit_type);

    if(Inpar->getPar("dataweightx", &dataweightx) == INPARSE_ERR) status = true;
    if(dataweightx){
        if(Inpar->getPar("Dataweightxfile", &Dataweightxfile) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("dataweightz", &dataweightz) == INPARSE_ERR) status = true;
    if(dataweightz){
        if(Inpar->getPar("Dataweightzfile", &Dataweightzfile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("mute", &mute) == INPARSE_ERR) status = true;
    if(mute){
        if(Inpar->getPar("Modmutefile", &Modmutefile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("vpregalpha", &vpregalpha) == INPARSE_ERR) status = true;
    if(Inpar->getPar("vsregalpha", &vsregalpha) == INPARSE_ERR) status = true;
    if(Inpar->getPar("rhoregalpha", &rhoregalpha) == INPARSE_ERR) status = true;
    if(Inpar->getPar("max_linesearch", &max_linesearch) == INPARSE_ERR) status = true;
    if(Inpar->getPar("max_iterations", &max_iterations) == INPARSE_ERR) status = true;

    if(Inpar->getPar("update_vp", &update_vp) == INPARSE_ERR) status = true;
    if(Inpar->getPar("update_vs", &update_vs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("update_rho", &update_rho) == INPARSE_ERR) status = true;
    if(Inpar->getPar("update_source", &update_source) == INPARSE_ERR) status = true;

    // Set scaling according to updates
    if(!update_vp) kvp = 0.0;
    if(!update_vs) kvs = 0.0;
    if(!update_rho) krho = 0.0;
    if(!update_source) ksource = 0.0;

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    inv->setOrder(order);
    inv->setLpml(lpml);
    inv->setFs(fs);
    inv->setSnapinc(snapinc);

    inv->setUxrecordfile(Uxrecordfile);
    inv->setUzrecordfile(Uzrecordfile);
    inv->setDataweightx(dataweightx);
    inv->setDataweightz(dataweightz);
    inv->setDataweightxfile(Dataweightxfile);
    inv->setDataweightzfile(Dataweightzfile);
    if(mute){
        inv->setModmutefile(Modmutefile);
    }
    inv->setVpfile(VPLSFILE);
    inv->setVsfile(VSLSFILE);
    inv->setRhofile(RHOLSFILE);
    inv->setWaveletfile(SOURCELSFILE);
    inv->setMisfitfile(MISFITFILE);
    inv->setUxmodelledfile(UXMODFILE);
    inv->setUxresidualfile(UXRESFILE);
    inv->setUzmodelledfile(UZMODFILE);
    inv->setUzresidualfile(UZRESFILE);
    inv->setSnapfile(Snapfile);
    inv->setApertx(apertx);
    inv->setSnapmethod(snapmethod);
    inv->setNsnaps(nsnaps);
    inv->setIncore(incore);
    inv->setMisfit_type(fwimisfit);

    inv->setVpgradfile(VPGRADFILE);
    inv->setVsgradfile(VSGRADFILE);
    inv->setRhogradfile(RHOGRADFILE);
    inv->setWavgradfile(SOURCEGRADFILE);
    inv->setKvp(kvp);
    inv->setKvs(kvs);
    inv->setKrho(krho);
    inv->setKsource(ksource);
    inv->setParamtype(paramtype);
    inv->setSourcetype(source_type);
    inv->setDtx(dtx);
    inv->setDtz(dtz);

    inv->setVpregalpha(0.0);
    inv->setVsregalpha(0.0);
    inv->setRhoregalpha(0.0);

    inv->setUpdates(update_vp, update_vs, update_rho, update_source);

    //MASTER
    if(mpi.getRank() == 0){
        // Create a sort class and map over shots
        std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
        Sort->setDatafile(Uxrecordfile);
        Sort->createShotmap(Uxrecordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        inv->writeLog("Starting finite difference gradient computation");
        inv->writeProgress("ITERATION\t   MISFIT               TIME");

        double f0,f,df;
        df= 1.0e10;
        int Nx=45; 
        int Nz=45;
        double *x = nullptr; 
        int N;
        int Ng=51;
        N = inv->setInitial(x, Vpfile, Vsfile, Rhofile, Waveletfile);
        x = (double *) calloc(N, sizeof(double));
        // Compute initial error
        evaluate(&f0, x);
        inv->setNoreverse(true);
        float *grad = (float *) calloc(Ng, sizeof(float));
        rockseis::Index I(Nx,Nz);
        for(int i=0; i < Ng; i++)
        {
            //x[I(22,i)] = df;
            x[i*20] = df;
            evaluate(&f, x);
            //x[I(22,i)]=0.0;
            x[i*20] = 0.0;
            grad[i]= (f-f0)/df;
        }

        // Misfit file creation
        std::shared_ptr<rockseis::File> Ffdgrad (new rockseis::File());
        Ffdgrad->output("Fdgrad.rss");
        Ffdgrad->setN(1,Ng);
        Ffdgrad->setD(1,1.0);
        Ffdgrad->setData_format(sizeof(float));
        Ffdgrad->write(grad, Ng, 0);
        Ffdgrad->close();

        // Send message for slaves to quit
        inv->writeLog("Optimisation algorithm finished");
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"Finite difference gradient finished.\n");
        inv->writeProgress(buffer);

        // Free initial model
        free(x);

    //SLAVE
    }else{
        bool stop = false;
        while(1)
        {
            MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
            switch(task)
            {
                case RUN_F_GRAD:
                    inv->runGrad();
                    inv->setNoreverse(true);
                    break;
                case BREAK_LOOP:
                    stop = true;
                    break;
            }
            if(stop){
                break;
            }
        }

    }

    return 0;
}

