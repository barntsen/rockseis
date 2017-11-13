#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "inversion.h"
#include "inparse.h"

// TODO Add regularization
//
using namespace rockseis;

/* Global variables */
std::shared_ptr<InversionAcoustic2D<float>> inv;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    fprintf(stderr, "Starting new evaluation\n");
    double *x = instance->x;
    double *g = instance->g;

    // Save linesearch model
    inv->saveLinesearch(x);

    // Start new gradient evaluation
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv->runGrad();

    // Apply mute to gradient
    inv->applyMute();

    if(inv->getParamtype() == PAR_BSPLINE)
    {
        task = RUN_BS_PROJ;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
        inv->runBsproj();
    }

    //Read gradient
    inv->readGrad(g);

    //Read misfit
    inv->readMisfit(&instance->f);
}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{

    double xnorm, gnorm, step;
    gnorm = opt->opt_vector_norm(instance->g, 2, instance->n);
    xnorm = opt->opt_vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    fprintf(stderr, "Iteration %d:\n", opt->getIter());
    fprintf(stderr, "  fx = %f \n", instance->f);
    fprintf(stderr, "  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
}

int main(int argc, char** argv) {

    // Initializing MPI
    MPImodeling mpi(&argc,&argv);
    int task; 

    // Initialize Inversion class
    inv = std::make_shared<rockseis::InversionAcoustic2D<float>>(&mpi);
	/* General input parameters */
    bool status;
	int lpml;
	bool fs;
    bool incore = false;
    bool dataweight;
    bool mute;
	int order;
	int snapinc;
	int nsnaps = 0;
	int _snapmethod;
    int _paramtype;
	int misfit_type;
    float apertx;
    float dtx=-1;
    float dtz=-1;
    float kvp, krho, ksource;
    std::string Waveletfile;
    std::string Vpfile;
    std::string Rhofile;
    std::string Vpgradfile;
    std::string Rhogradfile;
    std::string Wavgradfile;
    std::string Dataweightfile;
    std::string Misfitfile;
    std::string Psnapfile;
    std::string Precordfile;
    std::string Pmodelledfile;
    std::string Presidualfile;
    std::string Mutefile;


    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("lpml", &lpml) == INPARSE_ERR) status = true;
    if(Inpar->getPar("order", &order) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapinc", &snapinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("freesurface", &fs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rho", &Rhofile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavelet", &Waveletfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("kvp", &kvp) == INPARSE_ERR) status = true;
    if(Inpar->getPar("krho", &krho) == INPARSE_ERR) status = true;
    if(Inpar->getPar("ksource", &ksource) == INPARSE_ERR) status = true;
    if(Inpar->getPar("paramtype", &_paramtype) == INPARSE_ERR) status = true;
    rockseis::rs_paramtype paramtype = static_cast<rockseis::rs_paramtype>(_paramtype);
    if(paramtype == PAR_BSPLINE){
        if(Inpar->getPar("dtx", &dtx) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dtz", &dtz) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
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

    if(Inpar->getPar("dataweight", &dataweight) == INPARSE_ERR) status = true;
    if(dataweight){
        if(Inpar->getPar("Dataweightfile", &Dataweightfile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("mute", &mute) == INPARSE_ERR) status = true;
    if(mute){
        if(Inpar->getPar("Mutefile", &Mutefile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    inv->setOrder(order);
    inv->setLpml(lpml);
    inv->setFs(fs);
    inv->setSnapinc(snapinc);

    inv->setPrecordfile(Precordfile);
    inv->setDataweight(dataweight);
    inv->setDataweightfile(Dataweightfile);
    if(mute){
        inv->setMutefile(Mutefile);
    }

    inv->setVpfile(VPLSFILE);
    inv->setRhofile(RHOLSFILE);
    inv->setWaveletfile(SOURCELSFILE);
    inv->setMisfitfile(MISFITFILE);
    inv->setPmodelledfile(PMODFILE);
    inv->setPresidualfile(PRESFILE);
    inv->setPsnapfile(PSNAPFILE);
    inv->setApertx(apertx);
    inv->setSnapmethod(snapmethod);
    inv->setNsnaps(nsnaps);
    inv->setIncore(incore);
    inv->setMisfit_type(fwimisfit);

    inv->setVpgradfile(VPGRADFILE);
    inv->setRhogradfile(RHOGRADFILE);
    inv->setWavgradfile(SOURCEGRADFILE);
    inv->setKvp(kvp);
    inv->setKrho(krho);
    inv->setKsource(ksource);
    inv->setParamtype(paramtype);
    inv->setDtx(dtx);
    inv->setDtz(dtz);

    //MASTER
    if(mpi.getRank() == 0){
        // Create a sort class and map over shots
        std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
        Sort->setDatafile(Precordfile);
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // L-BFGS configuration
        /* Initialize the parameters for the L-BFGS optimization. */
        double *x = nullptr; 
        int N;
        N = inv->setInitial(x, Vpfile, Rhofile, Waveletfile);
        x = (double *) calloc(N, sizeof(double));
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));
        opt->opt_set_initial_guess(x);
        
        opt->setGtol(0.9);
        opt->setMax_linesearch(5);

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        //
        opt->opt_lbfgs(evaluate, progress);

        // Send message for slaves to quit
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
        fprintf(stderr, "%s", buffer);

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
                    break;
                case RUN_BS_PROJ:
                    inv->runBsproj();
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

