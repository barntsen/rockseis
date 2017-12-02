#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "inversion.h"
#include "inparse.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<InversionElastic3D<float>> inv;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    inv->writeLog("##### Starting new evaluation #####");
    double *x = instance->x;
    double *g = instance->g;

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

    // Compute regularization
    inv->writeLog("Computing regularisation");
    inv->computeRegularisation(x);

    // Combine data and model misfit gradients
    inv->writeLog("Combining gradients");
    inv->combineGradients();

    // Apply mute to gradient
    inv->writeLog("Muting gradients");
    inv->applyMute();

    // Project gradient to B-spline 
    if(inv->getParamtype() == PAR_BSPLINE)
    {
        inv->writeLog("Projecting gradient in B-spline grid");
        task = RUN_BS_PROJ;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
        inv->runBsproj();
    }

    //Read final gradient 
    inv->writeLog("Reading gradient into vector");
    inv->readGrad(g);

    //Read misfit
    inv->writeLog("Reading misfits");
    inv->readMisfit(&instance->f);
    inv->writeLog("##### Evaluation finished #####");

    // Normalize error and gradient
    if(inv->getFnorm() == 0.0){
        inv->setFnorm(instance->f);
    }
    inv->normalize(g, &instance->f, instance->n);

    // Writing progress information to log file
    double xnorm, gnorm, step;
    gnorm = inv->vector_norm(instance->g, 2, instance->n);
    xnorm = inv->vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    char buffer[512],c_iter[32],c_step[32],c_misfit[32],c_gnorm[32],c_mnorm[32];
    snprintf(c_iter,32,"Linesearch\t");
	snprintf(c_step,32,"%15.10e     ",step);
	snprintf(c_misfit,32,"%15.10e     ",instance->f);
	snprintf(c_gnorm,32,"%15.10e     ",gnorm);
	snprintf(c_mnorm,32,"%15.10e     ",xnorm);
    time_t tempo = time(NULL);

	// Creating string for file print
	strcpy(buffer,c_iter);
	strcat(buffer,c_step);
	strcat(buffer,c_misfit);
	strcat(buffer,c_gnorm);
	strcat(buffer,c_mnorm);
	strcat(buffer,ctime(&tempo));
    buffer[strlen(buffer)-1] ='\0';
	inv->writeProgress(buffer);

}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    // Copy new iteration files to results folder
    inv->saveResults(opt->getIter());

    double xnorm, gnorm, step;
    gnorm = inv->vector_norm(instance->g, 2, instance->n);
    xnorm = inv->vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    inv->writeLog("---------- > New iteration found: " + std::to_string(opt->getIter()));

    // Writing progress information to log file
    char buffer[512],c_iter[32],c_step[32],c_misfit[32],c_gnorm[32],c_mnorm[32];
    snprintf(c_iter,32,"Iteration %d\t",opt->getIter());
	snprintf(c_step,32,"%15.10e     ",step);
	snprintf(c_misfit,32,"%15.10e     ",instance->f);
	snprintf(c_gnorm,32,"%15.10e     ",gnorm);
	snprintf(c_mnorm,32,"%15.10e     ",xnorm);
    time_t tempo = time(NULL);

	// Creating string for file print
	strcpy(buffer,c_iter);
	strcat(buffer,c_step);
	strcat(buffer,c_misfit);
	strcat(buffer,c_gnorm);
	strcat(buffer,c_mnorm);
	strcat(buffer,ctime(&tempo));
    buffer[strlen(buffer)-1] ='\0';
	inv->writeProgress(buffer);
}

int main(int argc, char** argv) {

    // Initializing MPI
    MPImodeling mpi(&argc,&argv);
    int task; 


    if(mpi.getNrank() < 2){
        rs_error("This is a parallel program, it must run with at least 2 processors, use mpirun.");
    }

    if(argc < 2){
        if(mpi.getRank() == 0){
            PRINT_DOC(# MPI 3d elastic full-waveform inversion configuration file);
            PRINT_DOC();
            PRINT_DOC(# Modelling parameters);
            PRINT_DOC(freesurface = "true"; # True if free surface should be on);
            PRINT_DOC(order = "8"; # Order of finite difference stencil);
            PRINT_DOC(lpml = "8"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(snapinc = "1"; # Snap interval in multiples of modelling interval);
            PRINT_DOC(apertx = "900"; # Aperture for local model (source is in the middle));
            PRINT_DOC(aperty = "900"; # Aperture for local model (source is in the middle));
            PRINT_DOC(source_type = "3";);
            PRINT_DOC();
            PRINT_DOC(# Checkpointing parameters);
            PRINT_DOC(snapmethod = "1";  # 0- Full checkpointing; 1- Optimal checkpointing);
            PRINT_DOC(nsnaps = "11";);
            PRINT_DOC(incore = "true";);
            PRINT_DOC();
            PRINT_DOC(#Fwi parameters);
            PRINT_DOC(misfit_type = "0";  # 0- Difference; 1- Correlation);
            PRINT_DOC(dataweight = "false";);
            PRINT_DOC(Dataweightfile = "weights.rss";);
            PRINT_DOC(mute = "false";  # Mute gradient and updates);
            PRINT_DOC(Mutefile = "mute.rss"; # File with mute weights);
            PRINT_DOC(max_linesearch = "10"; # maximum number of linesearches);
            PRINT_DOC(max_iterations = "20"; # maximum number of iterations);
            PRINT_DOC(update_vp = "true"; # Update vp);
            PRINT_DOC(update_vs = "true"; # Update vs);
            PRINT_DOC(update_rho = "true"; # Update rho);
            PRINT_DOC(update_source = "false"; # Update source);
            PRINT_DOC();
            PRINT_DOC(# Diagonal scaling parameters);
            PRINT_DOC(kvp = "100.0";);
            PRINT_DOC(kvs = "100.0";);
            PRINT_DOC(krho = "100.0";);
            PRINT_DOC(ksource = "1.0";);
            PRINT_DOC();
            PRINT_DOC(#Parameterisation);
            PRINT_DOC(paramtype = "1";  # 0- grid; 1- B-spline;);
            PRINT_DOC(dtx = "25.0"; # knot sampling in B-spline);
            PRINT_DOC(dty = "25.0"; # knot sampling in B-spline);
            PRINT_DOC(dtz = "25.0"; # knot sampling in B-spline);
            PRINT_DOC();
            PRINT_DOC(#Regularisation);
            PRINT_DOC(vpregalpha = "1.0e-5";);
            PRINT_DOC(vsregalpha = "1.0e-5";);
            PRINT_DOC(rhoregalpha = "1.0e-5";);
            PRINT_DOC();
            PRINT_DOC(# Files);
            PRINT_DOC(Vp = "Vp3d.rss";);
            PRINT_DOC(Vs = "Vs3d.rss";);
            PRINT_DOC(Rho = "Rho3d.rss";);
            PRINT_DOC(Wavelet = "Wav3d.rss";);
            PRINT_DOC(Uxrecordfile = "Vxshot.rss";);
            PRINT_DOC(Uyrecordfile = "Vyshot.rss";);
            PRINT_DOC(Uzrecordfile = "Vzshot.rss";);
            PRINT_DOC(Snapfile = "Local/Snap.rss";);
        }
    }

    // Initialize Inversion class
    inv = std::make_shared<rockseis::InversionElastic3D<float>>(&mpi);
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
    int source_type;
	int misfit_type;
    float apertx;
    float aperty;
    float dtx=-1;
    float dty=-1;
    float dtz=-1;
    float kvp, kvs, krho, ksource;
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
    std::string Dataweightfile;
    std::string Misfitfile;
    std::string Snapfile;
    std::string Uxrecordfile;
    std::string Uxmodelledfile;
    std::string Uxresidualfile;
    std::string Uyrecordfile;
    std::string Uymodelledfile;
    std::string Uyresidualfile;
    std::string Uzrecordfile;
    std::string Uzmodelledfile;
    std::string Uzresidualfile;
    std::string Mutefile;


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
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("kvp", &kvp) == INPARSE_ERR) status = true;
    if(Inpar->getPar("kvs", &kvs) == INPARSE_ERR) status = true;
    if(Inpar->getPar("krho", &krho) == INPARSE_ERR) status = true;
    if(Inpar->getPar("ksource", &ksource) == INPARSE_ERR) status = true;
    if(Inpar->getPar("paramtype", &_paramtype) == INPARSE_ERR) status = true;
    if(Inpar->getPar("source_type", &source_type) == INPARSE_ERR) status = true;
    rockseis::rs_paramtype paramtype = static_cast<rockseis::rs_paramtype>(_paramtype);
    if(paramtype == PAR_BSPLINE){
        if(Inpar->getPar("dtx", &dtx) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dty", &dty) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dtz", &dtz) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Uxrecordfile", &Uxrecordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Uyrecordfile", &Uyrecordfile) == INPARSE_ERR) status = true;
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

    if(Inpar->getPar("dataweight", &dataweight) == INPARSE_ERR) status = true;
    if(dataweight){
        if(Inpar->getPar("Dataweightfile", &Dataweightfile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("mute", &mute) == INPARSE_ERR) status = true;
    if(mute){
        if(Inpar->getPar("Mutefile", &Mutefile) == INPARSE_ERR) status = true;
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
    inv->setUyrecordfile(Uyrecordfile);
    inv->setUzrecordfile(Uzrecordfile);
    inv->setDataweight(dataweight);
    inv->setDataweightfile(Dataweightfile);
    if(mute){
        inv->setMutefile(Mutefile);
    }
    inv->setVpfile(VPLSFILE);
    inv->setVsfile(VSLSFILE);
    inv->setRhofile(RHOLSFILE);
    inv->setWaveletfile(SOURCELSFILE);
    inv->setMisfitfile(MISFITFILE);
    inv->setUxmodelledfile(UXMODFILE);
    inv->setUxresidualfile(UXRESFILE);
    inv->setUymodelledfile(UXMODFILE);
    inv->setUyresidualfile(UXRESFILE);
    inv->setUzmodelledfile(UZMODFILE);
    inv->setUzresidualfile(UZRESFILE);
    inv->setSnapfile(Snapfile);
    inv->setApertx(apertx);
    inv->setAperty(aperty);
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
    inv->setDty(dty);
    inv->setDtz(dtz);

    inv->setVpregalpha(vpregalpha);
    inv->setVsregalpha(vsregalpha);
    inv->setRhoregalpha(rhoregalpha);

    inv->setUpdates(update_vp, update_vs, update_rho, update_source);

    //MASTER
    if(mpi.getRank() == 0){
        // Create a sort class and map over shots
        std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
        Sort->setDatafile(Uxrecordfile);
        Sort->createShotmap(Uxrecordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // L-BFGS configuration
        double *x = nullptr; 
        int N;
        N = inv->setInitial(x, Vpfile, Vsfile, Rhofile, Waveletfile);
        x = (double *) calloc(N, sizeof(double));
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));
        opt->opt_set_initial_guess(x);
        opt->setGtol(0.9);
        opt->setMax_linesearch(max_linesearch);
        opt->setMax_iterations(max_iterations);

        // Create results folder
        inv->createResult();

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        inv->writeLog("Starting optimisation algorithm");
        inv->writeProgress("Maximum number of iterations: " + std::to_string(opt->getMax_iterations()));
        inv->writeProgress("Maximum number of linesearches: " + std::to_string(opt->getMax_linesearch()));
        inv->writeProgress("ITERATION\t   STEP LENGTH            MISFIT               GNORM             MNORM                TIME");
        opt->opt_lbfgs(evaluate, progress);

        // Send message for slaves to quit
        inv->writeLog("Optimisation algorithm finished");
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
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

