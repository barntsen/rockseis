#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "lsmiginv.h"
#include "inparse.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<LsmiginvAcoustic2D<float>> lsmig;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    lsmig->writeLog("##### Starting new evaluation #####");
    double *x = instance->x;
    double *g = instance->g;

    lsmig->writeLog("Saving linesearch models");
    // Save linesearch model
    lsmig->saveLinesearch(x);
    lsmig->writeLog("Linesearch models saved");

    // Start new gradient evaluation
    lsmig->writeLog("Starting gradient computation");
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    lsmig->runGrad();
    lsmig->writeLog("Gradient computation finished");

    // Apply source ilumination
    lsmig->writeLog("Applying source ilumination correction");
    lsmig->applySrcilum();

    // Compute regularization
    lsmig->writeLog("Computing regularisation");
    lsmig->computeRegularisation(x);

    // Combine data and model misfit gradients
    lsmig->writeLog("Combining gradients");
    lsmig->combineGradients();

    // Apply mute to gradient
    lsmig->writeLog("Muting gradients");
    lsmig->applyMute();

    //Read final gradient 
    lsmig->writeLog("Reading gradient into vector");
    lsmig->readGrad(g);

    //Read misfit
    lsmig->writeLog("Reading misfits");
    lsmig->readMisfit(&instance->f);
    lsmig->writeLog("##### Evaluation finished #####");

    double xnorm, gnorm, step;

    // Set kvp and normalize error and gradient 
    if(lsmig->getFnorm() == 0.0){

        gnorm = lsmig->vector_norm(instance->g, 2, instance->n);
        // Reset kvp 
        lsmig->setKvp(lsmig->getKvp()*lsmig->getKvp()*instance->f/gnorm);
        // Re-read gradient with new kvp set
        lsmig->readGrad(g);
        lsmig->setFnorm(instance->f);
    }
    lsmig->normalize(g, &instance->f, instance->n);

    gnorm = lsmig->vector_norm(instance->g, 2, instance->n);
    xnorm = lsmig->vector_norm(instance->x, 2, instance->n);
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
	lsmig->writeProgress(buffer);
}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    // Copy new iteration files to results folder
    lsmig->saveResults(opt->getIter());

    double xnorm, gnorm, step;
    gnorm = lsmig->vector_norm(instance->g, 2, instance->n);
    xnorm = lsmig->vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    lsmig->writeLog("---------- > New iteration found: " + std::to_string(opt->getIter()));

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
	lsmig->writeProgress(buffer);
}

void finalize(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    // Do nothing
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
            PRINT_DOC(# MPI 2d acoustic full-waveform inversion configuration file);
            PRINT_DOC();
            PRINT_DOC(# Modelling parameters);
            PRINT_DOC(freesurface = "false"; # True if free surface should be on);
            PRINT_DOC(order = "8"; # Order of finite difference stencil);
            PRINT_DOC(lpml = "8"; # Size of pml absorbing boundary (should be larger than order + 5 ));
            PRINT_DOC(snapinc = "1"; # Snap interval in multiples of modelling interval);
            PRINT_DOC(apertx = "900"; # Aperture for local model (source is in the middle));
            PRINT_DOC();
            PRINT_DOC(# Checkpointing parameters);
            PRINT_DOC(snapmethod = "1";  # 0- Full checkpointing; 1- Optimal checkpointing);
            PRINT_DOC(nsnaps = "11";);
            PRINT_DOC(incore = "true";);
            PRINT_DOC();
            PRINT_DOC(#Fwi parameters);
            PRINT_DOC(misfit_type = "0";  # 0- Difference; 1- Correlation; 2- Adaptive with Gaussian; 3-Adaptive with linear);
            PRINT_DOC(dataweight = "false";);
            PRINT_DOC(Dataweightfile = "weights.rss";);
            PRINT_DOC(mute = "false";  # Mute gradient and updates);
            PRINT_DOC(Mutefile = "mute.rss"; # File with mute weights);
            PRINT_DOC(srcilum = "false"; # Correct gradient for source ilumination);
            PRINT_DOC(max_linesearch = "5"; # maximum number of linesearches);
            PRINT_DOC(max_iterations = "20"; # maximum number of iterations);

            PRINT_DOC(optmethod = "1"; # 1-L-BFGS; 2-CG_FR; 3-STEEPEST DESCENT; 4-CG_PR);
            PRINT_DOC(linesearch = "3"; # 1-Decrease; 2-Armijo; 3-Wolfe; 4-Strong Wolfe);
            PRINT_DOC();
            PRINT_DOC(#Filter parameters);
            PRINT_DOC(filter = "false"; # Apply 4 point filter to wavelet gradient and residuals);
            PRINT_DOC(f0 = "0"; # First point);
            PRINT_DOC(f1 = "2"; # Second point);
            PRINT_DOC(f2 = "10"; # Third point);
            PRINT_DOC(f3 = "12.5"; # Fourth point);
            PRINT_DOC();
            PRINT_DOC(#Preconditioning);
            PRINT_DOC(zder = "true"; # Z-derivative preconditioning);
            PRINT_DOC();
            PRINT_DOC(#Gradient scaling parameter);
            PRINT_DOC(kvp = "1.0";);
            PRINT_DOC();
            PRINT_DOC(#Regularisation);
            PRINT_DOC(vpregalpha = "0.0";);
            PRINT_DOC();
            PRINT_DOC(# Files);
            PRINT_DOC(Vp = "Vp2d.rss";);
            PRINT_DOC(Rho = "Rho2d.rss";);
            PRINT_DOC(Wavelet = "Wav2d.rss";);
            PRINT_DOC(Precordfile = "Pshot.rss";);
            PRINT_DOC(Psnapfile = "Local/Psnap.rss";);
        }
        exit(1);
    }

    // Initialize Lsmiginv class
    lsmig = std::make_shared<rockseis::LsmiginvAcoustic2D<float>>(&mpi);
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
	int misfit_type;
    float apertx;
    float kvp;
    float vpregalpha;
    float freqs[4];
    bool filter;
    int max_linesearch, max_iterations;
    int linesearch;
    int optmethod; 
    bool srcilum;
    bool zder;
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
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
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
    if(Inpar->getPar("max_linesearch", &max_linesearch) == INPARSE_ERR) status = true;
    if(Inpar->getPar("max_iterations", &max_iterations) == INPARSE_ERR) status = true;

    if(Inpar->getPar("srcilum", &srcilum) == INPARSE_ERR) status = true;
    if(Inpar->getPar("kvp", &kvp) == INPARSE_ERR) status = true;

    if(Inpar->getPar("linesearch", &linesearch) == INPARSE_ERR) status = true;
    if(Inpar->getPar("optmethod", &optmethod) == INPARSE_ERR) status = true;

    if(Inpar->getPar("filter", &filter) == INPARSE_ERR) status = true;
    if(filter){
        if(Inpar->getPar("f0", &freqs[0]) == INPARSE_ERR) status = true;
        if(Inpar->getPar("f1", &freqs[1]) == INPARSE_ERR) status = true;
        if(Inpar->getPar("f2", &freqs[2]) == INPARSE_ERR) status = true;
        if(Inpar->getPar("f3", &freqs[3]) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("zder", &zder) == INPARSE_ERR) status = true;

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    lsmig->setOrder(order);
    lsmig->setLpml(lpml);
    lsmig->setFs(fs);
    lsmig->setSnapinc(snapinc);

    lsmig->setPrecordfile(Precordfile);
    lsmig->setDataweight(dataweight);
    lsmig->setDataweightfile(Dataweightfile);
    if(mute){
        lsmig->setMutefile(Mutefile);
    }
    lsmig->setVpfile(VP0FILE);
    lsmig->setRhofile(RHO0FILE);
    lsmig->setWaveletfile(SOURCE0FILE);
    lsmig->setMisfitfile(MISFITFILE);
    lsmig->setPmodelledfile(PMODFILE);
    lsmig->setPresidualfile(PRESFILE);
    lsmig->setPsnapfile(Psnapfile);
    lsmig->setApertx(apertx);
    lsmig->setSnapmethod(snapmethod);
    lsmig->setNsnaps(nsnaps);
    lsmig->setIncore(incore);
    lsmig->setMisfit_type(fwimisfit);

    lsmig->setVpgradfile(VPGRADFILE);
    lsmig->setKvp(kvp);
    lsmig->setVpregalpha(vpregalpha);
    lsmig->setPimagefile(VPLSFILE);

    lsmig->setUpdate_vp(true);
    lsmig->setSrcilum(srcilum);
    lsmig->setSrcilumfile(SRCILUMFILE);
    lsmig->setFilter(filter);
    lsmig->setFreqs(&freqs[0]);
    lsmig->setZder(zder);

    //MASTER
    if(mpi.getRank() == 0){
        // Create a sort class and map over shots
        std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
        Sort->setDatafile(Precordfile);
        Sort->createShotmap(Precordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // L-BFGS configuration
        double *x = nullptr; 
        int N;
        N = lsmig->setInitial(x, Vpfile, Rhofile, Waveletfile);
        x = (double *) calloc(N, sizeof(double));
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));
        opt->opt_set_initial_guess(x);
        opt->setMax_linesearch(max_linesearch);
        opt->setMax_iterations(max_iterations);

        switch(optmethod) {
            case 1:
                opt->setGtol(0.9);
                break;
            case 4:
                opt->setGtol(0.1);
            default:
                break;
        }

        switch(linesearch) {
            case 1:
                opt->setLinesearch_condition(OPT_CONDITION_DECREASE);
                break;
            case 2:
                opt->setLinesearch_condition(OPT_CONDITION_ARMIJO);
                break;
            case 4:
                opt->setLinesearch_condition(OPT_CONDITION_STRONG_WOLFE);
                break;
            case 3:
            default:
                opt->setLinesearch_condition(OPT_CONDITION_WOLFE);
                break;
        }

        // Create results folder
        lsmig->createResult();

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        lsmig->writeLog("Starting optimisation algorithm");
        lsmig->writeProgress("Maximum number of iterations: " + std::to_string(opt->getMax_iterations()));
        lsmig->writeProgress("Maximum number of linesearches: " + std::to_string(opt->getMax_linesearch()));
        lsmig->writeProgress("ITERATION\t   STEP LENGTH            MISFIT               GNORM             MNORM                TIME");


        switch(optmethod) {
            case 4:
                opt->opt_conjugate_gradient_pr(evaluate,progress);
                break;
            case 3:
                opt->opt_steepest_descent(evaluate,progress);
                break;
            case 2:
                opt->opt_conjugate_gradient_fr(evaluate,progress);
                break;
            case 1:
            default:
                opt->opt_lbfgs(evaluate,progress,finalize);
                break;
        }

        // Send message for slaves to quit
        lsmig->writeLog("Optimisation algorithm finished");
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
        lsmig->writeProgress(buffer);

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
                    lsmig->runGrad();
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

