#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "tomo.h"
#include "inparse.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<TomoAcoustic2D<float>> fatt;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    fatt->writeLog("##### Starting new evaluation #####");
    double *x = instance->x;
    double *g = instance->g;

    fatt->writeLog("Saving linesearch models");
    // Save linesearch model
    fatt->saveLinesearch(x);
    fatt->writeLog("Linesearch models saved");

    // Start new gradient evaluation
    fatt->writeLog("Starting gradient computation");
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    fatt->runGrad();
    fatt->writeLog("Gradient computation finished");

    // Compute regularization
    fatt->writeLog("Computing regularisation");
    fatt->computeRegularisation(x);

    // Combine data and model misfit gradients
    fatt->writeLog("Combining gradients");
    fatt->combineGradients();

    // Apply mute to gradient
    fatt->writeLog("Muting gradients");
    fatt->applyMute();

    // Project gradient to B-spline 
    if(fatt->getParamtype() == PAR_BSPLINE)
    {
        fatt->writeLog("Projecting gradient in B-spline grid");
        task = RUN_BS_PROJ;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
        fatt->runBsproj();
    }

    //Read final gradient 
    fatt->writeLog("Reading gradient into vector");
    fatt->readGrad(g);

    //Read misfit
    fatt->writeLog("Reading misfits");
    fatt->readMisfit(&instance->f);
    fatt->writeLog("##### Evaluation finished #####");

    // Normalize error and gradient
    if(fatt->getFnorm() == 0.0){
        fatt->setFnorm(instance->f);
    }
    fatt->normalize(g, &instance->f, instance->n);

    double xnorm, gnorm, step;
    gnorm = fatt->vector_norm(instance->g, 2, instance->n);
    xnorm = fatt->vector_norm(instance->x, 2, instance->n);
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
	fatt->writeProgress(buffer);
}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    // Copy new iteration files to results folder
    fatt->saveResults(opt->getIter());

    double xnorm, gnorm, step;
    gnorm = fatt->vector_norm(instance->g, 2, instance->n);
    xnorm = fatt->vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    fatt->writeLog("---------- > New iteration found: " + std::to_string(opt->getIter()));

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
	fatt->writeProgress(buffer);
}

void finalize(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    if(opt->getCompdiaghessian())
    {
        fatt->writeLog("Saving diagonal Hessian");
        double *x = instance->diaghessian;
        fatt->un_normalize(x, instance->f, instance->n);
        fatt->saveHessian(x);
    }else{
        // Do nothing
    }
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
            PRINT_DOC(# MPI 2d first arrival tomography configuration file);
            PRINT_DOC();
            PRINT_DOC(#Fatt parameters);
            PRINT_DOC(dataweight = "false";);
            PRINT_DOC(Dataweightfile = "weights.rss";);
            PRINT_DOC(mute = "false";  # Mute gradient and updates);
            PRINT_DOC(Mutefile = "mute.rss"; # File with mute weights);
            PRINT_DOC(max_linesearch = "5"; # maximum number of linesearches);
            PRINT_DOC(max_iterations = "20"; # maximum number of iterations);
            PRINT_DOC(constrain = "false";  # Constrain inversion requires Lboundfile and Uboundfile containing models with the bounds);

            PRINT_DOC(optmethod = "1"; # 1-L-BFGS; 2-CG_FR; 3-STEEPEST DESCENT; 4-CG_PR);
            PRINT_DOC(linesearch = "3"; # 1-Decrease; 2-Armijo; 3-Wolfe; 4-Strong Wolfe);
            PRINT_DOC(update_vp = "true"; # Update vp);
            PRINT_DOC();
            PRINT_DOC(# Diagonal scaling parameters);
            PRINT_DOC(kvp = "100.0";);
            PRINT_DOC();
            PRINT_DOC(#Parameterisation);
            PRINT_DOC(paramtype = "1";  # 0- grid; 1- B-spline;);
            PRINT_DOC(dtx = "25.0"; # knot sampling in B-spline);
            PRINT_DOC(dtz = "25.0"; # knot sampling in B-spline);
            PRINT_DOC();
            PRINT_DOC(#Regularisation);
            PRINT_DOC(vpregalpha = "0.0";);
            PRINT_DOC();
            PRINT_DOC(# Files);
            PRINT_DOC(Vp = "Vp2d.rss";);
            PRINT_DOC(Trecordfile = "Tdata.rss";);
            PRINT_DOC(Lbound = "MinVp2d.rss";);
            PRINT_DOC(Ubound = "MaxVp2d.rss";);
        }
        exit(1);
    }

    // Initialize Inversion class
    fatt = std::make_shared<rockseis::TomoAcoustic2D<float>>(&mpi);
	/* General input parameters */
    bool status;
    bool dataweight;
    bool mute;
    bool constrain;
    int _paramtype;
    float dtx=-1;
    float dtz=-1;
    float kvp;
    float vpregalpha;
    int max_linesearch, max_iterations;
    int linesearch;
    int optmethod; 
    std::string Vpfile;
    std::string Vpgradfile;
    std::string Dataweightfile;
    std::string Misfitfile;
    std::string Trecordfile;
    std::string Tmodelledfile;
    std::string Tresidualfile;
    std::string Mutefile;
    std::string Lboundfile;
    std::string Uboundfile;


    /* Get parameters from configuration file */
    std::shared_ptr<rockseis::Inparse> Inpar (new rockseis::Inparse());
    if(Inpar->parse(argv[1]) == INPARSE_ERR) 
    {
        rs_error("Parse error on input config file", argv[1]);
    }
    status = false; 
    if(Inpar->getPar("Vp", &Vpfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("kvp", &kvp) == INPARSE_ERR) status = true;
    if(Inpar->getPar("paramtype", &_paramtype) == INPARSE_ERR) status = true;
    rockseis::rs_paramtype paramtype = static_cast<rockseis::rs_paramtype>(_paramtype);
    if(paramtype == PAR_BSPLINE){
        if(Inpar->getPar("dtx", &dtx) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dtz", &dtz) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("Trecordfile", &Trecordfile) == INPARSE_ERR) status = true;

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

    if(Inpar->getPar("linesearch", &linesearch) == INPARSE_ERR) status = true;
    if(Inpar->getPar("optmethod", &optmethod) == INPARSE_ERR) status = true;

    if(Inpar->getPar("constrain", &constrain) == INPARSE_ERR) status = true;
    if(constrain){
        if(Inpar->getPar("Lbound", &Lboundfile) == INPARSE_ERR) status = true;
        if(Inpar->getPar("Ubound", &Uboundfile) == INPARSE_ERR) status = true;
    }

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    fatt->setTrecordfile(Trecordfile);
    fatt->setDataweight(dataweight);
    fatt->setDataweightfile(Dataweightfile);
    if(mute){
        fatt->setMutefile(Mutefile);
    }
    fatt->setConstrain(constrain);
    if(constrain){
        fatt->setLboundfile(Lboundfile);
        fatt->setUboundfile(Uboundfile);
    }
    fatt->setVpfile(VPLSFILE);
    fatt->setMisfitfile(MISFITFILE);
    fatt->setTmodelledfile(TMODFILE);
    fatt->setTresidualfile(TRESFILE);

    fatt->setVpgradfile(VPGRADFILE);
    fatt->setKvp(kvp);
    fatt->setParamtype(paramtype);
    fatt->setDtx(dtx);
    fatt->setDtz(dtz);

    fatt->setVpregalpha(vpregalpha);

    //MASTER
    if(mpi.getRank() == 0){
        // Create a sort class and map over shots
        std::shared_ptr<rockseis::Sort<float>> Sort (new rockseis::Sort<float>());
        Sort->setDatafile(Trecordfile);
        Sort->createShotmap(Trecordfile); 
        Sort->writeKeymap();
        Sort->writeSortmap();

        // L-BFGS configuration
        double *x = nullptr; 
        int N;
        N = fatt->setInitial(x, Vpfile);
        x = (double *) calloc(N, sizeof(double));
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));
        opt->opt_set_initial_guess(x);
        opt->setMax_linesearch(max_linesearch);
        opt->setMax_iterations(max_iterations);
        opt->setCompdiaghessian(true);

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
        fatt->createResult();

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        fatt->writeLog("Starting optimisation algorithm");
        fatt->writeProgress("Maximum number of iterations: " + std::to_string(opt->getMax_iterations()));
        fatt->writeProgress("Maximum number of linesearches: " + std::to_string(opt->getMax_linesearch()));
        fatt->writeProgress("ITERATION\t   STEP LENGTH            MISFIT               GNORM             MNORM                TIME");


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
        fatt->writeLog("Optimisation algorithm finished");
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
        fatt->writeProgress(buffer);

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
                    fatt->runGrad();
                    break;
                case RUN_BS_PROJ:
                    fatt->runBsproj();
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

