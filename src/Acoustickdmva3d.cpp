#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "kdmva.h"
#include "inparse.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<KdmvaAcoustic3D<float>> kdmva;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    kdmva->writeLog("##### Starting new evaluation #####");
    double *x = instance->x;
    double *g = instance->g;

    kdmva->writeLog("Saving linesearch models");
    // Save linesearch model
    kdmva->saveLinesearch(x);
    kdmva->writeLog("Linesearch models saved");

    // Start new gradient evaluation
    kdmva->writeLog("Starting gradient computation");
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    kdmva->runGrad();
    kdmva->writeLog("Gradient computation finished");

    // Compute regularization
    kdmva->writeLog("Computing regularisation");
    kdmva->computeRegularisation(x);

    // Combine data and model misfit gradients
    kdmva->writeLog("Combining gradients");
    kdmva->combineGradients();

    // Apply Chain rule of logistic model
    if(kdmva->getConstrain()){
       kdmva->writeLog("Computing Chain rule");
       kdmva->applyChainrule(x);
    }

    // Apply mute to gradient
    kdmva->writeLog("Muting gradients");
    kdmva->applyMute();

    // Project gradient to B-spline 
    if(kdmva->getParamtype() == PAR_BSPLINE || kdmva->getParamtype() == PAR_AVG)
    {
        kdmva->writeLog("Projecting gradient in B-spline grid");
        task = RUN_BS_PROJ;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
        kdmva->runBsproj();
    }

    //Read final gradient 
    kdmva->writeLog("Reading gradient into vector");
    kdmva->readGrad(g);

    //Read misfit
    kdmva->writeLog("Reading misfits");
    kdmva->readMisfit(&instance->f);
    kdmva->writeLog("##### Evaluation finished #####");

    // Normalize error and gradient
    if(kdmva->getFnorm() == 0.0){
        kdmva->setFnorm(instance->f);
    }
    kdmva->normalize(g, &instance->f, instance->n);

    double xnorm, gnorm, step;
    gnorm = kdmva->vector_norm(instance->g, 2, instance->n);
    xnorm = kdmva->vector_norm(instance->x, 2, instance->n);
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
	kdmva->writeProgress(buffer);
}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    // Copy new iteration files to results folder
    kdmva->saveResults(opt->getIter());

    double xnorm, gnorm, step;
    gnorm = kdmva->vector_norm(instance->g, 2, instance->n);
    xnorm = kdmva->vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    kdmva->writeLog("---------- > New iteration found: " + std::to_string(opt->getIter()));

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
	kdmva->writeProgress(buffer);
}

void finalize(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{
    /*
    if(opt->getCompdiaghessian())
    {
        kdmva->writeLog("Saving diagonal Hessian");
        double *x = instance->diaghessian;
        kdmva->un_normalize(x, instance->f, instance->n);
        kdmva->saveHessian(x);
    }else{
        // Do nothing
    }
    */
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
            PRINT_DOC(#Mva parameters);
            PRINT_DOC(misfit_type = "0";  # 0- SI; 1- DS; 2- DS_PLUS_SI;);
            PRINT_DOC(modelmute = "false";  # Mute gradient and updates);
            PRINT_DOC(Modelmutefile = "Mmute.rss"; # File with mute weights);
            PRINT_DOC(residualmute = "false";  # Mute residual image);
            PRINT_DOC(Residualmutefile = "Rmute.rss"; # File with mute weights);
            PRINT_DOC(max_linesearch = "5"; # maximum number of linesearches);
            PRINT_DOC(max_iterations = "20"; # maximum number of iterations);
            PRINT_DOC(constrain = "false";  # Constrain inversion requires Lboundfile and Uboundfile containing models with the bounds);
            PRINT_DOC(incore = "false";  # Do all traveltime computation incore (No traveltime table storage));

            PRINT_DOC(optmethod = "1"; # 1-L-BFGS; 2-CG_FR; 3-STEEPEST DESCENT; 4-CG_PR);
            PRINT_DOC(linesearch = "3"; # 1-Decrease; 2-Armijo; 3-Wolfe; 4-Strong Wolfe);
            PRINT_DOC(update_vp = "true"; # Update vp);
            PRINT_DOC();
            PRINT_DOC(# Diagonal scaling parameters);
            PRINT_DOC(kvp = "100.0";);
            PRINT_DOC();
            PRINT_DOC(#Migration parameters);
            PRINT_DOC(apertx = "0"; # Aperture for local model (source is in the middle, 0 means automatically set to include all receivers));
            PRINT_DOC(aperty = "0"; # Aperture for local model (source is in the middle, 0 means automatically set to include all receivers));
            PRINT_DOC(Souinc = "1"; # Source decimation factor in traveltime table);
            PRINT_DOC(Recinc = "1"; # Receiver decimation factor in traveltime table);
            PRINT_DOC(radius = "50.0"; # Radius of traveltime interpolation);
            PRINT_DOC(nhx = "1"; # Subsurface offsets in x direction);
            PRINT_DOC(nhy = "1"; # Subsurface offsets in y direction);
            PRINT_DOC(nhz = "1"; # Subsurface offsets in z direction);
            PRINT_DOC();
            PRINT_DOC(#Parameterisation);
            PRINT_DOC(paramtype = "1";  # 0- grid; 1- B-spline; 3- Average);
            PRINT_DOC(dtx = "25.0"; # knot sampling in B-spline);
            PRINT_DOC(dty = "25.0"; # knot sampling in B-spline);
            PRINT_DOC(dtz = "25.0"; # knot sampling in B-spline);
            PRINT_DOC();
            PRINT_DOC(#Regularisation);
            PRINT_DOC(vpregalpha = "0.0";);
            PRINT_DOC();
            PRINT_DOC(# Uncertainty);
            PRINT_DOC(outputhess = "false"; # Output diagonal of L-BFGS inverse Hessian at last iteration ;)
            PRINT_DOC();
            PRINT_DOC(# Input files);
            PRINT_DOC(Vp = "Vp2d.rss";);
            PRINT_DOC(Precordfile = "Pdata.rss";);
            PRINT_DOC(Lbound = "MinVp2d.rss";);
            PRINT_DOC(Ubound = "MaxVp2d.rss";);
        }
        exit(1);
    }

    // Initialize Inversion class
    kdmva = std::make_shared<rockseis::KdmvaAcoustic3D<float>>(&mpi);
	/* General input parameters */
    bool status;
    bool modelmute;
    bool residualmute;
    bool constrain;
    bool outputhess;
    int _paramtype;
    int misfit_type;
    float dtx=-1;
    float dty=-1;
    float dtz=-1;
    float radius;
    int nhx, nhy, nhz;
    int souinc, recinc;
    float apertx;
    float aperty;
    float kvp;
    float vpregalpha;
    int max_linesearch, max_iterations;
    int linesearch;
    int optmethod; 
    bool incore;
    std::string Vpfile;
    std::string Misfitfile;
    std::string Precordfile;
    std::string Lboundfile;
    std::string Uboundfile;
    std::string Modelmutefile;
    std::string Residualmutefile;


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
    if(Inpar->getPar("outputhess", &outputhess) == INPARSE_ERR) status = true;
    rockseis::rs_paramtype paramtype = static_cast<rockseis::rs_paramtype>(_paramtype);
    if(paramtype == PAR_BSPLINE || paramtype == PAR_AVG){
        if(Inpar->getPar("dtx", &dtx) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dty", &dty) == INPARSE_ERR) status = true;
        if(Inpar->getPar("dtz", &dtz) == INPARSE_ERR) status = true;
    }
    if(Inpar->getPar("nhx", &nhx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhy", &nhy) == INPARSE_ERR) status = true;
    if(Inpar->getPar("nhz", &nhz) == INPARSE_ERR) status = true;
    if(Inpar->getPar("radius", &radius) == INPARSE_ERR) status = true;
    if(Inpar->getPar("incore", &incore) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Souinc", &souinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Recinc", &recinc) == INPARSE_ERR) status = true;
    if(Inpar->getPar("apertx", &apertx) == INPARSE_ERR) status = true;
    if(Inpar->getPar("aperty", &aperty) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;

    if(Inpar->getPar("modelmute", &modelmute) == INPARSE_ERR) status = true;
    if(modelmute){
        if(Inpar->getPar("Modelmutefile", &Modelmutefile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("residualmute", &residualmute) == INPARSE_ERR) status = true;
    if(residualmute){
        if(Inpar->getPar("Residualmutefile", &Residualmutefile) == INPARSE_ERR) status = true;
    }

    if(Inpar->getPar("vpregalpha", &vpregalpha) == INPARSE_ERR) status = true;
    if(Inpar->getPar("max_linesearch", &max_linesearch) == INPARSE_ERR) status = true;
    if(Inpar->getPar("max_iterations", &max_iterations) == INPARSE_ERR) status = true;

    if(Inpar->getPar("misfit_type", &misfit_type) == INPARSE_ERR) status = true;
    rockseis::rs_wemvamisfit wemvamisfit = static_cast<rockseis::rs_wemvamisfit>(misfit_type);

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

    kdmva->setPrecordfile(Precordfile);
    if(modelmute){
        kdmva->setModelmutefile(Modelmutefile);
    }
    if(residualmute){
        kdmva->setResidualmutefile(Residualmutefile);
    }
    kdmva->setConstrain(constrain);
    if(constrain){
        kdmva->setLboundfile(Lboundfile);
        kdmva->setUboundfile(Uboundfile);
    }
    kdmva->setVpfile(VPLSFILE);
    kdmva->setTtablefile(TTABLELSFILE);
    kdmva->setMisfitfile(MISFITFILE);
    kdmva->setApertx(apertx);
    kdmva->setAperty(aperty);

    if(wemvamisfit == DS_HMAX && nhz > 1) rs_error("DS_HMAX misfit is currently only implemented for the cases where nhz = 1");
    kdmva->setMisfit_type(wemvamisfit);

    kdmva->setVpgradfile(VPGRADFILE);
    kdmva->setPimagefile(PIMAGEFILE);
    kdmva->setKvp(kvp);
    kdmva->setParamtype(paramtype);

    kdmva->setNhx(nhx);
    kdmva->setNhy(nhy);
    kdmva->setNhz(nhz);

    kdmva->setDtx(dtx);
    kdmva->setDty(dty);
    kdmva->setDtz(dtz);

    kdmva->setVpregalpha(vpregalpha);

    kdmva->setRadius(radius);
    kdmva->setIncore(incore);
    kdmva->setSouinc(souinc);
    kdmva->setRecinc(recinc);

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
        N = kdmva->setInitial(x, Vpfile);
        x = (double *) calloc(N, sizeof(double));
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));
        opt->opt_set_initial_guess(x);
        opt->setMax_linesearch(max_linesearch);
        opt->setMax_iterations(max_iterations);
        opt->setCompdiaghessian(outputhess);

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
        kdmva->createResult();

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        kdmva->writeLog("Starting optimisation algorithm");
        kdmva->writeProgress("Maximum number of iterations: " + std::to_string(opt->getMax_iterations()));
        kdmva->writeProgress("Maximum number of linesearches: " + std::to_string(opt->getMax_linesearch()));
        kdmva->writeProgress("ITERATION\t   STEP LENGTH            MISFIT               GNORM             MNORM                TIME");


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
        kdmva->writeLog("Optimisation algorithm finished");
        task = BREAK_LOOP;
        MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Report the result. */
        char buffer[512];
        snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
        kdmva->writeProgress(buffer);

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
                    kdmva->runGrad();
                    break;
                case RUN_BS_PROJ:
                    kdmva->runBsproj();
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

