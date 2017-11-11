#include <iostream>
#include <mpi.h>
#include "opt.h"
#include "inversion.h"
#include "inparse.h"
#define K 400

using namespace rockseis;

/// TO DO: 
//Add mute 
//Add functions to save line search files 
//Add B-spline possibilites 
//Make parameter loading possible
//
/* Global variables */
std::shared_ptr<InversionAcoustic2D<float>> inv;

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    fprintf(stderr, "Starting new evaluation\n");
    double *x = instance->x;
    double *g = instance->g;
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> model0 (new rockseis::ModelAcoustic2D<float>(VP0FILE, RHO0FILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> source0 (new rockseis::Data2D<float>(SOURCE0FILE));
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> lsmodel (new rockseis::ModelAcoustic2D<float>(VPLSFILE, RHOLSFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> lssource (new rockseis::Data2D<float>(SOURCELSFILE));

    // Write linesearch model
    model0->readModel();
    lsmodel->readModel();
    float *vp0, *rho0, *vpls, *rhols;
    vp0 = model0->getVp(); 
    rho0 = model0->getR(); 
    vpls = lsmodel->getVp(); 
    rhols = lsmodel->getR(); 
    int i;
    float k=K;

    int N = (lsmodel->getGeom())->getNtot();

    for(i=0; i< N; i++)
    {
        vpls[i] = vp0[i] + x[i]*k;
        rhols[i] = rho0[i] + x[N+i]*k;
    }
    lsmodel->writeModel();


    // Start new gradient evaluation
    int task;
    task = RUN_F_GRAD;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv->runAcousticfwigrad2d();

    task = RUN_BS_PROJ;
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv->runBsprojection2d();


    //Read gradient
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> modelgrad (new rockseis::ModelAcoustic2D<float>(VPGRADFILE, RHOGRADFILE, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> sourcegrad (new rockseis::Data2D<float>(SOURCEGRADFILE));
    modelgrad->readModel();
    vp0 = modelgrad->getVp(); 
    rho0 = modelgrad->getR(); 
    for(i=0; i< N; i++)
    {
        g[i] = vp0[i]*k;
        g[N+i] = rho0[i]*k;
    }

    //Read misfit
    instance->f = 0.0;
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->input(MISFITFILE);
    float val;
    for(i=0; i<Fmisfit->getN(1); i++){
        Fmisfit->read(&val, 1); 
        instance->f += val;
    }
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
	int order;
	int snapinc;
	int nsnaps = 0;
	int snapmethod;
	int misfit_type;
    float apertx;
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
    if(Inpar->getPar("Psnapfile", &Psnapfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Vpgradfile", &Vpgradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Rhogradfile", &Rhogradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Wavgradfile", &Wavgradfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Misfitfile", &Misfitfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Precordfile", &Precordfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Presidualfile", &Presidualfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("Pmodelledfile", &Pmodelledfile) == INPARSE_ERR) status = true;
    if(Inpar->getPar("snapmethod", &snapmethod) == INPARSE_ERR) status = true;
    rockseis::rs_snapmethod checkpoint = static_cast<rockseis::rs_snapmethod>(snapmethod);
    switch(checkpoint){
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

	if(status == true){
		rs_error("Program terminated due to input errors.");
	}

    inv->setOrder(order);
    inv->setLpml(lpml);
    inv->setFs(fs);
    inv->setSnapinc(snapinc);

    inv->setDataweight(dataweight);
    inv->setDataweightfile(Dataweightfile);
    inv->setVpfile(Vpfile);
    inv->setRhofile(Rhofile);
    inv->setWaveletfile(Waveletfile);
    inv->setMisfitfile(Misfitfile);
    inv->setPsnapfile(Psnapfile);
    inv->setPrecordfile(Precordfile);
    inv->setPmodelledfile(Pmodelledfile);
    inv->setPresidualfile(Presidualfile);
    inv->setApertx(apertx);
    inv->setSnapmethod(checkpoint);
    inv->setNsnaps(nsnaps);
    inv->setIncore(incore);
    inv->setMisfit_type(fwimisfit);

    inv->setVpgradfile(Vpgradfile);
    inv->setRhogradfile(Rhogradfile);
    inv->setWavgradfile(Wavgradfile);

    float dtx = 50.0;
    float dtz = 50.0;

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

        // Get input model
        std::shared_ptr<rockseis::ModelAcoustic2D<float>> model0 (new rockseis::ModelAcoustic2D<float>(VP0FILE, RHO0FILE, 1 ,0));
        std::shared_ptr<rockseis::Data2D<float>> source0 (new rockseis::Data2D<float>(SOURCE0FILE));

        // L-BFGS parameters
        int N=(model0->getGeom())->getNtot();
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(2*N));

        /* Initialize the parameters for the L-BFGS optimization. */

        /* Initialize x */
        double *x = (double *) calloc(2*N, sizeof(double));
        //float *vp, *rho;
        opt->opt_set_initial_guess(x);
        // Copy initial model to linesearch model
        model0->readModel();
        model0->setVpfile(VPLSFILE);
        model0->setRfile(RHOLSFILE);
        model0->writeModel();
        source0->setFile(SOURCELSFILE);
        source0->write();
        opt->setGtol(0.9);
        opt->setMax_linesearch(1);

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

    //SLAVE
    }else{
        bool stop = false;
        while(1)
        {
            MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
            switch(task)
            {
                case RUN_F_GRAD:
                    inv->runAcousticfwigrad2d();
                    break;
                case RUN_BS_PROJ:
                    inv->runBsprojection2d();
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

