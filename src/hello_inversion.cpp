#include <iostream>
#include <mpi.h>
#include "inversion.h"
#include "opt.h"
#define K 400

using namespace rockseis;

/// TO DO: 
//Add data weight and mute 
//Add functions to save line search files 
//Add B-spline possibilites 
//Make parameter loading possible
/* Global variables */
//std::shared_ptr<MPImodeling> mpi;
MPImodeling *mpi;
std::shared_ptr<Inversion<float>> inv;
std::string vplsfile = "vp_ls.rss";
std::string rholsfile = "rho_ls.rss";
std::string sourcelsfile = "source_ls.rss";

std::string vp0file = "vp_0.rss";
std::string rho0file = "rho_0.rss";
std::string source0file = "source_0.rss";

std::string vpgradfile = "vp_grad.rss";
std::string rhogradfile = "rho_grad.rss";
std::string sourcegradfile = "source_grad.rss";

std::string misfitfile = "misfit.rss";

/* Global functions */
void evaluate(rockseis::OptInstancePtr instance)
{
    fprintf(stderr, "Starting new evaluation\n");
    double *x = instance->x;
    double *g = instance->g;
    // Models
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> model0 (new rockseis::ModelAcoustic2D<float>(vp0file, rho0file, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> source0 (new rockseis::Data2D<float>(source0file));
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> lsmodel (new rockseis::ModelAcoustic2D<float>(vplsfile, rholsfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> lssource (new rockseis::Data2D<float>(sourcelsfile));

    // Write linesearch model
    model0->readModel();
    lsmodel->readModel();
    float *vp0, *rho0, *vpls, *rhols;
    vp0 = model0->getVp(); 
    rho0 = model0->getR(); 
    vpls = lsmodel->getVp(); 
    rhols = lsmodel->getR(); 
    int i;
    int N = (lsmodel->getGeom())->getNtot();
    float k=K;
    std::cerr << "In evaluate: N: " << N << std::endl;
    for(i=0; i< N; i++)
    {
        vpls[i] = vp0[i] + x[i]*k;
        rhols[i] = rho0[i] + x[N+i]*k;
    }
    lsmodel->writeModel();


    // Start new gradient evaluation
    int task;
    task = RUN_F_GRAD;
    std::cerr << "Rank: " << mpi->getRank() << ": Giving order to run gradient." << std::endl; 
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv->runAcousticfwigrad2d(mpi);

    //Read gradient
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> modelgrad (new rockseis::ModelAcoustic2D<float>(vpgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> sourcegrad (new rockseis::Data2D<float>(sourcegradfile));
    modelgrad->readModel();
    vp0 = modelgrad->getVp(); 
    rho0 = modelgrad->getR(); 
    for(i=0; i< N; i++)
    {
        g[i] = vp0[i]*k;
        //g[N+i] = rho0[i]*k;
        g[N+i] = 0.0;
    }

    //Read misfit
    instance->f = 0.0;
    std::shared_ptr<rockseis::File> Fmisfit (new rockseis::File());
    Fmisfit->input(misfitfile);
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
    //mpi = std::make_shared<rockseis::MPImodeling>(&argc,&argv);
    mpi =  new rockseis::MPImodeling(&argc,&argv);
    int task; 

    inv = std::make_shared<rockseis::Inversion<float>>();

    //MASTER
    if(mpi->getRank() == 0){
        // Get input model
        std::shared_ptr<rockseis::ModelAcoustic2D<float>> model0 (new rockseis::ModelAcoustic2D<float>(vp0file, rho0file, 1 ,0));
        std::shared_ptr<rockseis::Data2D<float>> source0 (new rockseis::Data2D<float>(source0file));

        // L-BFGS parameters
        int N=(model0->getGeom())->getNtot();
        std::cerr << "In main: N: " << N << std::endl;
        std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(2*N));

        /* Initialize the parameters for the L-BFGS optimization. */

        /* Initialize x */
        double *x = (double *) calloc(2*N, sizeof(double));
        //float *vp, *rho;
        model0->readModel();
        //vp = model0->getVp(); 
        //rho = model0->getR(); 
        //int i;
        //float k=K;
        //for(i=0; i< N; i++)
        //{
         //   x[i] = vp[i]/k;
         //   x[N+i] = rho[i]/k;
        //}
        opt->opt_set_initial_guess(x);
        // Copy initial model to linesearch model
        model0->setVpfile(vplsfile);
        model0->setRfile(rholsfile);
        model0->writeModel();
        source0->setFile(sourcelsfile);
        source0->write();
        opt->setGtol(0.9);
        opt->setMax_linesearch(5);

        // Start the L-BFGS optimization; this will invoke the callback functions evaluate() and progress() when necessary.
        //
        opt->opt_lbfgs(evaluate, progress);

        // Send message for slaves to quit
        task = BREAK_LOOP;
        std::cerr << "Rank: " << mpi->getRank() << ": Giving order to break the loop." << std::endl; 
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
                    std::cerr << "Rank: " << mpi->getRank() << ": Running gradient." << std::endl; 
                    inv->runAcousticfwigrad2d(mpi);
                    break;
                case BREAK_LOOP:
                    std::cerr << "Rank: " << mpi->getRank() << ": Leaving the loop." << std::endl; 
                    stop = true;
                    break;
            }
            if(stop){
                break;
            }
        }

    }

    delete mpi;
    return 0;
}

