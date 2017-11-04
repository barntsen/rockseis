#include <iostream>
#include <mpi.h>
#include "inversion.h"
#include "opt.h"

using namespace rockseis;

/* Global variables */
std::shared_ptr<MPImodeling> mpi;
Inversion<float> inv = Inversion<float>();
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
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> lsmodel (new rockseis::ModelAcoustic2D<float>(vplsfile, rholsfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> lssource (new rockseis::Data2D<float>(sourcelsfile));

    // Write linesearch model
    lsmodel->readModel();
    float *vp, *rho;
    vp = lsmodel->getVp(); 
    rho = lsmodel->getR(); 
    int i;
    int N = (lsmodel->getGeom())->getNtot();
    std::cerr << "In evaluate: N: " << N << std::endl;
    for(i=0; i< N; i++)
    {
        vp[i] = x[i];
        rho[i] = x[N+i];
    }
    lsmodel->writeModel();


    // Start new gradient evaluation
    int task;
    task = RUN_F_GRAD;
    std::cerr << "Rank: " << mpi->getRank() << ": Giving order to run gradient." << std::endl; 
    MPI_Bcast(&task, 1, MPI_INT, 0, MPI_COMM_WORLD);
    inv.runAcousticfwigrad2d(mpi);

    //Read gradient
    std::shared_ptr<rockseis::ModelAcoustic2D<float>> modelgrad (new rockseis::ModelAcoustic2D<float>(vpgradfile, rhogradfile, 1 ,0));
    std::shared_ptr<rockseis::Data2D<float>> sourcegrad (new rockseis::Data2D<float>(sourcegradfile));
    modelgrad->readModel();
    vp = modelgrad->getVp(); 
    rho = modelgrad->getR(); 
    for(i=0; i< N; i++)
    {
        g[i] = vp[i];
        g[N+i] = rho[i];
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
    mpi = std::make_shared<rockseis::MPImodeling>(&argc,&argv);
    int task; 

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
        float *vp, *rho;
        model0->readModel();
        vp = model0->getVp(); 
        rho = model0->getR(); 
        int i;
        for(i=0; i< N; i++)
        {
            x[i] = vp[i];
            x[N+i] = rho[i];
        }
        opt->opt_set_initial_guess(x);
        // Copy initial model to linesearch model
        model0->setVpfile(vplsfile);
        model0->setRfile(rholsfile);
        model0->writeModel();
        source0->setFile(sourcelsfile);
        source0->write();
        opt->setGtol(0.9);
        opt->setMax_linesearch(1);

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
                    inv.runAcousticfwigrad2d(mpi);
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

        return 0;
    }
}

