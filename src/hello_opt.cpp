// Header inclusion 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <memory>
#include "opt.h"

/* Constants */
#define PI 3.14159265358979323846

/* Global variables */
void evaluate(rockseis::OptInstancePtr instance)
{
    fprintf(stderr, "Starting new evaluation\n");
    double *f = &instance->f;
    double *x = instance->x;
    double *g = instance->g;
	instance->f = 0.0;
    
	*f = ((x[0]*300)-100)*((x[0]*300)-100);
    g[0] = 2.0*((x[0]*300)-100)*300;
    fprintf(stderr, "  fx = %f, x[0] = %f\n", instance->f, instance->x[0]);
}

void progress(rockseis::Opt *opt, rockseis::OptInstancePtr instance)
{

    double xnorm, gnorm, step;
    gnorm = opt->opt_vector_norm(instance->g, 2, instance->n);
    xnorm = opt->opt_vector_norm(instance->x, 2, instance->n);
    step = instance->steplength;
    fprintf(stderr, "Iteration %d:\n", opt->getIter());
    fprintf(stderr, "  fx = %f, x[0] = %f\n", instance->f, instance->x[0]);
    fprintf(stderr, "  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
}

// Main 
int main(int argc, char* argv[]) {
	/* L-BFGS parameters */
    int N=1;

    std::shared_ptr<rockseis::Opt> opt (new rockseis::Opt(N));

	/* Initialize the parameters for the L-BFGS optimization. */

    /* Initialize x */
    double x[1];
    x[0] = 1;
    opt->opt_set_initial_guess(&x[0]);

    opt->setGtol(0.9);

    /*
       Start the L-BFGS optimization; this will invoke the callback functions
       evaluate() and progress() when necessary.
       */
    //opt->opt_steepest_descent(evaluate, progress);
    //opt->opt_conjugate_gradient_pr(evaluate, progress);
    opt->opt_lbfgs(evaluate, progress);
    //opt->opt_heat_bath(evaluate, progress);

    /* Report the result. */
    char buffer[512];
    snprintf(buffer,512,"L-BFGS algorithm finished with return status:\n%s\n",opt->getMsg());
    fprintf(stderr, "%s", buffer);

	//exit(0);
}
