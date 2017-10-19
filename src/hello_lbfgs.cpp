// Header inclusion 
#include <stdio.h>
#include <stdlib.h>
#include <lbfgs.h>

/* Constants */
#define PI 3.14159265358979323846

/* Global variables */
static lbfgsfloatval_t evaluate(
		void *instance,
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		)
{
    fprintf(stderr, "Starting new evaluation\n");
	lbfgsfloatval_t fx = 0.0;
	fx = ((x[0]*300)-100)*((x[0]*300)-100);
    g[0] = 2.0*((x[0]*300)-100)*300;

	return (fx);
}

static int progress(
		void *instance,
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
		)
{
    fprintf(stderr, "Iteration %d:\n", k);
    fprintf(stderr, "  fx = %f, x[0] = %f\n", fx, x[0]);
    fprintf(stderr, "  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    return 0;
}

// Main 
int main(int argc, char* argv[]) {
	/* L-BFGS parameters */
    int N=1;
	int ret = 0;
	lbfgsfloatval_t fx; /* Error value at current iteration*/
	lbfgsfloatval_t *x = lbfgs_malloc(N); /* Parameter vector */
	lbfgs_parameter_t param;  /* Step length at current iteration */

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
	param.orthantwise_c=0; /* Use L1 norm regularization */
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
	param.max_linesearch=10;
	param.max_iterations=10;
	param.m=10;

    /* Initialize x */
			x[0] = 1;

	/*
	   Start the L-BFGS optimization; this will invoke the callback functions
	   evaluate() and progress() when necessary.
	   */
	ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);

	/* Report the result. */
    fprintf(stderr, "L-BFGS optimization terminated with status code = %d\n", ret);
	fprintf(stderr, "  fx = %f, x[0] = %f\n", fx, x[0]);

	lbfgs_free(x);
	exit(0);
}
