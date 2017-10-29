// Header inclusion 
#include <stdio.h>
#include <stdlib.h>
#include <lbfgs.h>
#include <iostream>
#include "utils.h"
#include "data.h"
#include "file.h"
#include "image.h"
#include <memory>
#include <fstream>
#include <math.h>
#include <config4cpp/Configuration.h>

/* Global variables */
std::shared_ptr<rockseis::Image2D<float>> grad;
std::shared_ptr<rockseis::Data2D<float>> res;
int L, M;

static lbfgsfloatval_t evaluate(
		void *instance,
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		)
{
    fprintf(stderr, "Starting new evaluation\n");
    float *resdata = res->getData();
    float *graddata = grad->getImagedata();
    float Jd = 0;
	lbfgsfloatval_t fx = 0.0;
    for(int i=0; i < L; i++){
        Jd=0;
        for(int j=0; j < M; j++){
            Jd += x[M*i +j]*resdata[j];
        }
        for(int j=0; j < M; j++){
            g[M*i +j] = -2.0*(graddata[i] - Jd)*resdata[j];
        }
	    fx += (graddata[i] - Jd)*(graddata[i] - Jd);
    }

    fprintf(stderr, "  fx = %f, g[0] = %f\n", fx, g[0]);
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

	// Parameters
    std::string gradfile;
    std::string resfile;

	// Parse parameters from file
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	const char *     scope = "";
	const char *     configFile = "jacob.cfg";

    bool status = 0;
    try {
        cfg->parse(configFile);
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        return 1;
    }

    try {
        gradfile = cfg->lookupString(scope, "grad");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }

    try {
        resfile = cfg->lookupString(scope, "res");
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        status = 1;
    }
	// Destroy cfg
	cfg->destroy();

	if(status == 1){
		std::cerr << "Program terminated due to input errors." << std::endl;
		return 1;
	}


	// Create the classes 
    grad = std::make_shared<rockseis::Image2D<float>>(gradfile);
	res = std::make_shared<rockseis::Data2D<float>>(resfile);
    if(grad->read() == FILE_ERR) std::cerr << "Error reading gradient!" << std::endl;
    if(res->read() == FILE_ERR) std::cerr << "Error reading residual!" << std::endl;
    
    int nt = res->getNt();
    int ntr = res->getNtrace();
    int nx = grad->getNx();
    int nz = grad->getNz();

	/* L-BFGS parameters */
    int N=nx*nz*nt*ntr;
    L = nx*nz; 
    M = nt*ntr;

	int ret = 0;
	lbfgsfloatval_t fx; /* Error value at current iteration*/
	lbfgsfloatval_t *x = lbfgs_malloc(N); /* Parameter vector */
	lbfgs_parameter_t param;  /* Step length at current iteration */

	/* Initialize the parameters for the L-BFGS optimization. */
	lbfgs_parameter_init(&param);
	param.orthantwise_c=0; /* Use L1 norm regularization */
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
	param.max_linesearch=10;
	param.max_iterations=5;
	param.m=10;

    /* Initialize x */
    for(int i=0; i< N; i++) x[i] = 0;

	/*
	   Start the L-BFGS optimization; this will invoke the callback functions
	   evaluate() and progress() when necessary.
	   */
	ret = lbfgs(N, x, &fx, evaluate, progress, NULL, &param);

    float *graddata = grad->getImagedata();
    for(int i=0; i < L; i++){
        graddata[i] = 0.0;
        for(int j=0; j < M; j++){
            graddata[i] += x[M*i +j]*x[M*i + j];
        }
    }

    grad->setImagefile("Hess.rss");
    grad->write();


	/* Report the result. */
    fprintf(stderr, "L-BFGS optimization terminated with status code = %d\n", ret);
	fprintf(stderr, "  fx = %f, x[0] = %f\n", fx, x[0]);

	lbfgs_free(x);
	exit(0);
}
