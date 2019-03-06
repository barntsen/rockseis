/* Nonlinear optimization library
 
   Methods that are implemented
   1. Steepest-descent
   2. Conjugate gradient
   3. L-BFGS
   4. Heat bath simulated annealing
   5. Metropolis simulated annealing

   Author: Espen B. Raknes, espen.raknes@ntnu.no, NTNU, 2016
*/

/*
  Copyright (C) 2016 Norwegian University of Science and Technology 
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef OPT_H
#define OPT_H

/* Including libraries */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <memory>
#include <iostream>


namespace rockseis {

/* Constants */
#define OPT_STATUS_NOT_STARTED -1
#define OPT_STATUS_MINIMIZED 1
#define OPT_STATUS_MAXITERATIONS 2
#define OPT_STATUS_MAXLINESEARCH 3
#define OPT_STATUS_GRADIENT_DECREASE 4
#define OPT_STATUS_MAX_FUNCTIONAL_EVALUATIONS 5
#define OPT_STATUS_MINIMUM_TEMPERATURE 6
#define OPT_STATUS_PROBABILITY_PEAK 7
#define OPT_CONDITION_DECREASE 1
#define OPT_CONDITION_ARMIJO 2
#define OPT_CONDITION_WOLFE 3
#define OPT_CONDITION_STRONG_WOLFE 4
/*^*/

/* Defining structures */
typedef struct {
	int max_iterations;		// Maximum number of iterations
	int n;				// Length of vector
	int max_linesearch;		// Maximum number of linesearch trials accepted
	double xeps;			// Minimum decrease in model norm in linesearch update
	double geps;			// Minimum decrease in gradient norm
	int linesearch_condition;	// Which linesearch condition to use
	double fmin;			// The minimum value for f, default 1e-6
	double ftol;			// Factor to control decrease in f, default: 1e-4, see Nocedal and Wright, 2006, p33-34
	double gtol;			// Factor to control decrease in g, default: 0.9 for quasi-Newton methods and 0.1 for CG/steepest descent, see Nocedal and Wright, 2006, p33-34
	int ngradients;			// Number of gradients to keep for L-BFGS, default: 6
	double Tinit;			// Initial temperatur
	double Texp;			// Cooling scheme: T = Tinit*Texp^k (k is the iteration number)
	double Tmin;			// Minimum temperature where algorithm stops
	int maxFuncEval;		// Maximum number of functional values allowed (10000*n)
	int nTrials;			// The number of trials for each parameter
	int *xMask;			// Put to 0 if x[j] should not be optimized
	double *xMin;			// Minimum allowed value for x[j]
	double *xMax;			// Minimum allowed value for x[j]
	int nFuncEval;			// Number of functional evaluations performed
} OptParam;
/*^*/

typedef OptParam *OptParamPtr;
/*^*/

typedef struct {
	double *x;
	double f;
	int n;
} OptResult;
/*^*/

typedef OptResult *OptResultPtr;
/*^*/

typedef struct {
	double *x;
	double *g;
	double *tmp;
	double *pk;
	double beta;
	double f;
	double steplength;
	double T;
	int n;
} OptInstance;
/*^*/

typedef OptInstance *OptInstancePtr;
/*^*/

typedef struct {
	double *s;	// Array that holds all s vectors in ring storage
	double *y;	// Array that holds all y vectors in ring storage
	double *q;	// Temp array for computation of pk
	double *r;	// Temp array for computation of pk
	double *alpha;
	double *rho;
	int m;		// Number of previous iterations to store
	int n;		// Length of one vector
} OptLbfgsCont;
/*^*/

typedef OptLbfgsCont *OptLbfgsContPtr;
/*^*/

class Opt{
    public: 
        Opt(); ///< Default constructor
        Opt(const int n); ///<Constructor
        ~Opt(); ///< Destructor
        void opt_vector_copy(double *from, double *to, const int n);
        double opt_vector_norm(double *v, const int type, const int n);
        double opt_vector_diff_norm(double *v, double *u, const int type, const int n); 
        double opt_vector_dot(double *v, double *u, const int n); 
        OptResultPtr opt_result_init(const int n);
        void opt_result_free(OptResultPtr res);
        void opt_result_print(OptResultPtr res);
        OptParamPtr opt_param_init(const int n);
        void opt_param_free(OptParamPtr param);
        void opt_param_default(OptParamPtr param);
        OptInstancePtr opt_instance_init(const int n);
        void opt_instance_free(OptInstancePtr instance);
        void opt_instance_print(OptInstancePtr instance);
        void opt_instance_calculate_initial_cg_step(OptInstancePtr instance);
        void opt_instance_copy_g_to_pk(OptInstancePtr instance);
        void opt_instance_update_vector_g(OptInstancePtr current, OptInstancePtr next);
        void opt_instance_update_vector_pk(OptInstancePtr current, OptInstancePtr next);
        void opt_instance_update_cg_fr(OptInstancePtr current, OptInstancePtr next);
        void opt_instance_update_cg_pr(OptInstancePtr current, OptInstancePtr next);
        void opt_instance_set_steplength(OptInstancePtr previous, OptInstancePtr current);
        void opt_instance_set_steplength_cg(OptInstancePtr previous, OptInstancePtr current);
        void opt_instance_flip(OptInstancePtr current, OptInstancePtr next);
        void opt_copy_instance_to_result(OptResultPtr result, OptInstancePtr instance);
        OptLbfgsContPtr opt_lbfgs_container_init();
        void opt_lbfgs_container_free(OptLbfgsContPtr cont);
        double opt_lbfgs_calculate_gamma(OptLbfgsContPtr cont, const int iteration);
        void opt_lbfgs_calculate_initial_hessian(double *hessian, const int n, const double gamma); 
        void opt_lbfgs_calculate_pk(OptLbfgsContPtr cont, OptInstancePtr instance, double *hessian);
        void opt_lbfgs_container_update(OptLbfgsContPtr cont, OptInstancePtr current, OptInstancePtr next,  const int iteration);
        void opt_lbfgs_calculate_diag_hess(OptLbfgsContPtr cont, OptInstancePtr instance, double *hessian);
        void opt_set_initial_guess(double *x); 
        void opt_set_status_msg();
        int opt_linesearch_check(OptInstancePtr current, OptInstancePtr next);
        int opt_linesearch(OptInstancePtr current, OptInstancePtr next, void (*evaluate)(OptInstancePtr));
        int opt_check_optimization(OptInstancePtr current, OptInstancePtr next); 
        int opt_check_optimization_sa(OptInstancePtr next); 
        void opt_print_f();
        void opt_steepest_descent(void (*evaluate)(OptInstancePtr), void (*progress)(Opt*, OptInstancePtr));
        void opt_conjugate_gradient_fr(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr));
        void opt_conjugate_gradient_pr(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr));
        void opt_lbfgs(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr));
        void opt_heat_bath(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr));
        void opt_heat_bath_neighbour(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr), void (*neighbour)(OptInstancePtr, OptInstancePtr, OptParamPtr));
        void opt_metropolis(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr), void (*neighbour)(OptInstancePtr, OptInstancePtr, OptParamPtr));
        void opt_metropolis_random(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr));
        int getIter() { return iter; }
        char* getMsg() { return &msg[0]; }
        int  getStatus() { return status; }
        void setGtol(double gtol) { param->gtol = gtol;}
        void setLinesearch_condition(int condition) { param->linesearch_condition = condition;}
        void setMax_iterations(int max) {param->max_iterations = max;}
        void setMax_linesearch(int max) {param->max_linesearch = max;}

        int getMax_iterations() { return param->max_iterations; }
        int getMax_linesearch() { return param->max_linesearch; }

    private:
        OptParamPtr param;
        OptResultPtr result;
        double *init;
        int iter;
        int status;	// Status code
        char msg[128];	// Status message 
        double *f;	// Array with all
};

}
#endif // OPT_H
