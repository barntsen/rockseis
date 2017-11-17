// Include statements
#include "opt.h"

namespace rockseis {

//###### MISC functions
void Opt::opt_vector_copy(double *from, double *to, const int n)
/*<Copying two vectors>*/
{
	// Variables
	int i;

	for(i=0; i<n; i++) {
		to[i] = from[i];
	}
}

double Opt::opt_vector_norm(double *v, const int type, const int n) 
/*<Calculating norm of vector>*/
{
	// Variables
	int i;
	double norm;

	norm = 0.0;
	if(type == 2) {
		for(i=0; i<n; i++) {
			norm += v[i]*v[i];
		}
		norm = sqrtf(norm);
	}
	else if(type == 99999) {
		for(i=0; i<n; i++) {
			if(fabs(v[i]) >= norm) {
				norm = fabs(v[i]);
			}
		}
	}

	return norm;
}

double Opt::opt_vector_diff_norm(double *v, double *u, const int type, const int n) 
/*<Calculating norm of vector>*/
{
	// Variables
	int i;
	double norm;

	norm = 0.0;
	for(i=0; i<n; i++) {
		norm += (v[i]-u[i])*(v[i]-u[i]);
	}
	norm = sqrtf(norm);

	return norm;
}

double Opt::opt_vector_dot(double *v, double *u, const int n) 
/*<Calculating norm of vector>*/
{
	// Variables
	int i;
	double dot;

	dot = 0.0;
	for(i=0; i<n; i++) {
		dot += v[i]*u[i];
	}

	return dot;
}




//###### Results struct
OptResultPtr Opt::opt_result_init(const int n)
/*<Initializing the struct>*/
{
	// Variables
	OptResultPtr res;
	int i;

	// Malloc
	// Struct
	res = (OptResultPtr) malloc(sizeof(OptResult));
	if(res == NULL) {
		fprintf(stderr,"Error in memory allocation!\n");
		exit(0);
	}
	// Vectors
	res->x = (double *) malloc(sizeof(double)*n);
	
	for(i=0; i<n; i++) {
		res->x[i] = 0.0;
	}
	// Default values
	res->n = n;

	// Return
	return res;
}

void Opt::opt_result_free(OptResultPtr res)
/*<Cleaning memory>*/
{
	free(res->x);
	free(res);
}

void Opt::opt_result_print(OptResultPtr res)
/*<Printing to stderr of result struct>*/
{
	// Variables
	int i;

	fprintf(stderr,"******* RESULT INFO *******\n");
	fprintf(stderr,"n: %d\n",res->n);
	fprintf(stderr,"f: %10.8f\n",res->f);
	for(i=0; i<res->n; i++ ){
		if(((i) % 5) == 0 && i>1) fprintf(stderr,"\n");
		fprintf(stderr,"x[%03d]: %+6.6f   ",i,res->x[i]);
	}
	fprintf(stderr,"\n\n");
}



//###### Param struct
OptParamPtr Opt::opt_param_init(const int n)
/*<Initializing the struct>*/
{
	OptParamPtr param;

	// Malloc
	// Struct
	param = (OptParamPtr) malloc(sizeof(OptParam));
	if(param == NULL) {
		fprintf(stderr,"Error in memory allocation!\n");
		exit(0);
	}

	// Vectors
	param->xMask = (int *) malloc(sizeof(int)*n);
	param->xMin = (double *) malloc(sizeof(double)*n);
	param->xMax = (double *) malloc(sizeof(double)*n);

	// Default values
	param->n = n;

	// Return
	return param;
}

void Opt::opt_param_free(OptParamPtr param)
/*<Cleaning memory>*/
{
	free(param->xMask);
	free(param->xMin);
	free(param->xMax);
	free(param);
}

void Opt::opt_param_default(OptParamPtr param)
/*<Setting default values for the parameter struct>*/
{
	// Variables
	int i;

	param->max_iterations = 1000;
	param->max_linesearch = 20;
	param->xeps = 1e-10;
	param->geps = 1e-10;
	param->linesearch_condition = OPT_CONDITION_WOLFE;
	param->fmin = 1e-6;
	param->ftol = 1e-4;
	param->gtol = 0.9;
	param->ngradients = 6;
	
	param->Tinit = 1.0;
	param->Texp = 0.95;
	param->Tmin = 0.01;
	param->maxFuncEval = (param->n)*10000;
	param->nTrials = 10;
	param->nFuncEval = 0;
	
	for(i=0; i<param->n; i++) {
		param->xMask[i] = 1;
		param->xMin[i] = -10000.0;
		param->xMax[i] = 10000.0;
	}
}




//####### Instance struct
OptInstancePtr Opt::opt_instance_init(const int n)
/*<Initializing the struct>*/
{
	// Variables
	OptInstancePtr instance;
	int i;

	// Malloc
	// Struct
	instance = (OptInstancePtr) malloc(sizeof(OptInstance));
	if(instance == NULL) {
		fprintf(stderr,"Error in memory allocation!\n");
		exit(0);
	}
	// Vectors
	instance->x = (double *) malloc(sizeof(double)*n);
	instance->g = (double *) malloc(sizeof(double)*n);
	instance->pk = (double *) malloc(sizeof(double)*n);
	instance->tmp = (double *) malloc(sizeof(double)*n);
	
	for(i=0; i<n; i++) {
		instance->x[i] = 0.0;
		instance->g[i] = 0.0;
		instance->pk[i] = 0.0;
		instance->tmp[i] = 0.0;
	}
	// Default values
	instance->n = n;
	instance->f = 0.0;
	instance->steplength = 0.0;
	instance->T = 0.0;
	instance->beta = 0.0;

	// Return
	return instance;
}

void Opt::opt_instance_free(OptInstancePtr instance)
/*<Cleaning memory>*/
{
	free(instance->x);
	free(instance->g);
	free(instance->pk);
	free(instance->tmp);
	free(instance);
}

void Opt::opt_instance_print(OptInstancePtr instance)
/*<Print to stderr of instance>*/
{
	// Variables
	int i;

	fprintf(stderr,"******* OPT INSTANCE INFO *******\n");
	fprintf(stderr,"n: %d\n",instance->n);
	fprintf(stderr,"f: %10.8f\n",instance->f);
	fprintf(stderr,"steplength: %10.8f\n",instance->steplength);
	for(i=0; i<instance->n; i++ ){
		fprintf(stderr,"x[%d]: %6.6f g[%d]: %6.6f \t",i,instance->x[i],i,instance->g[i]);
		if((i%5) == 0 && i>0) fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n\n");
}

void Opt::opt_instance_calculate_initial_cg_step(OptInstancePtr instance)
/*<Calculates the initial search direction for CG.>*/
{
	int i;

	for(i=0; i<instance->n; i++) {
		instance->pk[i] = -(instance->g[i]);
	}
}

void Opt::opt_instance_copy_g_to_pk(OptInstancePtr instance)
/*<Copies the gradient to the search direction pk and scales with -1.0.>*/
{
	int i;

	for(i=0; i<instance->n; i++) {
		instance->pk[i] = -(instance->g[i]);
	}
}

void Opt::opt_instance_update_vector_g(OptInstancePtr current, OptInstancePtr next)
/*<Updating the vector using the gradient>*/
{
	// Variables
	int i;

	for(i=0; i<current->n; i++) {
		next->x[i] = current->x[i] - (current->steplength)*(current->g[i]);
	}
}

void Opt::opt_instance_update_vector_pk(OptInstancePtr current, OptInstancePtr next)
/*<Updating the vector using the search direction>*/
{
	// Variables
	int i;
	
	for(i=0; i<current->n; i++) {
		next->x[i] = current->x[i] + (current->steplength)*(current->pk[i]);
	}
	next->steplength = current->steplength;
}

void Opt::opt_instance_update_cg_fr(OptInstancePtr current, OptInstancePtr next)
/*<Updating the search directions using the Fletcher-Reeves formula>*/
{
	// Variables
	int i;
	double num = 0.0;
	double den = 0.0;
	
	// Find beta
	next->beta = 0.0;
	for(i=0; i<current->n; i++) {
		num += (next->g[i])*(next->g[i]);
		den += (current->g[i])*(current->g[i]);
	}
	next->beta = num/den;

	// Update search direction
	for(i=0; i<current->n; i++) {
		next->pk[i] = -(next->g[i]) + (next->beta)*(current->pk[i]);
	}
}

void Opt::opt_instance_update_cg_pr(OptInstancePtr current, OptInstancePtr next)
/*<Updating the search directions using the Polak-Ribiere formula>*/
{
	// Variables
	int i;
	double num = 0.0;
	double den = 0.0;
	
	// Find beta
	next->beta = 0.0;
	for(i=0; i<current->n; i++) {
		num += (next->g[i])*(next->g[i] - current->g[i]);
		den += (current->g[i])*(current->g[i]);
	}
	next->beta = num/den;

	// Update search direction
	for(i=0; i<current->n; i++) {
		next->pk[i] = -(next->g[i]) + (next->beta)*(current->pk[i]);
	}
}

void Opt::opt_instance_set_steplength(OptInstancePtr previous, OptInstancePtr current)
/*<Sets the next initial step length guess.>*/
{
	// Variables
	int i;
	double num,den;

	// Calculate numerator and denumerator
	num = 0.0;
	den = 0.0;
	for(i=0; i<current->n; i++) {
		num += (previous->g[i])*(previous->g[i]);
		den += (current->g[i])*(current->g[i]);
	}

	// Calculate next step
	current->steplength = previous->steplength*(num/den);
}

void Opt::opt_instance_set_steplength_cg(OptInstancePtr previous, OptInstancePtr current)
/*<Sets the next initial step length guess.>*/
{
	// Variables
	int i;
	double num,den;

	// Calculate numerator and denumerator
	num = 0.0;
	den = 0.0;
	for(i=0; i<current->n; i++) {
		num += (previous->g[i])*(previous->pk[i]);
		den += (current->g[i])*(current->pk[i]);
	}

	// Calculate next step
	current->steplength = previous->steplength*(num/den);
}

void Opt::opt_instance_flip(OptInstancePtr current, OptInstancePtr next)
/*<Flip values in instances.>*/
{
	// Variables
	double temp;

	opt_vector_copy(next->x,next->tmp,next->n);
	opt_vector_copy(current->x,next->x,next->n);
	opt_vector_copy(next->tmp,current->x,next->n);

	opt_vector_copy(next->pk,next->tmp,next->n);
	opt_vector_copy(current->pk,next->pk,next->n);
	opt_vector_copy(next->tmp,current->pk,next->n);
	
	opt_vector_copy(next->g,next->tmp,next->n);
	opt_vector_copy(current->g,next->g,next->n);
	opt_vector_copy(next->tmp,current->g,next->n);
	
	temp = current->f;
	current->f = next->f;
	next->f = temp;
	
	next->steplength = current->steplength;
}

void Opt::opt_copy_instance_to_result(OptResultPtr result, OptInstancePtr instance)
/*<Copying the instance to result>*/
{
	// Copying vector
	opt_vector_copy(instance->x,result->x,instance->n);

	// Misfit
	result->f = instance->f;
}



//####### Opt L-BFGS container
OptLbfgsContPtr Opt::opt_lbfgs_container_init()
/*<Initializing of the L-BFGS memory container>*/
{
	// Variables
	OptLbfgsContPtr cont;
	int i;
	
	// Malloc
	cont = (OptLbfgsContPtr) malloc(sizeof(OptLbfgsCont));
	if(cont == NULL) {
		fprintf(stderr,"Error in memory allocation!\n");
		exit(0);
	}

	// Malloc
	cont->s = (double *) malloc(sizeof(double)*(this->param->n)*(this->param->ngradients));
	cont->y = (double *) malloc(sizeof(double)*(this->param->n)*(this->param->ngradients));
	cont->q = (double *) malloc(sizeof(double)*(this->param->n));
	cont->r = (double *) malloc(sizeof(double)*(this->param->n));
	cont->alpha = (double *) malloc(sizeof(double)*(this->param->ngradients));
	cont->rho = (double *) malloc(sizeof(double)*(this->param->ngradients));
	if(cont->s == NULL || cont->y == NULL || cont->q == NULL || cont->r == NULL || cont->alpha == NULL || cont->rho == NULL) {
		fprintf(stderr,"Error in memory allocation!\n");
		exit(0);
	}

	for(i=0; i<(this->param->n)*(this->param->ngradients); i++) {
		cont->s[i] = 0.0;
		cont->y[i] = 0.0;
		if(i<this->param->n){
			cont->q[i] = 0.0;
			cont->r[i] = 0.0;
		}

		if(i<this->param->ngradients){
			cont->alpha[i] = 0.0;
			cont->rho[i] = 0.0;
		}
	}

	// Default values
	cont->m = this->param->ngradients;
	cont->n = this->param->n;

	// Return
	return cont;
}

void Opt::opt_lbfgs_container_free(OptLbfgsContPtr cont)
/*<Clearing memory>*/
{
	free(cont->s);
	free(cont->y);
	free(cont->q);
	free(cont->r);
	free(cont->alpha);
	free(cont->rho);
	free(cont);
}

double Opt::opt_lbfgs_calculate_gamma(OptLbfgsContPtr cont, const int iteration)
/*<Calculates the gamma variable for initial H.>*/
{
	// Variables
	double gamma,den,num;
	int prev;

	if(iteration == 0) {
		gamma = 1.0;
	}
	else {
		num = 0.0;
		den = 0.0;
		prev = (iteration%(cont->m))-1;
		if(prev < 0) prev += cont->m;

		
		num = opt_vector_dot(cont->s + prev*(cont->n),cont->y + prev*(cont->n),cont->n);
		den = opt_vector_dot(cont->y + prev*(cont->n),cont->y + prev*(cont->n),cont->n);
		gamma = num/den;
	}

	return gamma;
}

void Opt::opt_lbfgs_calculate_initial_hessian(double *hessian, const int n, const double gamma) 
/*<Calculates the initial Hessian in an iteration>*/
{
	// Variables
	int i;

	for(i=0; i<n; i++) {
		hessian[i] = gamma;
	}
}

void Opt::opt_lbfgs_calculate_pk(OptLbfgsContPtr cont, OptInstancePtr instance, double *hessian)
/*<Calculates the search direction pk for L-BFGS.>*/
{
	// Variables
	int i,j,mod;
	double dot, beta;


	// Save initial r (the gradient)
	opt_vector_copy(instance->g,cont->q,cont->n);
	
	// Calcuates q
	for(i=0; i<cont->m; i++) {
		// Calculates alpha_i
		mod = (cont->m-1-i)%cont->m;
		dot = opt_vector_dot(cont->s + mod*cont->n,cont->q,cont->n);
		cont->alpha[cont->m-1-i] = cont->rho[cont->m-1-i]*dot;
		
		// Updates q
		for(j=0; j<cont->n; j++) {
			cont->q[j] = cont->q[j] - cont->alpha[cont->m-1-i]*cont->y[mod*cont->n + j];
		}
	}
	
	// Saves initial r
	for(i=0; i<cont->n; i++) {
		cont->r[i] = hessian[i]*cont->q[i];
	}

	// Calculates r
	for(i=0; i<cont->m; i++) {
		// Calculates beta
		mod = i%cont->m;
		dot = opt_vector_dot(cont->y + mod*cont->n,cont->r,cont->n);
		beta = cont->rho[i]*dot;

		// Updates r
		for(j=0; j<cont->n; j++) {
			cont->r[j] = cont->r[j] + cont->s[mod*cont->n + j]*(cont->alpha[i]-beta);
		}
	}
	
	// Store to pk
	for(i=0; i<cont->n; i++) {
		instance->pk[i] = -1.0*cont->r[i];
	}
}

void Opt::opt_lbfgs_container_update(OptLbfgsContPtr cont, OptInstancePtr current, OptInstancePtr next,  const int iteration)
/*<Update the container with new information>*/
{
	// Variables
	int i,mod;
	double dot;

	// Get the position in the ring storage
	mod = iteration%cont->m;
	
	// Updates s and y
	for(i=0; i<cont->n; i++) {
		cont->y[mod*cont->n + i] = next->g[i] - current->g[i];
		cont->s[mod*cont->n + i] = next->x[i] - current->x[i];
	}
	
	// Updates rho
	dot = opt_vector_dot(cont->y + mod*cont->n,cont->s + mod*cont->n,cont->n);
	cont->rho[mod] = 1.0/dot;
}







//####### Opt Constructor
//
Opt::Opt() 
{
	// Variables
	int i;
    int n=1;

	this->param = opt_param_init(n);
	this->result = opt_result_init(n);
	this->init = (double *) malloc(sizeof(double)*n);
	for(i=0; i<n; i++) {
		this->init[i] = 0.0;
	}

	// Setting values
	opt_param_default(this->param); // Default values
	this->iter = 1;
	this->status = OPT_STATUS_NOT_STARTED;

	// Creating array for f-s
	this->f = (double *) malloc(sizeof(double)*(this->param->max_iterations));
	for(i=0; i<(this->param->max_iterations); i++) {
		this->f[i] = 0.0;
	}
}

//####### Opt Constructor
Opt::Opt(const int n) 
{
	// Variables
	int i;

	this->param = opt_param_init(n);
	this->result = opt_result_init(n);
	this->init = (double *) malloc(sizeof(double)*n);
	for(i=0; i<n; i++) {
		this->init[i] = 0.0;
	}

	// Setting values
	opt_param_default(this->param); // Default values
	this->iter = 1;
	this->status = OPT_STATUS_NOT_STARTED;

	// Creating array for f-s
	this->f = (double *) malloc(sizeof(double)*(this->param->max_iterations));
	for(i=0; i<(this->param->max_iterations); i++) {
		this->f[i] = 0.0;
	}
}

//####### Opt Destructor
Opt::~Opt()
{
    free(f);
	opt_param_free(param);
	opt_result_free(result);
}

void Opt::opt_set_initial_guess(double *x) 
/*<Setting initial guess>*/
{
	// Variables
	int i;

	for(i=0; i<this->param->n; i++) {
		this->init[i] = x[i];
	}
}

void Opt::opt_set_status_msg()
/*<Setting the status message>*/
{
	switch(this->status) {
		case OPT_STATUS_MAXLINESEARCH:
			snprintf(this->msg,128,"Maximum number of linesearches is performed.");
			break;
		case OPT_STATUS_MAXITERATIONS:
			snprintf(this->msg,128,"Maximum number of iterations is performed.");
			break;
		case OPT_STATUS_MINIMIZED:
			snprintf(this->msg,128,"Functional is minimized.");
			break;
		case OPT_STATUS_GRADIENT_DECREASE:
			snprintf(this->msg,128,"Too small change in gradient norm.");
			break;
		case OPT_STATUS_MAX_FUNCTIONAL_EVALUATIONS:
			snprintf(this->msg,128,"Maximum number of functional evalations is performed.");
			break;
		case OPT_STATUS_MINIMUM_TEMPERATURE:
			snprintf(this->msg,128,"Minimum temperature is achieved.");
			break;
		case OPT_STATUS_PROBABILITY_PEAK:
			snprintf(this->msg,128,"Probability peak observed. Functional is minimized.");
			break;
		case OPT_STATUS_NOT_STARTED:
		default:
			snprintf(this->msg,128,"Optimization method not started.");
			break;
	}
}

int Opt::opt_linesearch_check(OptInstancePtr current, OptInstancePtr next)
/*<Checks if we have a sufficient decrease in the functional.>*/
{
	// Variables
	int ret;
	double dot1,dot2;

	// Default
	ret = 0;

	if(this->param->linesearch_condition == OPT_CONDITION_DECREASE) {
		// Decrease condition
		if(next->f < current->f) {
			ret = 1;
		}		
	}
	else if(this->param->linesearch_condition == OPT_CONDITION_ARMIJO) {
		// Armijo condition
		dot1 = opt_vector_dot(current->g,current->pk,current->n);

		if(next->f <= (current->f + (this->param->ftol)*(current->steplength)*dot1)) {
			ret = 1;
		}		
	}
	else if(this->param->linesearch_condition == OPT_CONDITION_WOLFE) {
		// Wolfe condition
		dot1 = opt_vector_dot(current->g,current->pk,current->n);
		dot2 = opt_vector_dot(next->g,current->pk,current->n);
		
		if((next->f <= (current->f + (this->param->ftol)*(current->steplength)*dot1)) && (dot2 >= ((this->param->gtol)*dot1))) {
			ret = 1;
		}
	}
	else {
		// Wolfe condition
		dot1 = opt_vector_dot(current->g,current->pk,current->n);
		dot2 = opt_vector_dot(next->g,current->pk,current->n);
		
		if((next->f <= (current->f + (this->param->ftol)*(current->steplength)*dot1)) && (fabs(dot2) <= ((this->param->gtol)*fabs(dot1)))) {
			ret = 1;
		}
	}

	// Returning
	return ret;
}

int Opt::opt_linesearch(OptInstancePtr current, OptInstancePtr next, void (*evaluate)(OptInstancePtr))
/*<Performs the linesearch and finds a step length.>*/
{
	// Variables
	int ret,maxlinesearch;
	int i;

	// Setting default values
	ret = 0;
	maxlinesearch = this->param->max_linesearch;

	// Performing linesearch
	for(i=0; i<maxlinesearch; i++) {
		// Update to next vector
		opt_instance_update_vector_pk(current,next);

		// Evaluates the next vector
		evaluate(next);
		this->param->nFuncEval++; // Updating counter for number of functional evaluations
		
		// Check linesearch condition
		if(opt_linesearch_check(current,next) != 0) {
			ret = 1;
			break;
		}
		else {
			current->steplength = 0.5*current->steplength;
		}
	}

	// Return
	return ret;
}

int Opt::opt_check_optimization(OptInstancePtr current, OptInstancePtr next) 
/*<Checks if we should stop the optimization method.>*/
{
	// Variables
	int val;
	double norm;

	// Defaults
	val = 0;
	norm = 0.0;
	
	// Checking if vector updates are large enough
	norm = opt_vector_diff_norm(current->x,next->x,2,current->n);
	if(norm < this->param->xeps) {
		this->status = OPT_STATUS_MINIMIZED;
		val = 1;
	}

	// Checking if gradient norm is small (no change in updates)
	norm = opt_vector_norm(next->g,2,next->n);
	if(norm < this->param->geps) {
		this->status = OPT_STATUS_GRADIENT_DECREASE;
		val = 1;
	}

	// Checking if we have performed max number of iterations
	if(this->iter >= this->param->max_iterations) {
		this->status = OPT_STATUS_MAXITERATIONS;
		val = 1;
	}
	
	// Returning
	return val;
}

int Opt::opt_check_optimization_sa(OptInstancePtr next) 
/*<Checks if we should stop the simulated annealing optimization method.>*/
{
	fprintf(stderr,"Checking stopping criteria\n");
	// Variables
	int val = 0;
	
	// Checking if we have performed maximum number of function evaluations
	if(next->f < this->param->fmin) {
		this->status = OPT_STATUS_MINIMIZED;
		val = 1;
	}
	
	// Checking if we have performed maximum number of function evaluations
	if(this->param->nFuncEval > this->param->maxFuncEval) {
		this->status = OPT_STATUS_MAX_FUNCTIONAL_EVALUATIONS;
		val = 1;
	}
	
	// Checking if we have reached the minimum temperature
	if(next->T < this->param->Tmin) {
		this->status = OPT_STATUS_MINIMUM_TEMPERATURE;
		val = 1;
	}

	// Checking if we have performed max number of iterations
	if(this->iter >= this->param->max_iterations) {
		this->status = OPT_STATUS_MAXITERATIONS;
		val = 1;
	}
	
	// Returning
	return val;
}

void Opt::opt_print_f()
/*<Prints out all the values of f>*/
{
	// Variables
	int i;

	fprintf(stderr,"##### Values of the functional: #####\n");
	for(i=0; i<this->iter; i++) {
		if((i%6 == 0) && i>0) fprintf(stderr,"\n");

		fprintf(stderr,"f[%03d]: %+15.10e  ",i,this->f[i]);
	}
	fprintf(stderr,"\n\n");
}



//####### Optimization methods
void Opt::opt_steepest_descent(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running steepest descent optimization>*/
{
	// Variables
	int i,ls,status;
	OptInstancePtr current,next;
	double steplength;

	// Initializing instances
	current = opt_instance_init(this->param->n);
	next = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Initial evaluation
	evaluate(current);
	this->param->nFuncEval++;
	opt_instance_copy_g_to_pk(current);
	this->f[0] = current->f;

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		// Find initial step length for line search
		if(i == 0) {
			steplength = opt_vector_norm(current->g,2,current->n);
			current->steplength = 1.0/steplength;
		}
		else {	
			opt_instance_set_steplength(next,current);
		}

		// Perform linesearch;
		ls = opt_linesearch(current,next,evaluate);
		if(ls == 0) {
			// Stopping if maximum linesearch trials have been performed
			opt_copy_instance_to_result(this->result,current);
			this->status = OPT_STATUS_MAXLINESEARCH;
			break;
		}
		
		// Taking the step
		progress(this, next);
		this->f[this->iter] = next->f;
		this->iter++;
		
		// Checking if we should stop
		status = opt_check_optimization(current,next);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,next);
			break;
		}

		// Update search direction
		opt_instance_copy_g_to_pk(next);

		// Copying information from next to current to avoid tmp structs
		opt_instance_flip(current,next);
	}

	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(current);
	opt_instance_free(next);
}

void Opt::opt_conjugate_gradient_fr(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running non-linear CG optimization using the Fletcher-Reeves method>*/
{
	// Variables
	int i,ls,status;
	OptInstancePtr current,next;
	double steplength;

	// Initializing instances
	current = opt_instance_init(this->param->n);
	next = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Initial evaluation
	evaluate(current);
	this->param->nFuncEval++;
	opt_instance_calculate_initial_cg_step(current);
	this->f[0] = current->f;

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		// Find initial step length for line search
		if(i == 0) {
			steplength = opt_vector_norm(current->g,2,current->n);
			current->steplength = 1.0/(steplength);
		}
		else {	
			opt_instance_set_steplength_cg(next,current);
		}

		// Perform linesearch
		ls = opt_linesearch(current,next,evaluate);
		if(ls == 0) {
			// Stopping if maximum linesearch trials have been performed
			opt_copy_instance_to_result(this->result,current);
			this->status = OPT_STATUS_MAXLINESEARCH;
			break;
		}

		// Taking the step
		progress(this, next);
		this->f[this->iter] = next->f;
		this->iter++;
		
		// Checking if we should stop
		status = opt_check_optimization(current,next);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,next);
			break;
		}
		
		// Update beta, pk
		opt_instance_update_cg_fr(current,next);

		// Copying information from next to current to avoid tmp structs
		opt_instance_flip(current,next);
	}

	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(current);
	opt_instance_free(next);
}

void Opt::opt_conjugate_gradient_pr(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running non-linear CG optimization using the Polak-Ribiere method>*/
{
	// Variables
	int i,ls,status;
	OptInstancePtr current,next;
	double steplength;

	// Initializing instances
	current = opt_instance_init(this->param->n);
	next = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Initial evaluation
	evaluate(current);
	this->param->nFuncEval++;
	opt_instance_calculate_initial_cg_step(current);
	this->f[0] = current->f;

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		// Find initial step length for line search
		if(i == 0) {
			steplength = opt_vector_norm(current->g,2,current->n);
			current->steplength = 1.0/(steplength);
		}
		else {	
			opt_instance_set_steplength_cg(next,current);
		}

		// Perform linesearch
		ls = opt_linesearch(current,next,evaluate);
		if(ls == 0) {
			// Stopping if maximum linesearch trials have been performed
			opt_copy_instance_to_result(this->result,current);
			this->status = OPT_STATUS_MAXLINESEARCH;
			break;
		}

		// Taking the step
		progress(this, next);
		this->f[this->iter] = next->f;
		this->iter++;
		
		// Checking if we should stop
		status = opt_check_optimization(current,next);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,next);
			break;
		}
		
		// Update beta, pk
		opt_instance_update_cg_pr(current,next);

		// Copying information from next to current to avoid tmp structs
		opt_instance_flip(current,next);
	}

	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(current);
	opt_instance_free(next);
}

void Opt::opt_lbfgs(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running L-BFGS optimization>*/
{
	// Variables
	int i,ls,status;
	OptInstancePtr current,next;
	OptLbfgsContPtr cont;
	double steplength,gamma;
	double *hessian;

	// Initializing structs and vectors
	current = opt_instance_init(this->param->n);
	next = opt_instance_init(this->param->n);
	cont = opt_lbfgs_container_init();
	hessian = (double *) malloc(sizeof(double*)*(this->param->n));
	if(hessian == NULL) {
		fprintf(stderr,"Error in malloc!\n");
		exit(0);
	}

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Initial evaluation
	evaluate(current);
	this->param->nFuncEval++;
	this->f[0] = current->f;

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		// Find initial step length for line search
		if(i == 0) {
			steplength = opt_vector_norm(current->g,2,current->n);
			current->steplength = 1.0/(steplength);
		}
		else {	
			current->steplength = 1.0;
		}

		// Calculates the initial Hessian
		gamma = opt_lbfgs_calculate_gamma(cont,i);
		opt_lbfgs_calculate_initial_hessian(hessian,this->param->n,gamma);

		// Get search direction
		opt_lbfgs_calculate_pk(cont,current,hessian);

		// Perform linesearch
		ls = opt_linesearch(current,next,evaluate);
		if(ls == 0) {
			// Stopping if maximum linesearch trials have been performed
			opt_copy_instance_to_result(this->result,current);
			this->status = OPT_STATUS_MAXLINESEARCH;
			break;
		}
		
		// Taking the step
		progress(this, next);
		this->f[this->iter] = next->f;
		this->iter++;
		
		// Checking if we should stop
		status = opt_check_optimization(current,next);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,next);
			break;
		}
		
		// Update L-BFGS container
		opt_lbfgs_container_update(cont,current,next,i);
		
		// Copying information from next to current to avoid tmp structs
		opt_instance_flip(current,next);
	}
	
	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(current);
	opt_instance_free(next);
	opt_lbfgs_container_free(cont);
	free(hessian);
}

void Opt::opt_heat_bath(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running heat bath simulated annealing optimiziation>*/
{
	// Variables
	int i,j,k,l,status;
	OptInstancePtr next;

	// Initializing instances
	next = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,next->x,this->param->n);

	// Initial evaluation
	evaluate(next);
	this->param->nFuncEval++;
	this->f[0] = next->f;

	// Get variables from container
	int nTrials = this->param->nTrials;
	double T = this->param->Tinit;
	double dec  = this->param->Texp;
	
	// Malloc
	double *P = (double *) malloc(sizeof(double)*nTrials); // The probability function
	double *C = (double *) malloc(sizeof(double)*nTrials); // The cumulative probability function
	double *fTrial = (double *) malloc(sizeof(double)*nTrials); // The functional values
	double *xTrial = (double *) malloc(sizeof(double)*nTrials); // The model values at each trial

	for(k=0; k<nTrials; k++) {
		P[k] = 0.0;
		C[k] = 0.0;
		fTrial[k] = 0.0;
		xTrial[k] = 0.0;
	}

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		bool PPeak = 0; // If we have a peak in the probability and should stop

		for(j=0; j<(this->param->n); j++) {
			if(this->param->xMask[j] == 0) {
				continue;
			}
			else {
				double Psum = 0.0;
				double dx = (this->param->xMax[j] - this->param->xMin[j])/nTrials;
				
				fprintf(stderr,"\n\nLIBOPT:\n");
				for(k=0; k<nTrials; k++) {
					// Computing new value
					xTrial[k] = this->param->xMin[j] + dx*k;

					fprintf(stderr,"Index: %d, xVal: %f\n",k,xTrial[k]);

					// Updating trial
					next->x[j] = xTrial[k];

					// Evaluate
					evaluate(next);
					this->param->nFuncEval++;

					// Get functional value
					fTrial[k] = next->f;

					// Compute probability
					P[k] = expf(-(fTrial[k]/T));
					Psum += P[k];
				}

				// Updating probability functions
				for(k=0; k<nTrials; k++) {
					P[k] /= Psum;
				}

				for(k=0; k<nTrials; k++){
					for(l=0; l<k+1; l++) {
						C[k] += P[l];
					}
				}
				
				for(k=0; k<nTrials; k++) {
					fprintf(stderr,"C[%d]: %f,  ",k,C[k]);
				}
				fprintf(stderr,"\n");

				// Draw a random number
				srand(time(NULL) + j + i + i*j); // Use time right now as seed
				double r =((double) rand()/(RAND_MAX));
				
				// Finding which value we should use
				int ind = 0;
				for(k=0; k<nTrials; k++) {
					if(C[k] > r) {
						ind=k;
						if(ind < 0) ind = 0;
						break;
					}
				}
				
				fprintf(stderr,"r: %f, index: %d\n\n",r,ind);
				
				// Updating parameter
				next->x[j] = xTrial[ind];
				next->f = fTrial[ind];

				// Checking if we should stop
				// If enough function evaluations
				if(this->param->nFuncEval > this->param->maxFuncEval) {
					break;
				}
				// if some of the probabilites is 1.0 or nan
				for(k=0; k<nTrials; k++) {
					if(P[k] == 1.0 || isnan(P[k])) {
						PPeak = true;
						break;
					}
				}
				if(PPeak) break;

				// Zeroing vectors	
				for(k=0; k<nTrials; k++) {
					P[k] = 0.0;
					C[k] = 0.0;
					fTrial[k] = 0.0;
					xTrial[k] = 0.0;
				}
			}
		}
		
		// Taking the step
		next->T = T;
		progress(this, next);
		this->f[this->iter] = next->f;
		this->iter++;
		
		// Checking if we should stop
		status = opt_check_optimization_sa(next);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,next);
			break;
		}
		if(PPeak) {
			opt_copy_instance_to_result(this->result,next);
			this->status = OPT_STATUS_PROBABILITY_PEAK;
			break;
		}
		
		// Updating temperature
		T = T*dec;
	}

	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(next);
	free(P);
	free(C);
	free(fTrial);
	free(xTrial);
}

void Opt::opt_heat_bath_neighbour(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr), void (*neighbour)(OptInstancePtr, OptInstancePtr, OptParamPtr))
/*<Running heat bath simulated annealing optimiziation where the perturbation is performed by input function>*/
{
	// Variables
	int i,k,l,status;
	OptInstancePtr current;

	// Initializing instances
	current = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Initial evaluation
	evaluate(current);
	this->param->nFuncEval++;
	this->f[0] = current->f;

	// Get variables from container
	int nTrials = this->param->nTrials;
	double T = this->param->Tinit;
	double dec  = this->param->Texp;
	
	// Malloc
	double *P = (double *) malloc(sizeof(double)*nTrials); // The probability function
	double *C = (double *) malloc(sizeof(double)*nTrials); // The cumulative probability function
	OptInstancePtr *trial = (OptInstancePtr *) malloc(sizeof(OptInstancePtr)*nTrials); // The model values at each trial

	for(k=0; k<nTrials; k++) {
		P[k] = 0.0;
		C[k] = 0.0;
		trial[k] = opt_instance_init(this->param->n);
	}

	// Performing iterations
	for(i=0; i<this->param->max_iterations; i++) {
		bool PPeak = 0; // If we have a peak in the probability and should stop
		double Psum = 0.0;

		for(k=0; k<this->param->nTrials; k++) {
			// Finding next neighbour
			neighbour(trial[k],current,this->param);

			// Evaluate functional
			trial[k]->T = T;
			evaluate(trial[k]);
			this->param->nFuncEval++;

			// Compute probability
			P[k] = expf(-((trial[k]->f)/T));
			Psum += P[k];

		}
				
		// Updating probability functions
		fprintf(stderr,"\n\nLIBOPT:\n");
		for(k=0; k<nTrials; k++) {
			P[k] /= Psum;
			fprintf(stderr,"P[%d]: %f  f[%d]: %f, T: %f \n",k,P[k],k,trial[k]->f,T);
		}
		fprintf(stderr,"\n");
				
		for(k=0; k<nTrials; k++){
			for(l=0; l<k+1; l++) {
				C[k] += P[l];
			}
			fprintf(stderr,"C[%d]: %f  ",k,C[k]);
		}

		// Draw a random number
		srand(time(NULL) + i ); // Use time right now as seed
		double r =((double) rand()/(RAND_MAX));
		
		// Finding which value we should use
		int ind = 0;
		for(k=0; k<nTrials; k++) {
			if(C[k] > r) {
				ind=k;
				if(ind < 0) ind = 0;
				break;
			}
		}
		fprintf(stderr,"r: %f, ind: %d\n",r,ind);
		
		// Checking if we should stop
		// If enough function evaluations
		if(this->param->nFuncEval > this->param->maxFuncEval) {
			break;
		}
		// if some of the probabilites is 1.0 or nan
		for(k=0; k<nTrials; k++) {
			if(P[k] == 1.0 || isnan(P[k])) {
				PPeak = true;
				break;
			}
		}

		// Zeroing vectors	
		for(k=0; k<nTrials; k++) {
			P[k] = 0.0;
			C[k] = 0.0;
		}

		// Checking if we should stop
		status = opt_check_optimization_sa(trial[ind]);
		if(status != 0) {
			opt_copy_instance_to_result(this->result,trial[ind]);
			break;
		}
		if(PPeak) {
			opt_copy_instance_to_result(this->result,trial[ind]);
			this->status = OPT_STATUS_PROBABILITY_PEAK;
			break;
		}
		
		// Taking the step
		progress(this, trial[ind]);
		this->f[this->iter] = trial[ind]->f;
		this->iter++;

		// Copying the result to current value
		opt_vector_copy(trial[ind]->x,current->x,current->n);
		
		// Updating temperature
		T = T*dec;
	}

	// Setting status message
	opt_set_status_msg();

	// Clearing results
	for(k=0; k<nTrials; k++) opt_instance_free(trial[k]);
	free(P);
	free(C);
}

void Opt::opt_metropolis(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr), void (*neighbour)(OptInstancePtr, OptInstancePtr, OptParamPtr))
/*<Running Metropolis simulated annealing optimiziation>*/
{
	// Variables
	int i,l,status;
	OptInstancePtr next,current;

	// Initializing instances
	next = opt_instance_init(this->param->n);
	current = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,current->x,this->param->n);

	// Get variables from container
	double T = this->param->Tinit;
	double dec  = this->param->Texp;
	
	// Initial evaluation
	current->T = T;
	evaluate(current);
	this->param->nFuncEval++;
	this->f[0] = current->f;
	double fCurr = current->f;
	double fNext = 0.0;

	// Performing iterations
	l=0; 		// Counter for performed iterations
	for(i=0; i<this->param->max_iterations; i++) {
		for(int iTrials=0; iTrials<this->param->nTrials; iTrials++) {
			// Finding next neighbour
			neighbour(next,current,this->param);
			
			// Evaluate functional
			next->T = T;
			evaluate(next);
			this->param->nFuncEval++;

			// Find difference
			fNext = next->f;
			double fDiff = fNext - fCurr;

			// Compute probability
			double P = expf(-fDiff/T);

			fprintf(stderr,"LIBOPT: fcurr:%f, fNext: %f, fdiff:%f, T: %f, P:%f\n",fCurr,fNext,fDiff,T,P);

			// Make decision
			if(fDiff <= 0.0) {
				fprintf(stderr,"\n\nDECISION MAKING: Take new step\n\n");
				progress(this, next);
				this->f[this->iter] = next->f;
				this->iter++;
				fCurr = fNext;		
				l++;
				
				// Copy to the returning instance
				opt_vector_copy(next->x,current->x,next->n);
				current->f = fCurr;
				current->T = T;
				break;
			}
			else {
				// Draw a random number
				srand(time(NULL)); // Random number
				double r = ((double) rand()/(RAND_MAX));
				
				fprintf(stderr,"\n\nDECISION MAKING: Probability step\n");

				if(P > r) {
					fprintf(stderr,"DECISION MAKING: Yes\n\n");
					// We accept the guess
					progress(this, next);
					this->f[this->iter] = next->f;
					this->iter++;
					fCurr = fNext;		
					l++;

					// Copy to the returning instance
					opt_vector_copy(next->x,current->x,next->n);
					current->f = fCurr;
					current->T = T;
					break;
				}
				else {
					fprintf(stderr,"DECISION MAKING: No\n\n");
					// Do not accept guess
					opt_vector_copy(current->x,next->x,next->n);
				}
			}

			// Checking if we should stop
			status = opt_check_optimization_sa(next);
			if(status != 0) {
				opt_copy_instance_to_result(this->result,current);
				break;
			}
		}
		
		// Updating temperature
		fprintf(stderr,"\n\nLIBOPT: Lowering temperature from %f to %f.\n\n",T,T*dec);
		T = T*dec;
	}

	if(i >= this->param->max_iterations) {
		this->status = OPT_STATUS_MAXITERATIONS;
		opt_copy_instance_to_result(this->result,current);
	}


	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(next);
	opt_instance_free(current);
}

void Opt::opt_metropolis_random(void (*evaluate)(OptInstancePtr), void (*progress)(Opt *, OptInstancePtr))
/*<Running Metropolis simulated annealing optimiziation where the model parameters are perturbated randomly>*/
{
	// Variables
	int i,j,l,k,status;
	OptInstancePtr next,current;

	// Initializing instances
	next = opt_instance_init(this->param->n);
	current = opt_instance_init(this->param->n);

	// Copying initial model into current instance
	opt_vector_copy(this->init,next->x,this->param->n);

	// Initial evaluation
	evaluate(next);
	this->param->nFuncEval++;
	this->f[0] = next->f;

	// Get variables from container
	double T = this->param->Tinit;
	double dec  = this->param->Texp;
	double fCurr = next->f;
	double fNext = 0.0;
	double *curr = (double *) malloc(sizeof(double)*this->param->n);

	for(i=0; i<this->param->n; i++) {
		curr[i] = 0.0;
	}

	// Performing iterations
	l=0; 		// Counter for performed iterations
	int seed = 0; 	// For getting a "random" number if iterations go very quickly 
	for(i=0; i<this->param->max_iterations; i++) {
		for(int iTrials=0; iTrials<this->param->nTrials; iTrials++) {
			for(j=0; j<(this->param->n); j++) {
				if(this->param->xMask[j] == 0) {
					curr[j] = next->x[j]; // Current value
				}
				else {
					// Creating a new guess for parameter j
					srand(time(NULL) + j + i + i*j + seed); // Random number
					seed++;
					double r =((double) rand()/(RAND_MAX))*2-1.0;
					
					// Current value to replace if not accepted
					curr[j] = next->x[j];
					
					// Finding maximum distance allowed to walk
					double maxDist = fmax(fabs(this->param->xMax[j]-curr[j]),fabs(this->param->xMin[j]-curr[j]));

					double xNew = 0.0;
					for(k=0; k<100; k++) {
						xNew = curr[j] + maxDist*r;

						if((xNew < this->param->xMin[j]) || (xNew > this->param->xMax[j])) {
							maxDist /= 2;
						}
						else {
							break;
						}
					}
					
					// New guess
					next->x[j] = xNew; 
				}
			}

			// Evaluate functional
			next->T = T;
			evaluate(next);
			this->param->nFuncEval++;

			// Find difference
			fNext = next->f;
			double fDiff = fNext - fCurr;

			// Compute probability
			double P = expf(-fDiff/T);

			// Make decision
			if(fDiff <= 0.0) {
				progress(this, next);
				this->f[this->iter] = next->f;
				this->iter++;
				fCurr = fNext;		
				l++;
				
				// Copy to the returning instance
				opt_vector_copy(next->x,current->x,next->n);
				current->f = fCurr;
				break;
			}
			else {
				// Draw a random number
				srand(time(NULL)); // Random number
				double r = ((double) rand()/(RAND_MAX));

				if(P > r) {
					// We accept the guess
					progress(this, next);
					this->f[this->iter] = next->f;
					this->iter++;
					fCurr = fNext;		
					l++;

					// Copy to the returning instance
					opt_vector_copy(next->x,current->x,next->n);
					current->f = fCurr;
					break;
				}
				else {
					// Do not accept guess
					for(j=0; j<(this->param->n); j++) {
						next->x[j] = curr[j]; // Restoring value
					}
				}
			}

			// Checking if we should stop
			status = opt_check_optimization_sa(next);
			if(status != 0) {
				opt_copy_instance_to_result(this->result,current);
				break;
			}
		}
		
		// Updating temperature
		T = T*dec;
	}

	if(i >= this->param->max_iterations) {
		this->status = OPT_STATUS_MAXITERATIONS;
		opt_copy_instance_to_result(this->result,current);
	}


	// Setting status message
	opt_set_status_msg();

	// Clearing results
	opt_instance_free(next);
	opt_instance_free(current);
	free(curr);
}

}
