// TODO: implement filtering by using sequential monte carlo
// just based on eouModel

//[[Rcpp::depends("Rcpp")]]
//[[Rcpp::plugins("cpp11")]]
//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("msde")]]
// #include <RcppArmadillo.h>
#include <smctc.h>
#include <Rcpp.h>
#include <sdeUtils.h>

#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>

#include "eouSMC.hh"
#include "eouModel.h"

using namespace Rcpp;


// Use SMCTC to create a particle filtering function for SDE models 
// return a list of log weights, values of particles, approximation of log p(y)
// \param obsData The observed data input, a vector or matrix with dimension nObs x nobsVar
// \param Params The SDE parameter settings, we assume parameters are known
// \param dT  A vector of the time interval between two observations
// \param lNumber The number of partilces
// [[Rcpp::export("pf_eval")]]
List sde_pf(NumericMatrix obsData, NumericVector Params, 
			NumericVector dT, long lNumber) {
// Load the observations/data


long lIterates = something; // number of iterations determined by the input observations
int nDims = sMod::nDims; // number of dimensions
int nComp = initData.ncol();
int nPart = lNumber;

NumericMatrix dataOut(nDims*nPart, nComp);
NumericMatrix LogWeightOut(nPart, nComp);
double *yOut = REAL(dataOut);
double *lwgt = REAL(LogWeightOut);

try {
	// Create the smc::sampler<sdeModel> object Sampler
	smc::sampler<particle<sdeModel> > Sampler(lNumber, SMC_HISTORY_NONE);
	// Create the smc::sampler<sdeModel> object Moveset 
	// arguments: initialization, proposal, MCMC move
	smc::moveset<particle<sdeModel> > Moveset(fInitialise, fMove, NULL);

	Rprintf("Sampler and Moveset created.\n");

	// Set resampling strategy
	Sampler.SetResampleParams(SMC_RESAMPLE_MULTINOMIAL, 0.2);
	// Tell the sampler object all the information about 
	// initialization, proposals, weighting and MCMC moves to 
	// move the particle system
	Sampler.SetMoveSet(Moveset);

	// Initialize the sampler
	Sampler.Initialise();

	Rprintf("Sampler initialized.\n");

	// Start the iteration
	Rprintf("Iteration Starts.\n");
	for(auto ii = 1; ii != lIterates; ++ii) {
		Rprintf("lTime = %i\n", ii);
		// Move the particle system to the next iteration
		Sampler.Iterate();
	}
	// get the final log-weights and approximation of marginal density
	return List::create(Rcpp::Named("X") = dataOut, Rcpp::Named("lwgt") = LogWeightOut);
}
catch(smc::exception e) 
  {
	Rcpp::Rcout << e;
	//exit(e.lCode);
  }
}

// Define the loglikelihood function
// \param lTime The current time (i.e. the index of the current distribution)
// \param X     The sdeModel to consider 
double loglikelihood(long lTime, const sdeParticle<sMod> &X) {

}

// Define the initialization function
// \param pRng A pointer to an smc:rng class which serves as random number generators
smc::particle<sdeModel> fInitialise(smc::rng *pRng) {

}

// Define the proposal function
// \param lTime The current time (i.e. the index of the current distribution)
// \param pFrom The reference to the particle to be updated
// \param pRng A pointer to an smc:rng class which serves as random number generators
void fMove(long lTime, smc::particle<sdeParticle<sMod> > &pFrom, smc::rng *pRng) {

}

// Define the MCMC move
// \param lTime The current time (i.e. the index of the current distribution)
// \param pFrom The reference to the particle to be updated
// \param pRng A pointer to an smc:rng class which serves as random number generators
int fMCMC(long lTime, smc:particle<sdeParticle<sMod> > &pFrom, smc::rng *pRng) {

}
