#ifndef eouSMC
#define eouSMC 1

#include <smctc.hh>
#include "eouModel.h"

// Define a template class to create SDE partilce object
// or we can define a derived class based on eouModel
template <class sMod> class sdeParticle {
	private:
  		static const int nDims = sMod::nDims;
  	public:
  		double yObs; // observations of SDE
  		double xObs; // unobservable states

  			
};

// Declare the loglikelihood function at each iteration time lTime
// If the proposal is f(x_n | x_n-1) then the likelihood here should be g(y_n | x_n)
double loglikelihood(long lTime, const sdeParticle<sMod> &X);

// Declare the initialization function
smc::particle<sdeParticle<sMod> > fInitialise(smc::rng *pRng);

// Declare fSelect to select different proposal functions
// We don't need this for eouModel
long fSelect(long lTime, const smc::particle<sdeParticle<sMod> > &p, smc::rng *pRng);

// Declare the proposal function at each iteration time lTime
// By Euler approximation, we can just use the normal density f(x_n | x_n-1)
void fMove(long lTime, smc::particle<sdeParticle<sMod> > &pFrom, smc::rng *pRng);

// MCMC move proposal function
// return 0 if a move is rejected; return 1 if a move is accepted
int fMCMC(long lTime, smc:particle<sdeParticle<sMod> > &pFrom, smc::rng *pRng);

#endif
