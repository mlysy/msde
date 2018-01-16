//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
// #include <RcppArmadillo.h>
// #include <smctc.h>
// #include <sdeUtils.h>
#include "sdeSMC.h"
namespace eou {
#include "eouModel.h"
}
using namespace Rcpp;

// typedef Rcpp::LogicalVector Logical;
// typedef Rcpp::NumericVector Numeric;
// typedef Rcpp::NumericMatrix Matrix;
// typedef Rcpp::IntegerVector Integer;
// typedef Rcpp::List List;

// particle filter with SMCTC
// requires precomputed normal draws, i.e., deterministic output

// initialize and move methods for Sampler.
// note that the actual computations are done with the sdeFilter object,
// i.e., the algorithm parameters, which efficiently handles memory allocation.

template <class sMod>
void fInitialise(sdeParticle<sMod>& value, double& logweight,
		 sdeFilter<sMod> & pf_calcs) {
  //Rprintf("made it to fInitialise.\n");
  // set particle
  logweight = pf_calcs.init(value.yObs);
  pf_calcs.increase_counter();
  return;
}

template <class sMod>
void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
	   sdeFilter<sMod> & pf_calcs) {
  int iCore = 0;
  //Rprintf("counter = %i\n", pf_calcs.get_counter());
  logweight += pf_calcs.update(value.yObs, value.yObs, lTime,
			       pf_calcs.get_counter(), iCore);
  pf_calcs.increase_counter();
  return;
}

// extract double* representation of each particle from the sdeParticle representation used in the Sampler, i.e., the final output of the PF is a vector of logweights, and a nPart x nDims matrix, where each row is a samples from p(X_T, Y_T | X_0:T, theta).
template <class sMod>
void save_state(double *yOut, double *lwgt,
		smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> > & Sampler,
		sdeParticle<sMod> & pTmp) {
  for(long iPart=0; iPart<Sampler.GetNumber(); iPart++) {
    pTmp = Sampler.GetParticleValueN(iPart);
    pTmp.get_yObs(&yOut[iPart*sMod::nDims]);
    lwgt[iPart] = Sampler.GetParticleLogWeightN(iPart);
  }
  return;
}

// ok next test: use SMCTC to create a log-likelihood evaluator.
// returns a list with the log-weight and updated lastmiss.
// for now same signature as pf_update, but this will be changed as we move
// out of debugging phase.
//[[Rcpp::export(".pf_eval")]]
List particleEval(NumericVector initParams, NumericMatrix initData,
		  NumericVector dT, IntegerVector nDimsPerObs,
		  NumericMatrix NormalDraws) {
  typedef eou::sdeModel sMod;
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  int nPart = NormalDraws.nrow()/nDims;
  int nStride = nDims*nPart;
  // for debugging purposes, output whole history
  NumericMatrix dataOut(nDims*nPart, nComp);
  NumericMatrix LogWeightOut(nPart, nComp);
  // pointers to Rcpp memory, i.e., double* representation
  double *yOut = REAL(dataOut); 
  double *lwgt = REAL(LogWeightOut);
  // Sampler can only deep-copy particles out from its storage
  sdeParticle<sMod> pTmp;
  // this is for eventual parallel implementation.
  smc::adaptMethods<sdeParticle<sMod>, sdeFilter<sMod> > *Adapt;
  Adapt = new sdeAdapt<sMod>;
  //Rprintf("before SMC.\n");
  // SMC
  try {
    smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> >
      Sampler((long)nPart, HistoryType::NONE);
    smc::moveset<sdeParticle<sMod>, sdeFilter<sMod> >
      Moveset(fInitialise<sMod>, fMove<sMod>, NULL);
    //Rprintf("Sampler and Moveset created.\n");

    // AlgParam needs to be deep-copied into Sampler
    //Rprintf("right before SetAlgParams.\n");
    Sampler.SetAlgParam(sdeFilter<sMod>(nPart, nComp, REAL(dT),
					INTEGER(nDimsPerObs),
					REAL(initData), REAL(initParams),
					REAL(NormalDraws)));
    //Rprintf("algParams passed in.\n");
    // no resampling to give same output as R code (for debugging only)
    Sampler.SetResampleParams(ResampleType::RESIDUAL, -1.0);
    Sampler.SetMoveSet(Moveset);
    Sampler.Initialise();
    // extract particle from Sampler
    save_state<sMod>(yOut, lwgt, Sampler, pTmp); 
    Sampler.SetAdaptMethods(Adapt); // TBD for parallel processing
    //Rprintf("Sampler initialized.\n");
    // long lTime;
    for(int ii=1; ii<nComp; ii++) {
      //Rprintf("lTime = %i\n", ii);
      Sampler.Iterate();
      save_state<sMod>(&yOut[ii*(nPart*nDims)],&lwgt[ii*nPart], Sampler, pTmp);
    }
    //delete Sampler;
    delete Adapt;
    return List::create(Rcpp::Named("X") = dataOut,
			Rcpp::Named("lwgt") = LogWeightOut);
  }
  catch(smc::exception e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue;
}
