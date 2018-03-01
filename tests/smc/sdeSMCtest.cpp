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

// better test this before we get too far...
//[[Rcpp::export("pf_test")]]
int Test() {
  sdeParticle<eou::sdeModel> particle;
  return particle.get_nDims();
}

// particle filter without SMCTC
// requires precomputed normal draws, i.e., deterministic output
//[[Rcpp::export("pf_update")]]
List particleUpdate(NumericVector initParams, NumericMatrix initData,
		    NumericVector dT, IntegerVector nDimsPerObs,
		    NumericMatrix NormalDraws) {
  typedef eou::sdeModel sMod;
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  int nPart = NormalDraws.nrow()/nDims;
  int nStride = nDims*nPart;
  NumericMatrix dataOut(nDims*nPart, nComp);
  NumericMatrix LogWeightOut(nPart, nComp);
  double *yOut = REAL(dataOut);
  double *lwgt = REAL(LogWeightOut);
  sdeFilter<sMod> pf(nPart, nComp, REAL(dT), INTEGER(nDimsPerObs),
		     REAL(initData), REAL(initParams), REAL(NormalDraws));
  // initialize
  int ii,jj;
  for(jj=0; jj<nPart; jj++) {
    lwgt[jj] = pf.init(&yOut[nDims*jj]);
  }
  for(ii=1; ii<nComp; ii++) {
    for(jj=0; jj<nPart; jj++) {
      lwgt[ii*nPart + jj] = pf.update(&yOut[ii*nDims*nPart + jj*nDims],
				      &yOut[(ii-1)*nDims*nPart + jj*nDims],
				      ii, jj, 0);
    }
  }
  return List::create(Rcpp::Named("X") = dataOut,
		      Rcpp::Named("lwgt") = LogWeightOut);
}

// experiment with "algorithm" parameter to SMC constructor
template <class sMod>
class testAlg {
private:
  static const int nDims = sMod::nDims;
public:
  std::vector<double> yObs;
  testAlg() {
    yObs.resize(nDims);
    for(int ii=0; ii<nDims; ii++) {
      yObs[ii] = 1.0*ii;
    }
  }
  ~testAlg() {
  }
};


// -----------------------------------------------------------------------------

// initialize and move methods for Sampler.
// note that the actual computations are done with the sdeFilter object,
// i.e., the algorithm parameters, which efficiently handles memory allocation.

template <class sMod>
void fInitialise(sdeParticle<sMod>& value, double& logweight,
		 sdeFilter<sMod> & pf_calcs) {
  Rprintf("made it to fInitialise.\n");
  // set particle
  logweight = pf_calcs.init(value.yObs);
  pf_calcs.increase_counter();
  return;
}

template <class sMod>
void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
	   sdeFilter<sMod> & pf_calcs) {
  int iCore = 0;
  Rprintf("counter = %i\n", pf_calcs.get_counter());
  logweight += pf_calcs.update(value.yObs, value.yObs, lTime,
			       pf_calcs.get_counter(), iCore);
  pf_calcs.increase_counter();
  return;
}


//[[Rcpp::export("pf_eval_test")]]
int particleEvalTest(NumericVector initParams, NumericMatrix initData,
		     NumericVector dT, IntegerVector nDimsPerObs,
		     NumericMatrix NormalDraws) {
  typedef eou::sdeModel sMod;
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  int nPart = NormalDraws.nrow()/nDims;
  int nStride = nDims*nPart;
  NumericMatrix dataOut(nDims*nPart, nComp);
  NumericMatrix LogWeightOut(nPart, nComp);
  double *yOut = REAL(dataOut);
  double *lwgt = REAL(LogWeightOut);
  sdeParticle<sMod> pOut;
  Rprintf("before SMC.\n");
  // SMC
  try {
    smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> >
      Sampler((long)nPart, HistoryType::NONE);
    smc::moveset<sdeParticle<sMod>, sdeFilter<sMod> >
      Moveset(fInitialise<sMod>, fMove<sMod>, NULL);

    sdeFilter<sMod> sf;
    Rprintf("created empty sdeFilter.\n");
    sdeFilter<sMod> sf2(nPart, nComp, REAL(dT),
			INTEGER(nDimsPerObs),
			REAL(initData), REAL(initParams),
			REAL(NormalDraws));
    Rprintf("created full sdeFilter.\n");
    sf = sf2;
    Rprintf("assigned sdeFilter.\n");
    // Sampler.SetAlgParam(sdeFilter<sMod>(nPart, nComp, REAL(dT),
    // 					INTEGER(nDimsPerObs),
    // 					REAL(initData), REAL(initParams),
    // 					REAL(NormalDraws)));
    return 0;
  }
  catch(smc::exception e) {
    Rcpp::Rcout << e;
  }
  return 0;
}
