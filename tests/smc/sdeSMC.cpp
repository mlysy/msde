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

/*
template <class sMod>
void fInitialise(sdeParticle<sMod>& value, double& logweight,
		 testAlg<sMod> & pf_calcs) {
  //Rprintf("made it to fInitialise.\n");
  logweight = 0.0;
  //for(int ii=0; ii<sMod::nDims; ii++) {
    //Rprintf("value.yObs.size() = %i\n", value.yObs.size());
    //Rprintf("value.yObs[%i] = %f\n", ii, value.yObs[ii]);
    //value.yObs[ii] = ii*1.0;
  //}
  value.set_yObs(&pf_calcs.yObs[0]);
  Rprintf("made it through fInitialize.\n");
  return;
}

template <class sMod>
void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
	   testAlg<sMod> & pf_calcs) {
  int iCore = 0;
  logweight = lTime * 1.0;
  for(int ii=0; ii<sMod::nDims; ii++) {
    Rprintf("value.yObs[%i] = %f\n", ii, value.yObs[ii]);
    //value.yObs[ii] = lTime*ii*1.0;
  }
  return;
}
*/

// -----------------------------------------------------------------------------

// same thing, but using SMCTC

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
//[[Rcpp::export("pf_eval")]]
List particleEval(NumericVector initParams, NumericMatrix initData,
		  NumericVector dT, IntegerVector nDimsPerObs,
		  NumericMatrix NormalDraws) {
  typedef eou::sdeModel sMod;
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  int nPart = NormalDraws.nrow()/nDims;
  int nStride = nDims*nPart;
  NumericMatrix dataOut(nDims*nPart, nComp);
  NumericMatrix LogWeightOut(nPart, nComp);
  NumericMatrix unLogWeightOut(nPart, nComp);
  double *yOut = REAL(dataOut);
  double *lwgt = REAL(LogWeightOut);
  double *unlwgt = REAL(unLogWeightOut);
  sdeParticle<sMod> pTmp;
  smc::adaptMethods<sdeParticle<sMod>, sdeFilter<sMod> > *Adapt;
  Adapt = new sdeAdapt<sMod>;
  Rprintf("before SMC.\n");
  // SMC
  try {
    smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> >
      Sampler((long)nPart, HistoryType::NONE);
    smc::moveset<sdeParticle<sMod>, sdeFilter<sMod> >
      Moveset(fInitialise<sMod>, fMove<sMod>, NULL);
    Rprintf("Sampler and Moveset created.\n");

    Rprintf("right before SetAlgParams.\n");
    Sampler.SetAlgParam(sdeFilter<sMod>(nPart, nComp, REAL(dT),
					INTEGER(nDimsPerObs),
					REAL(initData), REAL(initParams),
					REAL(NormalDraws)));
    Rprintf("algParams passed in.\n");
    // no resampling to give same output as R code (for debugging only)
    Sampler.SetResampleParams(ResampleType::RESIDUAL, -1.0);
    Sampler.SetMoveSet(Moveset);
    Sampler.Initialise();
    save_state<sMod>(yOut, lwgt, Sampler, pTmp);
    Sampler.SetAdaptMethods(Adapt);
    Rprintf("Sampler initialized.\n");
    // long lTime;
    for(int ii=1; ii<nComp; ii++) {
      Rprintf("lTime = %i\n", ii);
      Sampler.Iterate();
      // to get unormalized log-weights from smctc member function
      for(int jj=0; jj != nPart; ++jj) {
          int tmp = jj + (ii-1)*nPart;
          unlwgt[tmp] = Sampler.GetParticleLogWeightN(jj);
      }
      save_state<sMod>(&yOut[ii*(nPart*nDims)],&lwgt[ii*nPart], Sampler, pTmp);
    }
    //delete Sampler;
    delete Adapt;
    return List::create(Rcpp::Named("X") = dataOut,
			Rcpp::Named("lwgt") = LogWeightOut,
			Rcpp::Named("unlwgt") = unLogWeightOut);
  }
  catch(smc::exception e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue;
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
