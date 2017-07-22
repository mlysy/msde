/*

OK interface considerations for the implementation.

1. sde.make.model should create an Xptr instead of returning R/C++ entrypoints
   for all functions.  Therefore, it should trigger the construction of an
   sdeInterface object of which the members interface between R/C++ for the
   necessary functions: drift, diff, sim, post, etc.
2. must be able to hold different models in one R session.  This means that Xptr
   should be a templated class.  That is, each class/function making use of sdeModel and sdePrior should accept templates of these.
3. should have pre-compiled versions of some models.  For this, the src of the package should have pointers to different models in its shared object...

OK:

Abstract base, template derived class
sdeCobj : sdeRobj<sdeModel, sdePrior>

R constructor:
sdeCobj *sde = new sdeRobj<sdeModel, sdePrior>;
XPtr<sdeRobj> sdeptr(sde, true);
return sdeptr;

R methods:
sde.drift/diff
sde.loglik
sde.prior
sde.valid.params/data
sde.sim
sde.post

*/

#ifndef sdeInterface_h
#define sdeInterface_h

#include <Rcpp.h>
using namespace Rcpp;
//#include <sdeMCMC.h>
//#include <mcmcUtils.h>
//#include "sdeModel.h"

// Base class for R object
class sdeCobj {
 public:
  virtual int get_nDims() = 0;
  virtual int get_nParams() = 0;
  /*
  virtual LogicalVector isData(NumericVector xIn, NumericVector thetaIn,
			       bool singleX, bool singleTheta, int nReps) = 0;
  virtual LogicalVector isParams(NumericVector thetaIn, int nReps) = 0;
  virtual NumericVector Drift(NumericVector xIn, NumericVector thetaIn,
			      bool singleX, bool singleTheta, int nReps) = 0;
  virtual NumericVector Diff(NumericVector xIn, NumericVector thetaIn,
			     bool singleX, bool singleTheta, int nReps) = 0;
  virtual NumericVector LogLik(NumericVector xIn, NumericVector dTIn,
			       NumericVector thetaIn,
			       int nComp, int nReps,
			       bool singleX, bool singleTheta, int nCores) = 0;
  virtual NumericVector Prior(NumericVector thetaIn, NumericVector xIn,
			      bool singleTheta, bool singleX,
			      int nReps, List phiIn) = 0;
  virtual List Sim(int nDataOut, int N, int burn, int reps, int r, double dT,
		   int MAXBAD, NumericVector initData, NumericVector params,
		   bool singleX, bool singleTheta) = 0;
  virtual List Post(NumericVector initParams, NumericVector initData,
		    NumericVector dT, IntegerVector nDimsPerObs,
		    LogicalVector fixedParams, int nSamples, int burn,
		    int nParamsOut, int nDataOut, IntegerVector dataOutSmp,
		    IntegerVector dataOutComp, IntegerVector dataOutDims,
		    double updateParams, double updateData, List priorArgs,
		    List tunePar, int updateLogLik, int nLogLikOut,
		    int updateLastMiss, int nLastMissOut, int nCores) = 0;
  */
  virtual ~sdeCobj() = 0;
};

// Derived class template for R functions
template <class sMod, class sPi>
  class sdeRobj : public sdeCobj {
 public:
  virtual int get_nDims();
  virtual int get_nParams();
  /*
  virtual LogicalVector isData(NumericVector xIn, NumericVector thetaIn,
			       bool singleX, bool singleTheta, int nReps);
  virtual LogicalVector isParams(NumericVector thetaIn, int nReps);
  virtual NumericVector Drift(NumericVector xIn, NumericVector thetaIn,
			      bool singleX, bool singleTheta, int nReps);
  virtual NumericVector Diff(NumericVector xIn, NumericVector thetaIn,
			     bool singleX, bool singleTheta, int nReps);
  virtual NumericVector LogLik(NumericVector xIn, NumericVector dTIn,
			       NumericVector thetaIn,
			       int nComp, int nReps,
			       bool singleX, bool singleTheta, int nCores);
  virtual NumericVector Prior(NumericVector thetaIn, NumericVector xIn,
			      bool singleTheta, bool singleX,
			      int nReps, List phiIn);
  virtual List Sim(int nDataOut, int N, int burn, int reps, int r, double dT,
		   int MAXBAD, NumericVector initData, NumericVector params,
		   bool singleX, bool singleTheta);
  virtual List Post(NumericVector initParams, NumericVector initData,
		    NumericVector dT, IntegerVector nDimsPerObs,
		    LogicalVector fixedParams, int nSamples, int burn,
		    int nParamsOut, int nDataOut, IntegerVector dataOutSmp,
		    IntegerVector dataOutComp, IntegerVector dataOutDims,
		    double updateParams, double updateData, List priorArgs,
		    List tunePar, int updateLogLik, int nLogLikOut,
		    int updateLastMiss, int nLastMissOut, int nCores);
  */
  virtual ~sdeRobj() {}; // since base has virtual destructor
};

#include <sdeRUtils.h>
/* #include <sdeSim.h> */
/* #include <sdePost.h> */


#endif
