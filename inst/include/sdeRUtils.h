#ifndef sdeRUtils_h
#define sdeRUtils_h

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeUtils.h"
#include "sdeLogLik.h"
#include "sdeInterface.h"

template <class sMod, class sPi>
  inline int sdeRobj<sMod, sPi>::get_nDims() {
  return sMod::nDims;
}


template <class sMod, class sPi>
  inline int sdeRobj<sMod, sPi>::get_nParams() {
  return sMod::nParams;
}

template <class sMod, class sPi>
  inline LogicalVector sdeRobj<sMod, sPi>::isData(NumericVector xIn,
						  NumericVector thetaIn,
						  bool singleX,
						  bool singleTheta, int nReps) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  LogicalVector validOut(nReps);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidData(&x[ii*(!singleX)*nDims],
				   &theta[ii*(!singleTheta)*nParams]);
  }
  return validOut;
}

template <class sMod, class sPi>
  inline LogicalVector sdeRobj<sMod, sPi>::isParams(NumericVector thetaIn, int nReps) {
  int nParams = sMod::nParams;
  double *theta = REAL(thetaIn);
  LogicalVector validOut(nReps);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidParams(&theta[ii*nParams]);
  }
  return validOut;
}

// Best to parallelize these from within R: easier and prob faster since
// memory can be allocated in parallel there
template <class sMod, class sPi>
  inline NumericVector sdeRobj<sMod, sPi>::Drift(NumericVector xIn,
						 NumericVector thetaIn,
						 bool singleX,
						 bool singleTheta,
						 int nReps) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    sde.sdeDr(&dr[ii*nDims],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
  }
  return drOut;
}

template <class sMod, class sPi>
  inline NumericVector sdeRobj<sMod, sPi>::Diff(NumericVector xIn,
						NumericVector thetaIn,
						bool singleX,
						bool singleTheta, int nReps) {
  int nDims = sMod::nDims;
  int nDims2 = nDims*nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims2);
  double *df = REAL(dfOut);
  sMod sde;
  int ii,jj;
  for(ii=0; ii<nReps; ii++) {
    sde.sdeDf(&df[ii*nDims2],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
    if(sMod::diagDiff) {
      if(sMod::sdDiff) {
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = df[ii*nDims2+jj];
	}
      } else {
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = sqrt(df[ii*nDims2+jj]);
	}
      }
    } else {
      if(!sMod::sdDiff) {
	chol_decomp(&df[ii*nDims2], &df[ii*nDims2], nDims);
      }
    }
  }
  return dfOut;
}

// SDE log-likelihood evaluation.
template <class sMod, class sPi>
  inline NumericVector sdeRobj<sMod, sPi>::LogLik(NumericVector xIn,
						  NumericVector dTIn,
						  NumericVector thetaIn,
						  int nComp, int nReps,
						  bool singleX,
						  bool singleTheta,
						  int nCores) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  sdeLogLik<sMod> sdeLL(nComp, REAL(dTIn), nCores);
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*(!singleTheta)*nParams],
			  &x[ii*(!singleX)*nDims*nComp]);
  }
  return llOut;
}

template <class sMod, class sPi>
  inline NumericVector sdeRobj<sMod, sPi>::Prior(NumericVector thetaIn,
						 NumericVector xIn,
						 bool singleTheta,
						 bool singleX,
						 int nReps, List phiIn) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  int ii;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  PriorArgs priorArgs(phiIn);
  /* int nArgs = phiIn.length(); */
  /* double **phi = new double*[nArgs]; */
  /* int *nEachArg = new int[nArgs]; */
  /* for(ii=0; ii<nArgs; ii++) { */
  /*   if(Rf_isNull(phiIn[ii])) { */
  /*     nEachArg[ii] = 0; */
  /*   } else { */
  /*     nEachArg[ii] = as<NumericVector>(phiIn[ii]).length(); */
  /*     phi[ii] = REAL(phiIn[ii]); */
  /*   } */
  /* } */
  sPi prior(priorArgs.phi, priorArgs.nArgs, priorArgs.nEachArg);
  NumericVector lpOut(nReps);
  double *lp = REAL(lpOut);
  // NOTE: this can't be parallelized because private storage is common
  // to parallelize need array of sdePrior objects
  for(ii=0; ii<nReps; ii++) {
    lp[ii] = prior.logPrior(&theta[ii*(!singleTheta)*nParams],
			    &x[ii*(!singleX)*nDims]);
  }
  /* delete [] phi; */
  /* delete [] nEachArg; */
  return lpOut;
}

#endif
