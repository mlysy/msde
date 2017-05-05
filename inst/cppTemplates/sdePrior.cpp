#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("msdeHeaders")]]
#include <mvnPrior.h>

//[[Rcpp::export("logprior")]]
NumericVector sde_Prior(NumericVector thetaIn, NumericVector xIn,
			int nReps, List phiIn) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  int ii;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  int nArgs = phiIn.length();
  double **phi = new double*[nArgs];
  int *nEachArg = new int[nArgs];
  for(ii=0; ii<nArgs; ii++) {
    if(Rf_isNull(phiIn[ii])) {
      nEachArg[ii] = 0;
    } else {
      nEachArg[ii] = as<NumericVector>(phiIn[ii]).length();
      phi[ii] = REAL(phiIn[ii]);
    }
  }
  bool *fixedParams;
  int nMiss0 = 0;
  Prior prior(phi, nArgs, nEachArg, fixedParams, nMiss0);
  NumericVector lpOut(nReps);
  double *lp = REAL(lpOut);
  // NOTE: this can't be parallelized because private storage is common
  // to parallelize need array of Prior objects
  for(ii=0; ii<nReps; ii++) {
    lp[ii] = prior.logPrior(&theta[ii*nParams], &x[ii*nDims]);
  }
  delete [] phi;
  delete [] nEachArg;
  return lpOut;
}
