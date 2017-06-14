#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::depends("msdeHeaders")]]
#include <sdeMCMC.h>
#include <mcmcUtils.h>

// some utilities to be exported, namely for checking the c++ code
//
// the function are: sde's drift, diffusion, log-likelihood, prior.

// //[[Rcpp::export(".get.sde.consts")]]
// List sde_getConsts() {
//   return List::create(_["ndims"] = sdeModel::nDims,
// 		      _["nparams"] = sdeModel::nParams,
// 		      _["sd.diff"] = sdeModel::sdDiff,
// 		      _["diag.diff"] = sdeModel::diagDiff);
// }

//[[Rcpp::export("ndims")]]
int sde_getNDims() {
  return sdeModel::nDims;
}

//[[Rcpp::export("nparams")]]
int sde_getNParams() {
  return sdeModel::nParams;
}

//[[Rcpp::export("sde.model$is.data")]]
LogicalVector sde_isData(NumericVector xIn, NumericVector thetaIn,
			 bool singleX, bool singleTheta, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  LogicalVector validOut(nReps);
  sdeModel sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidData(&x[ii*(!singleX)*nDims],
				   &theta[ii*(!singleTheta)*nParams]);
  }
  return validOut;
}

//[[Rcpp::export("sde.model$is.params")]]
LogicalVector sde_isParams(NumericVector thetaIn, int nReps) {
  int nParams = sdeModel::nParams;
  double *theta = REAL(thetaIn);
  LogicalVector validOut(nReps);
  sdeModel sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidParams(&theta[ii*nParams]);
  }
  return validOut;
}

// Best to parallelize these from within R: easier and prob faster since
// memory can be allocated in parallel there
//[[Rcpp::export("sde.model$drift")]]
NumericVector sde_Drift(NumericVector xIn, NumericVector thetaIn,
			bool singleX, bool singleTheta, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sdeModel sde;
  for(int ii = 0; ii < nReps; ii++) {
    sde.sdeDr(&dr[ii*nDims],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
  }
  return drOut;
}

//[[Rcpp::export("sde.model$diff")]]
NumericVector sde_Diff(NumericVector xIn, NumericVector thetaIn,
		       bool singleX, bool singleTheta, int nReps) {
  int nDims = sdeModel::nDims;
  int nDims2 = nDims*nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims2);
  double *df = REAL(dfOut);
  sdeModel sde;
  int ii,jj;
  for(ii=0; ii<nReps; ii++) {
    //Rprintf("ii = %i\n", ii);
    sde.sdeDf(&df[ii*nDims2],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
    // for(jj=0; jj<nDims2; jj++) {
    //   Rprintf("df[%i] = %f\n", jj, df[ii*nDims2+jj]);
    // }
    if(sdeModel::diagDiff) {
      //Rprintf("diag scale\n");
      if(sdeModel::sdDiff) {
	//Rprintf("sd scale\n");
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = df[ii*nDims2+jj];
	}
      } else {
	//Rprintf("var scale");
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = sqrt(df[ii*nDims2+jj]);
	}
      }
    } else {
      //Rprintf("matrix\n");
      if(!sdeModel::sdDiff) {
	//Rprintf("var\n");
	chol_decomp(&df[ii*nDims2], &df[ii*nDims2], nDims);
      }
    }
  }
  return dfOut;
}

// SDE log-likelihood evaluation.
//[[Rcpp::export("sde.model$loglik")]]
NumericVector sde_LogLik(NumericVector xIn, NumericVector dTIn,
			 NumericVector thetaIn,
			 int nComp, int nReps,
			 bool singleX, bool singleTheta, int nCores) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  sdeLogLik sdeLL(nComp, REAL(dTIn), nCores);
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*(!singleTheta)*nParams],
			  &x[ii*(!singleX)*nDims*nComp]);
  }
  return llOut;
}

//[[Rcpp::export("sde.model$logprior")]]
NumericVector sde_Prior(NumericVector thetaIn, NumericVector xIn,
			bool singleTheta, bool singleX,
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
  Prior prior(phi, nArgs, nEachArg);
  NumericVector lpOut(nReps);
  double *lp = REAL(lpOut);
  // NOTE: this can't be parallelized because private storage is common
  // to parallelize need array of Prior objects
  for(ii=0; ii<nReps; ii++) {
    lp[ii] = prior.logPrior(&theta[ii*(!singleTheta)*nParams],
			    &x[ii*(!singleX)*nDims]);
  }
  delete [] phi;
  delete [] nEachArg;
  return lpOut;
}
