#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::depends("msde")]]
#include <sdeMCMC.h>
#include <mcmcUtils.h>

// NOTE: doesn't output initial data value to simplify algorithm
// easier to do this at R level
//[[Rcpp::export("sde.model$sim")]]
List sdeEulerSim(int nDataOut,
		 int N, int burn, int reps, int r, double dT,
		 int MAXBAD,
		 NumericVector initData, NumericVector params,
		 bool singleX, bool singleTheta) {
  RNGScope scope;

  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double sqrtDT = sqrt(dT);
  int bad = 0;
  // output
  NumericVector dataOut(nDataOut);
  int nBadDraws;

  // storage
  sdeModel *sde = new sdeModel; // sde model
  double *mean = new double[nDims]; // mean
  double *sd = new double[nDims*nDims]; // cholesky factor
  double *X = new double[nDims]; // current value
  double *tmpX = new double[nDims]; // proposed value
  double *Z = new double[nDims]; // random draw
  double *theta; // pointer to parameter
  int ii,jj,kk;

  for(ii = 0; ii < reps; ii++) {
    // initialize chains
    for(kk = 0; kk < nDims; kk++) {
      X[kk] = initData[ii*(!singleX)*nDims + kk];
    }
    theta = &params[ii*(!singleTheta)*nParams];
    // loop through
    for(jj = -burn*r; jj < N*r; jj++) {
      // repeatedly draw from Euler until proposal is valid
      mvEuler(mean, sd, X, dT, sqrtDT, theta, sde);
      do {
	for(kk = 0; kk < nDims; kk++) {
	  Z[kk] = norm_rand();
	}
	xmvn(X, Z, mean, sd, nDims);
	// validate draw
      } while(!sde->isValidData(X, theta) && bad++ < MAXBAD);
      if (bad == MAXBAD) {
	goto stop;
      }
      // store
      if(jj >= 0 && (jj+1) % r == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  dataOut[ii*N*nDims + (jj/r)*nDims + kk] = X[kk];
	}
      }
    }
  }

 stop:
  nBadDraws = bad;

  delete [] X;
  delete [] tmpX;
  delete [] Z;
  delete [] mean;
  delete [] sd;
  delete sde;

  return List::create(_["dataOut"] = dataOut, _["nBadDraws"] = nBadDraws);
}
