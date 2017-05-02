/////////////////////////////////////////////////////

// Euler Simulation of Multivariate SDEs

/////////////////////////////////////////////////////

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("msdeHeaders")]]
#include <sdeLogLik.h>

//[[Rcpp::export("sde.model$sim")]]
List sdeEulerSim(int nDataOut,
		 int N, int reps, int r, double dT,
		 int MAXBAD, NumericVector initData, NumericVector params) {
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
  double *Z = new double[nDims*nCores]; // random draw

  int ii,jj,kk;
  for(ii = 0; ii < reps; ii++) {
    for(jj = 0; jj < N*r; jj++) {
      // initialize chains
      if(jj == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = initData[ii*nDims + kk];
	}
      }
      else {
	mvEuler(mean, sd, X, dT, sqrtDT, &params[ii*nParams], sde);
	// repeatedly draw from Euler until proposal is valid
	do {
	  for(kk = 0; kk < nDims; kk++) {
	    Z[kk] = norm_rand();
	  }
	  xmvn(tmpX, Z, mean, sd, nDims);
	  // validate draw
	} while(!sdeModel::isValidData(tmpX, params) && bad++ < MAXBAD);
	if (bad == MAXBAD) {
	  goto stop;
	}
	for(kk = 0; kk < nDims; kk++) {
	  X[kk] = tmpX[kk];
	}
      }
      // store
      if(jj % r == 0) {
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
