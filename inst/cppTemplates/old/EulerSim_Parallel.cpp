/////////////////////////////////////////////////////

// Euler Simulation of Multivariate SDEs

/////////////////////////////////////////////////////

// this is a multi-core implementation, but not RNG thread-safe.
// probably better to parallelize from R instead.

#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("msdeHeaders")]]
#include <sdeLogLik.h>

//[[Rcpp::export("sde.model$sim")]]
List sdeEulerSim(int nDataOut,
		 int N, int reps, int r,
		 double tStart, double dT,
		 int MAXBAD, NumericVector initData, NumericVector params) {
  RNGScope scope;
  const int nCores = 1; // for future parallel implementation

  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double sqrtDT = sqrt(dT);
  int bad = 0;
  // output
  NumericVector dataOut(nDataOut);
  int nBadDraws;

  // storage
  sdeLogLik data(nCores); // sde mean variance, proposals, etc.
  double *X = new double[nDims*nCores]; // current value
  double *tmpX = new double[nDims*nCores]; // proposed value
  double *Z = new double[nDims*nCores]; // random draw

  // *** PARALLELIZABLE FOR-LOOP ***
  int nBlock = reps/nCores+1; // number of replicates per core
  for(int mm=0; mm<nCores; mm++) {
    int IIstart = mm*nBlock;
    int IIend = (mm+1)*nBlock;
    IIend = IIend < reps ? IIend : reps; // last cores have less
    int ii, jj, kk;
    int II = mm*nDims;
    double currT;
    for(ii=IIstart; ii<IIend; ii++) {
      for(jj=0; jj<N*r; jj++) {
	if(jj == 0) {
	  // initialize chains
	  currT = tStart;
	  for(kk=0; kk<nDims; kk++) {
	    X[II + kk] = initData[ii*nDims + kk];
	  }
	} else {
	  // create draw
	  data.mvEuler(&X[II],
		       currT, dT, sqrtDT, &params[ii*nParams], mm);
	  do {
	    // repeatedly draw until draw is valid
	    for(kk = 0; kk < nDims; kk++) {
	      Z[II + kk] = norm_rand();
	    }
	    xmvn(&tmpX[II], &Z[II],
		 data.getMean(mm), data.getSd(mm), nDims);
	    // validate draw
	  } while(!sdeModel::isValidData(tmpX) && bad++ < MAXBAD);
	  if (bad >= MAXBAD) {
	    goto stop;
	  }
	  for(kk=0; kk<nDims; kk++) {
	    X[II + kk] = tmpX[II + kk];
	  }
	  currT += dT;
	}
	// store
	if(jj % r == 0) {
	  for(kk = 0; kk < nDims; kk++) {
	    dataOut[ii*N*nDims + (jj/r)*nDims + kk] = X[II + kk];
	  }
	}
      }
    }
  }

 stop:
  nBadDraws = bad;

  delete [] X;
  delete [] tmpX;
  delete [] Z;
  delete data;

  return List::create(_["dataOut"] = dataOut, _["nBadDraws"] = nBadDraws);
}
