#ifndef sdePrior_h
#define sdePrior_h 1

#include "sdeModel.h"

// Multivariate normal prior
// three parameters: nProperRV, mean, cholSD
// 0 <= nProperRV <= nDims + nParams
// mvn prior on _last_ nProperRV elements of (x, theta)
// flat prior on remaining elements
// R end should contain basic name and dimension checks

// this is specific to mvnPrior
#include "mvnUtils.h"

class Prior {
 private:
  static const int nDims = sdeModel::nDims;
  static const int nParams = sdeModel::nParams;
  static const int nTotalRV = nDims + nParams;
  int nActiveRV, nSkippedRV;
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
 public:
  double logPrior(double *theta, double *x);
  Prior(double **phi, int nArgs, int *nEachArg, bool *fixedParams, int nMiss0);
  ~Prior();
};

inline Prior::Prior(double **phi, int nArgs, int *nEachArg,
		    bool *fixedParams, int nMiss0) {
  nActiveRV = nEachArg[0];
  nSkippedRV = nTotalRV-nActiveRV;
  if(nActiveRV > 0) {
    mean = new double[nActiveRV];
    cholSd = new double[nActiveRV*nActiveRV];
    tmpX = new double[nTotalRV];
    tmpZ = new double[nTotalRV];
    for(int ii=0; ii<nActiveRV; ii++) {
      mean[ii] = phi[0][ii];
    }
    for(int ii=0; ii<nActiveRV*nActiveRV; ii++) {
      cholSd[ii] = phi[1][ii];
    }
  }
}

inline Prior::~Prior() {
  if(nActiveRV > 0) {
    delete [] mean;
    delete [] cholSd;
    delete [] tmpX;
    delete [] tmpZ;
  }
}

inline double Prior::logPrior(double *theta, double *x) {
  if(nActiveRV == 0) return(0.0);
  double lp;
  int ii;
  for(ii=0; ii<nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii=0; ii<nParams; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lmvn(&tmpX[nSkippedRV], &tmpZ[nSkippedRV], mean, cholSd, nActiveRV);
  return(lp);
}

#endif
