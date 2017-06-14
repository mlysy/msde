////////////////////////////////////////////////////////////////////////////

#ifndef Prior_h
#define Prior_h 1

#include <Rcpp.h>
using namespace Rcpp;
#include "mvnUtils.h"
#include "sdeModel.h"

// prior classes
//
// there are four types of priors:
// (1) flat
// (2) normal
// (3) gaussian copula
// (4) custom
// the priors are constructed from a Prior base class,
// from which virtual method
// logPrior(theta, x0)
// is defined in each of the derived classes.  The prior is on theta and the
// first observation x0, of which some components may be latent.
//
// there are nParams parameters and nDims components of x0,
// of which nMiss0 are unobserved.
//
// flat prior: do nothing
//
// normal prior: tuning parameters are mean and variance.  the dimension of
// these is 0 <= nProperRV <= nDims + nTheta.  By default, this corresponds
// to the _last_ nProperRV components of (x, theta), with a flat prior
// assumed on the remaining components.
//
// gaussian copula prior: more tuning parameters, but same defaults
// regarding dimensions.  untested, since overhauling gcop interface.
//
// custom prior: constructor takes argument "Rcpp::List priorParams"
// from which named elements can be extracted and nRV = nMiss0 + nParams.
// currently unimplemented.

//base prior class
class Prior {
 protected:
  static const int nDims = sdeModel::nDims;
  static const int nParams = sdeModel::nParams;
 public:
  // the type of prior
  enum Type {
    Flat = 1,
    Normal = 2,
    GCop = 3,
    Custom = 4
  };

  //pure virtual function must be defined in each derived class
  virtual double logPrior(double *theta, double *x)=0;
};

// flat prior
class FlatPrior : public Prior {
 public:
  FlatPrior();
  ~FlatPrior();
  double logPrior(double *theta, double *x);
};

// normal prior
class NormalPrior : public Prior {
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
  int nActiveRV, nSkippedRV, nTotalRV;
 public:
  NormalPrior(List priorParams, int nRv);
  ~NormalPrior();
  double logPrior(double *theta, double *x);
};

// gaussian copula prior
class GCopPrior : public Prior {
  // gaussian copula parameters
  int *nBreaks;
  double *range, *dx, *pdf, *lpdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  double *tmpX, *tmpQ;
  int nActiveRV, nSkippedRV, nTotalRV;
  double lgcop(double *x, double *qNorm, int *nBreaks, double *range,
	       double *dx, double *pdf, double *lpdf, double *cdf,
	       double *mean, double *sd, double *RhoCholSd, int n);
 public:
  GCopPrior(List priorParams, int nRv);
  ~GCopPrior();
  // parameter + unobserved first states gaussian copula prior.
  double logPrior(double *theta, double *x);
};

///////////////////////////////////////

// constructor, destructor, and logPrior for Flat Prior

inline FlatPrior::FlatPrior() {
}
inline FlatPrior::~FlatPrior() {
}
inline double FlatPrior::logPrior(double *theta, double *x) {
  return 0.0;
}


///////////////////////////////////////

// constructor, destructor, and logPrior for Normal Prior

inline NormalPrior::NormalPrior(List priorParams, int nRv) {
  nTotalRV = nDims + nParams;
  nActiveRV = nRv;
  nSkippedRV = nTotalRV - nActiveRV;
  mean = REAL(priorParams["Mu"]);
  cholSd = REAL(priorParams["V"]);
  tmpX = new double[nTotalRV];
  tmpZ = new double[nTotalRV];
}

inline NormalPrior::~NormalPrior() {
  delete [] tmpX;
  delete [] tmpZ;
}

inline double NormalPrior::logPrior(double *theta, double *x) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nParams; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lmvn(&tmpX[nSkippedRV], &tmpZ[nSkippedRV], mean, cholSd, nActiveRV);
  return(lp);
}


///////////////////////////////////////

// constructor, destructor, and logPrior for Gaussian Copula Prior

inline GCopPrior::GCopPrior(List priorParams, int nRv) {
  nTotalRV = nDims + nParams;
  nActiveRV = nRv;
  nSkippedRV = nTotalRV - nActiveRV;
  nBreaks = INTEGER(priorParams["nbreaks"]);
  range = REAL(priorParams["rx"]);
  dx = REAL(priorParams["dx"]);
  pdf = REAL(priorParams["dens.y"]);
  lpdf = REAL(priorParams["ldens.y"]);
  cdf = REAL(priorParams["Dens.y"]);
  mean = REAL(priorParams["mean"]);
  sd = REAL(priorParams["sd"]);
  RhoCholSd = REAL(priorParams["Rho"]);
  tmpX = new double[nTotalRV];
  tmpQ = new double[nTotalRV];
}

inline GCopPrior::~GCopPrior() {
  delete [] tmpX;
  delete [] tmpQ;
}

inline double GCopPrior::logPrior(double *theta, double *x) {
  double lp;
  int ii;
  for(ii = 0; ii < nDims; ii++) {
    tmpX[ii] = x[ii];
  }
  for(ii = 0; ii < nParams; ii++) {
    tmpX[nDims+ii] = theta[ii];
  }
  lp = lgcop(&tmpX[nSkippedRV], &tmpQ[nSkippedRV],
	     nBreaks, range, dx, pdf, lpdf, cdf,
	     mean, sd, RhoCholSd, nActiveRV);
  return(lp);
}

// log-density of a Gaussian copula
// TODO: include pi factor
inline double GCopPrior::lgcop(double *x, double *qNorm, int *nBreaks, double *range,
	     double *dx, double *pdf, double *lpdf, double *cdf,
	     double *mean, double *sd, double *RhoCholSd, int n) {
  int ii, jj, colI, densElt, start;
  double lp = 0.0;
  double tmpSum = 0.0;
  // normal quantiles and marginal components
  start = 0;
  for(ii = 0; ii < n; ii ++) {
    densElt = (int) floor((x[ii]-range[2*ii])/dx[ii]);
    if((densElt >= 0) & (densElt < nBreaks[ii])) {
      lp += lpdf[densElt + start];
      qNorm[ii] = x[ii] - (range[2*ii] + densElt * dx[ii]);
      qNorm[ii] *= pdf[densElt + start];
      qNorm[ii] = Rf_qnorm5(cdf[densElt + start + ii] +  qNorm[ii], 0.0, 1.0, 1, 0);
    }
    else {
      lp += Rf_dnorm4(x[ii], mean[ii], sd[ii], 1);
      qNorm[ii] = (x[ii] - mean[ii])/sd[ii];
    }
    start += nBreaks[ii];
  }
  // copula components
  // iid standard normal densities
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  lp += 0.5 * tmpSum;
  // multivariate normal density
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) {
      tmpSum += RhoCholSd[colI + jj] * qNorm[jj];
    }
    qNorm[ii] = (qNorm[ii] - tmpSum)/RhoCholSd[colI + ii];
  }
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
  }
  tmpSum *= 0.5;
  for(ii = 0; ii < n; ii++) {
    tmpSum += log(RhoCholSd[n*ii + ii]);
  }
  lp -= tmpSum;
  return(lp);
}


#endif
