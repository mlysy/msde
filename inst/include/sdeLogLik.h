#ifndef sdeLogLik_h
#define sdeLogLik_h 1

/* -------------------------------------------------------------------------
/*
this is a class which defines the complete data and whatever it needs
to evaluate the log-density.
the log-density of an sde is essentially the sum of conditional normals,
each having its own mean and variance for every point.
i.e.,
logDens(x | theta) = sum dnorm(x_n+1 | mu(x_n, theta), sd(x_n, theta)).

so the constructor should create an array of sdeModels, each with enough
storage to compute log-densities.

public members:
* nComp, tSeq, dT, sqrtDT
* Drift(x, t, theta, i), Diff(x, t, theta, i)
* EulerMV(x, theta, i)
* loglik(x, theta)
*/

#include "LinAlgUtils.h"
#include "mvnUtils.h"
#include "sdeModel.h"
#include "sdeUtils.h"

class sdeLogLik {
 protected:
  int nDims2;
  //void init(int n); // initialize temporary storage
 public:
  double *propMean, *propSd, *propZ; // for storing normal calculations
  sdeModel *sde; // for storing drift and diffusion calculations
  int nDims, nParams; // internal representations taken from sdeModel
  int nComp; // number of observations INCLUDING first one
  double *dT, *sqrtDT; // times
  /* // observation-specific drift and diffusions */
  /* void drift(double *dr, double *x, double *theta, int iObs); */
  /* void diff(double *df, double *x, double *theta, int iObs); */
  /* // observation-specific mean and variance */
  /* // in case time needs to be specified */
  /* void mvEuler(double *x, double t, double dT, double sqrtDT, */
  /* 	       double *theta, int iObs); */
  /* void mvEuler(double *x, double *theta, int iObs); */
  /* // returns pointer */
  /* double* getMean(int iObs) { */
  /*   return &propMean[iObs*nDims]; */
  /* } */
  /* double* getSd(int iObs) { */
  /*   return &propSd[iObs*nDims2]; */
  /* } */
  // log-density
  double loglik(double *theta, double *x);
  // constructor and destructor
  //sdeLogLik(int n);
  sdeLogLik(int n, double *dt);
  ~sdeLogLik();
};


// constructors ------------------------------------------------------------

inline sdeLogLik::sdeLogLik(int n, double *dt) {
  nComp = n;
  nDims = sdeModel::nDims;
  nDims2 = nDims*nDims;
  nParams = sdeModel::nParams;
  // create storage space
  sde = new sdeModel[nComp];
  propMean = new double[nComp*nDims];
  propSd = new double[nComp*nDims*nDims];
  propZ = new double[nComp*nDims];
  dT = new double[nComp];
  sqrtDT = new double[nComp];
  // timing
  for(int ii=0; ii<nComp-1; ii++) {
    dT[ii] = dt[ii];
    sqrtDT[ii] = sqrt(dT[ii]);
  }
}

/* // initialize to avoid c++11 */
/* inline void sdeLogLik::init(int n) { */
/*   nComp = n; */
/*   nDims = sdeModel::nDims; */
/*   nDims2 = nDims*nDims; */
/*   nParams = sdeModel::nParams; */
/*   // create storage space */
/*   sde = new sdeModel[nComp]; */
/*   propMean = new double[nComp*nDims]; */
/*   propSd = new double[nComp*nDims*nDims]; */
/*   propZ = new double[nComp*nDims]; */
/*   tSeq = new double[nComp]; */
/*   dT = new double[nComp]; */
/*   sqrtDT = new double[nComp]; */
/* } */
/* inline sdeLogLik::sdeLogLik(int n) { */
/*   init(n); */
/* } */
/* inline sdeLogLik::sdeLogLik(int n, double *t) { */
/*   init(n); */
/*   // timing */
/*   tSeq[0] = t[0]; */
/*   for(int ii=0; ii<nComp-1; ii++) { */
/*     tSeq[ii+1] = t[ii+1]; */
/*     dT[ii] = tSeq[ii+1]-tSeq[ii]; */
/*     sqrtDT[ii] = sqrt(dT[ii]); */
/*   } */
/* } */

// destructor
inline sdeLogLik::~sdeLogLik() {
  delete [] sde;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] dT;
  delete [] sqrtDT;
}

/* // drift and diffusion functions */
/* inline void sdeLogLik::drift(double *dr, double *x, double *theta, */
/* 			     int iObs) { */
/*   sde[iObs].sdeDr(dr, x, tSeq[iObs], theta); */
/*   return; */
/* } */

/* inline void sdeLogLik::diff(double *df, double *x, double *theta, */
/* 			    int iObs) { */
/*   sde[iObs].sdeDf(df, x, tSeq[iObs], theta); */
/*   return; */
/* } */

/* // mean and variance of Euler step */

/* // general case */
/* inline void sdeLogLik::mvEuler(double *x, */
/* 			       double t, double dT, double sqrtDT, */
/* 			       double *theta, int iObs) { */
/*   // mean = x + drift(x,t,theta)*dT */
/*   int II = iObs*nDims; */
/*   sde[iObs].sdeDr(&propMean[II], x, t, theta); */
/*   v_mult(&propMean[II], dT, nDims); */
/*   for(int jj = 0; jj < nDims; jj++) { */
/*     propMean[II + jj] += x[jj]; */
/*   } */
/*   // sd = diff(x,t,theta)*sqrt(dT) */
/*   II *= nDims; */
/*   sde[iObs].sdeDf(&propSd[II], x, t, theta); */
/*   U_mult(&propSd[II], sqrtDT, nDims); */
/*   return; */
/* } */

/* // time-assigned case */
/* inline void sdeLogLik::mvEuler(double *x, double *theta, int iObs) { */
/*   mvEuler(x, tSeq[iObs], dT[iObs], sqrtDT[iObs], theta, iObs); */
/*   return; */
/* } */

// full log-likelihood evaluation
inline double sdeLogLik::loglik(double *theta, double *x) {
  double ll = 0;
  // *** PARALLELIZABLE FOR-LOOP ***
  for(int ii = 0; ii < nComp-1; ii++) {
    mvEuler(&propMean[ii*nDims], &propSd[ii*nDims2],
	    &x[ii*nDims], dT[ii], sqrtDT[ii], theta, &sde[ii]);
    ll += lmvn(&x[(ii+1)*nDims], &propZ[ii*nDims],
	       &propMean[ii*nDims], &propSd[ii*nDims2], nDims);
  }
  return(ll);
}

#endif
