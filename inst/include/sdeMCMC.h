#ifndef sdeMCMC_h
#define sdeMCMC_h 1

// posterior inference for msde's

// for speed considerations, and since it is unlikely that users will
// be building their own MCMC routines, all transitions will be members
// of the same object.
// the object will inherit from sdeLogLik.

#include "mvnUtils.h"
#include "sdeLogLik.h"
#include "sdePrior.h"

class sdeMCMC : public sdeLogLik {
  int *missInd;
  int nMiss, nMiss0;
  Prior *prior;
  void mvEraker(double *mean, double *sd,
		double *x0, double *x2,
		double b, double b2, double *theta,
		sdeModel *sde);
 public:
  //int nComp, nDims, nParams; coming from sdeLogLik
  double *currFull, *propFull, *propAccept;
  double *currX, *propX, *currTheta, *propTheta;
  //double *dT, *sqrtDT; coming from sdeLogLik
  double *B, *sqrtB;
  int *nObsComp;
  bool *fixedTheta;
  // propMV **mvX; coming from sdeLogLik
  //propMV *mvTheta; don't need
  //sdeModel *sde; coming from sdeLogLik
  void missGibbsUpdate(double *jumpSd, int *gibbsAccept, int *paramAccept);
  void paramVanillaUpdate(double *jumpSd, int *paramAccept);
  sdeMCMC(int n, double *dt, double *xInit, double *thetaInit,
	  int *xIndex, int *thetaIndex,
	  double **phi, int nArgs, int *nEachArg);
  ~sdeMCMC();
};

inline sdeMCMC::sdeMCMC(int n, double *dt,
			double *xInit, double *thetaInit,
			int *xIndex, int *thetaIndex, double **phi,
			int nArgs, int *nEachArg) : sdeLogLik(n, dt) {
  int ii, jj;
  // problem dimensions
  //nComp = N;
  //nDims = sdeModel::nDims;
  //nParams = sdeModel::nParams;
  // times
  //dT = dt;
  //sqrtDT = new double[nComp];
  B = new double[nComp];
  sqrtB = new double[nComp];
  for(ii=1; ii<nComp-1; ii++) {
    //sqrtDT[ii] = sqrt(dT[ii]);
    //if(ii > 0) {
    B[ii] = dT[ii]/(dT[ii-1] + dT[ii]);
    sqrtB[ii] = sqrt((1-B[ii]) * dT[ii]);
    //}
  }
  // data
  currFull = new double[nComp*nDims + nParams];
  propFull = new double[nComp*nDims + nParams];
  propAccept = new double[nComp];
  currX = currFull + nParams;
  propX = propFull + nParams;
  nObsComp = new int[nComp];
  //mvX = new propMV*[nComp];
  // initialize
  for(ii=0; ii<nComp; ii++) {
    propAccept[ii] = 0.0;
    nObsComp[ii] = xIndex[ii];
    //mvX[ii] = new propMV(nDims);
    for(jj=0; jj<nDims; jj++) {
      currX[ii*nDims + jj] = xInit[ii*nDims + jj];
      propX[ii*nDims + jj] = currX[ii*nDims + jj];
    }
  }
  // missing data
  int nMiss0 = nDims-nObsComp[0]; // unobserved states in first observation
  int nMissN = nDims-nObsComp[nComp-1]; // unobserved states in last observation
  // identify missing data indices, i.e. at least one component to update
  nMiss = 0;
  for(ii = 0; ii < nComp; ii++) {
    nMiss += (nObsComp[ii] < nDims);
  }
  missInd = new int[nMiss + (nMiss == 0)];
  jj = 0;
  for(ii = 0; ii < nComp; ii++) {
    if(nObsComp[ii] < nDims) {
      missInd[jj++] = ii;
    }
  }
  // parameters
  fixedTheta = new bool[nParams];
  currTheta = currFull;
  propTheta = propFull;
  for(ii=0; ii<nParams; ii++) {
    fixedTheta[ii] = (thetaIndex[ii] != 0); // R logicals stored as C++ ints
    currTheta[ii] = thetaInit[ii];
    propTheta[ii] = currTheta[ii];
  }
  // prior
  // any advantage to passing pointer?
  prior = new Prior(phi, nArgs, nEachArg);
}

inline sdeMCMC::~sdeMCMC() {
  delete [] B;
  delete [] sqrtB;
  delete [] currFull;
  delete [] propFull;
  delete [] propAccept;
  delete [] missInd;
  delete [] nObsComp;
  delete [] fixedTheta;
  delete prior;
}

#include "MissGibbsUpdate.h"
#include "ParamVanillaUpdate.h"

#endif
