#ifndef MissGibbsUpdate_h
#define MissGibbsUpdate_h 1

// componentwise data updates.

#include "mvnUtils.h"
#include "sdeLogLik.h"

class MissGibbs {
 private:
  int nComp, nParams, nDims;
  int *missInd;
  double *propMean, *propSd, *propZ;
  double *B, sqrtB;
  sdeLogLik *data;
 public:
  void mvEraker(double *mean, double *sd, double *x, int iObs);
  void update(double *x, double *theta,
	      double *jumpSd, double *missAccept, double *paramAccept);
  MissGibbs(sdeLogLik *LL, int *xIndex);
  ~MissGibbs;
};

inline MissGibbs::MissGibbs(int *xIndex, sdeLogLik *LL) {
  data = LL;
  nComp = data->nObs;
  nDims = sdeModel::nDims;
  nParams = sdeModel::nParams;
  nObsComp = xIndex; // TODO: copy memory?
  // unobserved states in first observation
  int nMiss0 = nDims-nObsComp[0];
  // unobserved states in last observation
  int nMissN = nDims-nObsComp[nComp-1]; 
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

}

inline void MissGibbs::update(double *x, double *theta, double *jumpSd,
			      double *missAccept, double *paramAccept) {
  int ii, II, jj, JJ;
  int startII, endII;
  //int foundNaN = 0;
  // only elements in missInd are updated.
  // first and last components are handled separately.
  startII = (missInd[0] == 0) ? 1 : 0;
  endII = (missInd[nMiss-1] == nComp-1) ? nMiss-1 : nMiss;
  // Markov chain elements are conditionally independent,
  // so every other can be updated in parallel
  for(JJ = 0; JJ < 2; JJ++) {
    // *** PARALLELIZABLE FOR-LOOP ***
    for(II = startII+JJ; II < endII; II = II+2) {
      ii = missInd[II];
      // intermediate data points
      if((0 < ii) && (ii < nComp-1)) {
	mvEraker(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
		 &currX[(ii+1)*nDims], B[ii], sqrtB[ii], currTheta, &sde[ii]);
	// partial observations
	if(nObsComp[ii] == 0) {
	  for(jj = 0; jj < nDims; jj++) mvX[ii]->z[jj] = norm_rand();
	}
	else {
	  zmvn(mvX[ii]->z, &currX[ii*nDims], mvX[ii]->mean,
	       mvX[ii]->sd, nDims, nObsComp[ii]);
	  for(jj = nObsComp[ii]; jj < nDims; jj++) {
	    mvX[ii]->z[jj] = norm_rand();
	  }
	}
	// proposals
	xmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
	     mvX[ii]->sd, nDims);
	//if(isNaN(&propX[ii*nDims], nDims)) {
	//  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
	//  foundNaN = 1;
	//  goto stop;
	//}
	// only calculate acceptance rate if proposal is valid
	if(sdeModel::isValidData(&propX[ii*nDims])) {
	  // acceptance rate
	  // proposal
	  propAccept[ii] = lmvn(&currX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				mvX[ii]->sd, nDims);
	  propAccept[ii] -= lmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // target 1
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
		  dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  propAccept[ii] -= lmvn(&currX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // target 2
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &propX[ii*nDims], dT[ii],
		  sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[ii*nDims], dT[ii],
		  sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
				 mvX[ii]->sd, nDims);
	  // evaluate mh ratio
	  if(exp(propAccept[ii]) >= unif_rand()) {
	    for(jj = 0; jj < nDims; jj++) {
	      currX[ii*nDims + jj] = propX[ii*nDims + jj];
	    }
	    gibbsAccept[ii]++;
	  }
	}
      }
    }
  }
  // special treatment for last datapoint
  if(missInd[nMiss-1] == nComp-1) {
    ii = nComp-1;
    mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
    // partial observations
    if(nObsComp[ii] == 0) {
      for(jj = 0; jj < nDims; jj++) {
	mvX[ii]->z[jj] = norm_rand();
      }
    }
    else {
      zmvn(mvX[ii]->z, &currX[ii*nDims], mvX[ii]->mean,
	   mvX[ii]->sd, nDims, nObsComp[ii]);
      for(jj = nObsComp[ii]; jj < nDims; jj++) {
	mvX[ii]->z[jj] = norm_rand();
      }
    }
    // proposals
    xmvn(&propX[ii*nDims], mvX[ii]->z, mvX[ii]->mean,
	 mvX[ii]->sd, nDims);
    //if(isNaN(&propX[ii*nDims], nDims)) {
    //  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
    //  foundNaN = 1;
    //  goto stop;
    //}
    // acceptance is 100% as long as the proposal is valid
    if(sdeModel::isValidData(&propX[ii*nDims])) {
      for(jj = 0; jj < nDims; jj++) {
	currX[ii*nDims + jj] = propX[ii*nDims + jj];
      }
      gibbsAccept[ii]++;
    }
  }
  // special treatment for first datapoint
  if(missInd[0] == 0) {
    ii = 0;
    // initialize
    for(jj = 0; jj < nDims; jj++) {
      propX[jj] = currX[jj];
    }
    // random walk metropolis
    for(jj = 0; jj < nMiss0; jj++) {
      // proposal
      propX[nObsComp[0]+jj] = currX[nObsComp[0]+jj] + jumpSd[nParams+jj] * norm_rand();
      if(sdeModel::isValidData(&propX[ii*nDims])) {
	// acceptance rate.
	// target 1
	propAccept[ii] = prior->logPrior(currTheta, propX);
	propAccept[ii] -= prior->logPrior(currTheta, currX);
	// target 2
	mvEuler(mvX[ii]->mean, mvX[ii]->sd, &propX[ii*nDims], dT[ii],
		sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] += lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
			       mvX[ii]->sd, nDims);
	mvEuler(mvX[ii]->mean, mvX[ii]->sd, &currX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], mvX[ii]->z, mvX[ii]->mean,
			       mvX[ii]->sd, nDims);
	// evaluate mh ratio
	if(exp(propAccept[ii]) >= unif_rand()) {
	  currX[ii*nDims + nObsComp[0]+jj] = propX[ii*nDims + nObsComp[0]+jj];
	  paramAccept[nParams + jj]++;
	}
	else {
	  propX[ii*nDims + nObsComp[0]+jj] = currX[ii*nDims + nObsComp[0]+jj];
	}
      }
    }
  }
}

#endif
