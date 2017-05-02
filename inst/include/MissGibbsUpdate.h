#ifndef MissGibbsUpdate_h
#define MissGibbsUpdate_h 1

#include "mvnUtils.h"
#include "sdeMCMC.h"

// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
inline void sdeMCMC::mvEraker(double *mean, double *sd,
	      double *x0, double *x2,
	      double b, double b2, double *theta,
	      sdeModel *sde) {
  double b1 = 1.0-b;
  for(int ii=0; ii<sdeModel::nDims; ii++) {
    mean[ii] = x0[ii] * b + x2[ii] * b1;
  }
  sde->sdeDf(sd, x0, theta);
  U_mult(sd, b2, sdeModel::nDims);
  return;
}

inline void sdeMCMC::missGibbsUpdate(double *jumpSd, int *gibbsAccept,
				     int *paramAccept) {
  int ii, II, jj, JJ;
  int startII, endII;
  double *mean, *sd, *Z, *tmpZ;
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
      // these variables will need to be privatized
      mean = &propMean[ii*nDims];
      sd = &propSd[ii*nDims2];
      Z = &propZ[ii*nDims];
      // intermediate data points
      if((0 < ii) && (ii < nComp-1)) {
	mvEraker(mean, sd,
		 &currX[(ii-1)*nDims], &currX[(ii+1)*nDims],
		 B[ii], sqrtB[ii], currTheta,
		 &sde[ii]);
	// partial observations
	if(nObsComp[ii] == 0) {
	  for(jj=0; jj<nDims; jj++) Z[jj] = norm_rand();
	}
	else {
	  zmvn(Z, &currX[ii*nDims], mean, sd, nDims, nObsComp[ii]);
	  for(jj=nObsComp[ii]; jj<nDims; jj++) {
	    Z[jj] = norm_rand();
	  }
	}
	// proposals
	xmvn(&propX[ii*nDims], Z, mean, sd, nDims);
	//if(isNaN(&propX[ii*nDims], nDims)) {
	//  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
	//  foundNaN = 1;
	//  goto stop;
	//}
	// only calculate acceptance rate if proposal is valid
	if(sde[ii].isValidData(&propX[ii*nDims], currTheta)) {
	  // acceptance rate
	  // proposal
	  propAccept[ii] = lmvn(&currX[ii*nDims], Z, mean, sd, nDims);
	  propAccept[ii] -= lmvn(&propX[ii*nDims], Z, mean, sd, nDims);
	  // target 1
	  mvEuler(mean, sd, &currX[(ii-1)*nDims],
		  dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&propX[ii*nDims], Z, mean, sd, nDims);
	  propAccept[ii] -= lmvn(&currX[ii*nDims], Z, mean, sd, nDims);
	  // target 2
	  mvEuler(mean, sd, &propX[ii*nDims],
		  dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] += lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	  mvEuler(mean, sd, &currX[ii*nDims],
		  dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	  propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
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
    mean = &propMean[ii*nDims];
    sd = &propSd[ii*nDims2];
    Z = &propZ[ii*nDims];
    mvEuler(mean, sd, &currX[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currTheta, &sde[ii]);
    // partial observations
    if(nObsComp[ii] == 0) {
      for(jj = 0; jj < nDims; jj++) {
	Z[jj] = norm_rand();
      }
    }
    else {
      zmvn(Z, &currX[ii*nDims], mean, sd, nDims, nObsComp[ii]);
      for(jj=nObsComp[ii]; jj<nDims; jj++) {
	Z[jj] = norm_rand();
      }
    }
    // proposals
    xmvn(&propX[ii*nDims], Z, mean, sd, nDims);
    //if(isNaN(&propX[ii*nDims], nDims)) {
    //  Rprintf("found NaN in missGibbsUpdate, ii = %i.\n", ii);
    //  foundNaN = 1;
    //  goto stop;
    //}
    // acceptance is 100% as long as the proposal is valid
    if(sde[ii].isValidData(&propX[ii*nDims], currTheta)) {
      for(jj = 0; jj < nDims; jj++) {
	currX[ii*nDims + jj] = propX[ii*nDims + jj];
      }
      gibbsAccept[ii]++;
    }
  }
  // special treatment for first datapoint
  if(missInd[0] == 0) {
    ii = 0;
    mean = &propMean[ii*nDims];
    sd = &propSd[ii*nDims2];
    Z = &propZ[ii*nDims];
    // initialize
    for(jj = 0; jj < nDims; jj++) {
      propX[jj] = currX[jj];
    }
    // random walk metropolis
    for(jj = 0; jj < nMiss0; jj++) {
      // proposal
      propX[nObsComp[0]+jj] = currX[nObsComp[0]+jj] + jumpSd[nParams+jj] * norm_rand();
      if(sde[ii].isValidData(&propX[ii*nDims], currTheta)) {
	// acceptance rate.
	// target 1
	propAccept[ii] = prior->logPrior(currTheta, propX);
	propAccept[ii] -= prior->logPrior(currTheta, currX);
	// target 2
	mvEuler(mean, sd, &propX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] += lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	mvEuler(mean, sd, &currX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[ii]);
	propAccept[ii] -= lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
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
  //stop:
  //return foundNaN;
}


#endif
