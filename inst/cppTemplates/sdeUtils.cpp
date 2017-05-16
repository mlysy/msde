#include <Rcpp.h>
using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

//[[Rcpp::depends("msdeHeaders")]]
#include <sdeMCMC.h>
#include <mcmcUtils.h>

// some utilities to be exported, namely for checking the c++ code
//
// the function are: sde's drift, diffusion, log-likelihood,
// and eventually custom prior.

//[[Rcpp::export("ndims")]]
int sde_getNDims() {
  return sdeModel::nDims;
}

//[[Rcpp::export("nparams")]]
int sde_getNParams() {
  return sdeModel::nParams;
}

// Best to parallelize these from within R: easier and prob faster since
// memory can be allocated in parallel there
//[[Rcpp::export("sde.model$drift")]]
NumericVector sde_Drift(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sdeModel *sde = new sdeModel[1];
  for(int ii = 0; ii < nReps; ii++) {
    sde[0].sdeDr(&dr[ii*nDims], &x[ii*nDims], &theta[ii*nParams]);
  }
  delete [] sde;
  return drOut;
}

//[[Rcpp::export("sde.model$diff")]]
NumericVector sde_Diff(NumericVector xIn, NumericVector thetaIn, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector dfOut(nReps*nDims*nDims);
  double *df = REAL(dfOut);
  sdeModel *sde = new sdeModel[0];
  for(int ii = 0; ii < nReps; ii++) {
    sde[0].sdeDf(&df[ii*nDims*nDims], &x[ii*nDims], &theta[ii*nParams]);
  }
  delete [] sde;
  return dfOut;
}

// SDE log-likelihood evaluation.
//[[Rcpp::export("sde.model$loglik")]]
NumericVector sde_LogLik(NumericVector xIn, NumericVector dTIn,
			NumericVector thetaIn,
			int nComp, int nReps) {
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  NumericVector llOut(nReps);
  double *ll = REAL(llOut);
  sdeLogLik sdeLL(nComp, REAL(dTIn));
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*nParams], &x[ii*nDims*nComp]);
  }
  return llOut;
}

//[[Rcpp::export("sde.model$logprior")]]
NumericVector sde_Prior(NumericVector thetaIn, NumericVector xIn,
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
    lp[ii] = prior.logPrior(&theta[ii*nParams], &x[ii*nDims]);
  }
  delete [] phi;
  delete [] nEachArg;
  return lpOut;
}


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
  double *Z = new double[nDims]; // random draw

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
	} while(!sde->isValidData(tmpX, &params[ii*nParams]) &&
		bad++ < MAXBAD);
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

//[[Rcpp::export("sde.model$post")]]
List sdeEulerMCMC(NumericVector initParams, NumericVector initData,
		  NumericVector dT, IntegerVector nDimsPerObs,
		  LogicalVector fixedParams,
		  int nSamples, int burn,
		  int nParamsOut, int nDataOut,
		  IntegerVector dataOutRow, IntegerVector dataOutCol,
		  double updateParams, double updateData,
		  List priorArgs, List tunePar,
		  int updateLogLik, int nLogLikOut,
		  int updateLastMiss, int nLastMissOut) {
  RNGScope scope;
  int ii, jj, kk;

  // problem dimensions
  int nDims = sdeModel::nDims;
  int nParams = sdeModel::nParams;
  int nComp = initData.length()/nDims;
  int nCompOut = dataOutCol.length();
  int nMiss0 = nDims-nDimsPerObs[0]; // unobserved states in first observation
  int nMissN = nDims-nDimsPerObs[nComp-1]; // unobserved states in last observation

  // output variables
  NumericVector paramsOut(nParamsOut);
  NumericVector dataOut(nDataOut);
  IntegerVector paramAcceptOut(nParams + nMiss0);
  IntegerVector gibbsAcceptOut(nComp);
  NumericVector logLikOut(nLogLikOut);
  NumericVector lastMissOut(nLastMissOut);
  NumericVector lastIter(nParams + nComp*nDims);
  NumericVector mwgSdOut(nParams + nDims);
  // pointers to acceptance rate counters for internal use
  int *paramAccept = INTEGER(paramAcceptOut);
  int *gibbsAccept = INTEGER(gibbsAcceptOut);
  double *mwgSd = REAL(mwgSdOut);

  // MCMC tuning parameters
  for(ii=0; ii<nParams+nDims; ii++) {
    mwgSd[ii] = REAL(tunePar["sd"])[ii];
  }
  mwgAdapt tuneMCMC(REAL(tunePar["max"]), REAL(tunePar["rate"]),
		    LOGICAL(tunePar["adapt"]), nParams+nDims);
  // double *jumpSd = REAL(rwJumpSd); // random walk jump sizes
  // double *mwgSd = new double[nParams + nDims];
  // double *adaptMax = new double[nParams + nDims];
  // double *adaptRate = new double[nParams + nDims];
  // bool *doAdapt = new bool[nParams + nDims];
  // for(ii=0; ii<nParams+nDims; ii++) {
  //   mwgSd[ii] = REAL(tunePar["sd"])[ii];
  //   adaptMax[ii] = REAL(tunePar["max"])[ii];
  //   adaptRate[ii] = REAL(tunePar["rate"])[ii];
  //   doAdapt[ii] = (LOGICAL(tunePar["adapt"])[ii] != 0); // explicit logical cast
  // }

  // prior specification
  // hyper parameters: actual prior gets constructed inside MCMC object
  int nArgs = priorArgs.length();
  double **phi = new double*[nArgs];
  int *nEachArg = new int[nArgs];
  for(ii=0; ii<nArgs; ii++) {
    if(Rf_isNull(priorArgs[ii])) {
      nEachArg[ii] = 0;
    } else {
      nEachArg[ii] = as<NumericVector>(priorArgs[ii]).length();
      phi[ii] = REAL(priorArgs[ii]);
    }
  }

  // initialize MCMC
  // prior gets constructed inside of object -- is this really beneficial?
  sdeMCMC mcmc(nComp, REAL(dT), REAL(initData), REAL(initParams),
	       INTEGER(nDimsPerObs), LOGICAL(fixedParams),
	       phi, nArgs, nEachArg);

  // main MCMC loop
  jj = 0;
  for(int smp = -burn; smp < nSamples; smp++) {
    // user interrupt
    if(smp % (int) 1e3) {
      Rcpp::checkUserInterrupt();
    }
    // missing data update
    if(updateComponent(updateData, smp)) {
      mcmc.missGibbsUpdate(mwgSd, gibbsAccept, paramAccept);
    }
    // parameter update
    if(updateComponent(updateParams, smp)) {
      mcmc.paramVanillaUpdate(mwgSd, paramAccept);
    }
    // adaptive MCMC
    tuneMCMC.adapt(mwgSd, paramAccept, burn+smp+1);
    // log-likelihood
    // TODO: keep track of this interally after every MCMC step
    if(smp >= 0) {
      if(updateLogLik) logLikOut[smp] = mcmc.loglik(mcmc.currTheta, mcmc.currX);
    }
    // storage
    if(smp == dataOutRow[jj]) {
      if(updateData > 0.0) {
	for(ii=0; ii<nCompOut; ii++) {
	  for(kk=0; kk<nDims; kk++) {
	    dataOut[jj*nDims*nCompOut+ii*nDims+kk] = mcmc.currX[dataOutCol[ii]*nDims+kk];
	  }
	}
      }
      jj++;
    }
    if((updateParams > 0.0) && (smp >= 0)) {
      for(ii=0; ii<nParams; ii++) paramsOut[smp*nParams + ii] = mcmc.currTheta[ii];
    }
    if(updateLastMiss && smp >= 0) {
      for(ii=0; ii<nMissN; ii++) {
	lastMissOut[smp*nMissN+ii] = mcmc.currX[(nComp-1)*nDims+nDimsPerObs[nComp-1]+ii];
      }
    }
    // NaN breakpoint
    //if(foundNaN) {
    //  Rprintf("smp = %i\n", smp);
    //  goto stop;
    //}
  }

  // store last iteration, to resume MCMC later if needed
  for(ii=0; ii<nParams; ii++) {
    lastIter[ii] = mcmc.currTheta[ii];
  }
  for(ii=0; ii<nDims*nComp; ii++) {
    lastIter[nParams+ii] = mcmc.currX[ii];
  }

  //stop:
  // delete dynamic variables
  delete [] phi;
  delete [] nEachArg;

  return List::create(_["paramsOut"] = paramsOut,
		      _["dataOut"] = dataOut,
		      _["paramAccept"] = paramAcceptOut,
		      _["gibbsAccept"] = gibbsAcceptOut,
		      _["logLikOut"] = logLikOut,
		      _["lastMissOut"] = lastMissOut,
		      _["lastIter"] = lastIter,
		      _["mwgSd"] = mwgSdOut);
}
