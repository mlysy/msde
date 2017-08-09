// Gaussian Copula prior

#ifndef sdePrior_h
#define sdePrior_h 1

#include <Rcpp.h>
#include <mvnUtils.h>

class sdePrior {
  private:
  int nRV, nParamRV, nDataRV; // number of active variables of each type
  int *paramId, *dataId; // index vectors
  int *nBreaks;
  double *range, *dx, *pdf, *logPdf, *cdf;
  double *mean, *sd, *RhoCholSd;
  double *qNorm, *tmpX, *mean0;
 public:
  double logPrior(double *theta, double *x); // log-prior function
  sdePrior(double **phi, int nArgs, int *nEachArg); // constructor
  ~sdePrior(); // destructor
};

// argument order is:
// 0. nBreaks, 1. range, 2. dx, 3. pdf, 4. logPdf, 5. cdf,
// 6. mean, 7. sd, 8. RhoCholSd
// 9. paramId, 10. dataId
inline sdePrior::sdePrior(double **phi, int nArgs, int *nEachArg) {
  int ii,jj;
  nRV = nEachArg[0];
  nParamRV = nEachArg[9];
  nDataRV = nEachArg[10];
  if(nRV > 0) {
    // allocate memory
    nBreaks = new int[nRV];
    range = new double[2*nRV];
    dx = new double[nRV];
    pdf = new double[nEachArg[3]];
    logPdf = new double[nEachArg[4]];
    cdf = new double[nEachArg[5]];
    mean = new double[nRV];
    sd = new double[nRV];
    RhoCholSd = new double[nRV*nRV];
    qNorm = new double[nRV];
    tmpX = new double[nRV];
    mean0 = new double[nRV];
    // assign values
    for(ii=0; ii<nRV; ii++) {
      nBreaks[ii] = (int) phi[0][ii];
      range[2*ii] = phi[1][2*ii];
      range[2*ii+1] = phi[1][2*ii+1];
      dx[ii] = phi[2][ii];
      mean[ii] = phi[6][ii];
      sd[ii] = phi[7][ii];
      mean0[ii] = 0.0;
    }
    jj = 3;
    for(ii=0; ii<nEachArg[jj]; ii++) {
      pdf[ii] = phi[jj][ii];
    }
    jj = 4;
    for(ii=0; ii<nEachArg[jj]; ii++) {
      logPdf[ii] = phi[jj][ii];
    }
    jj = 5;
    for(ii=0; ii<nEachArg[jj]; ii++) {
      cdf[ii] = phi[jj][ii];
    }
    jj = 8;
    for(ii=0; ii<nEachArg[jj]; ii++) {
      RhoCholSd[ii] = phi[jj][ii];
    }
    // index arrays for each type of variable
    jj = 9;
    if(nParamRV > 0) {
      paramId = new int[nParamRV];
      for(ii=0; ii<nParamRV; ii++) {
	paramId[ii] = (int) phi[jj][ii] - 1;
	//Rprintf("paramId[%i] = %i\n", ii, paramId[ii]);
      }
    }
    jj = 10;
    if(nDataRV > 0) {
      dataId = new int[nDataRV];
      for(ii=0; ii<nDataRV; ii++) {
	dataId[ii] = (int) phi[jj][ii] - 1;
	//Rprintf("dataId[%i] = %i\n", ii, dataId[ii]);
      }
    }
  }
}

// destructor
inline sdePrior::~sdePrior() {
  if(nRV > 0) {
    delete [] nBreaks;
    delete [] range;
    delete [] dx;
    delete [] pdf;
    delete [] logPdf;
    delete [] cdf;
    delete [] mean;
    delete [] sd;
    delete [] RhoCholSd;
    delete [] qNorm;
    delete [] tmpX;
    delete [] mean0;
    if(nParamRV > 0) {
      delete [] paramId;
    }
    if(nDataRV > 0) {
      delete [] dataId;
    }
  }
}


// log-prior function
inline double sdePrior::logPrior(double *theta, double *x) {
  if(nRV == 0) return(0.0);
  double lp = 0.0;
  int n = nRV;
  int ii, jj, colI, densElt, start;
  double tmpSum;
  // store
  for(ii=0; ii<nParamRV; ii++) {
    tmpX[ii] = theta[paramId[ii]];
  }
  for(ii=0; ii<nDataRV; ii++) {
    tmpX[nParamRV+ii] = x[dataId[ii]];
  }
  // normal quantiles and marginal components
  start = 0;
  //Rprintf("marginal component:\n");
  for(ii = 0; ii < n; ii ++) {
    densElt = (int) floor((tmpX[ii]-range[2*ii])/dx[ii]);
    if((densElt >= 0) & (densElt < nBreaks[ii])) {
      lp += logPdf[densElt + start];
      qNorm[ii] = tmpX[ii] - (range[2*ii] + densElt * dx[ii]);
      qNorm[ii] *= pdf[densElt + start];
      qNorm[ii] += cdf[densElt + start + ii];
      qNorm[ii] = Rf_qnorm5(qNorm[ii], 0.0, 1.0, 1, 0);
      /* qNorm[ii] = Rf_qnorm5(cdf[densElt + start + ii] +  qNorm[ii], 0.0, 1.0, 1, 0); */
    }
    else {
      lp += Rf_dnorm4(tmpX[ii], mean[ii], sd[ii], 1);
      qNorm[ii] = (tmpX[ii] - mean[ii])/sd[ii];
    }
    start += nBreaks[ii];
    //Rprintf("lp[%i] = %f\n", ii, lp);
    //Rprintf("qNorm[%i] = %f\n", ii, qNorm[ii]);
  }
  // copula components
  // iid standard normal densities
  //Rprintf("standard normals:\n");
  tmpSum = 0.0;
  for(ii = 0; ii < n; ii++) {
    tmpSum += qNorm[ii] * qNorm[ii];
    //Rprintf("tmpSum[%i] = %f\n", ii, tmpSum);
  }
  lp += 0.5 * tmpSum;
  // multivariate normal density
  //Rprintf("mvn density:\n");
  lp += lmvn_chol(qNorm, tmpX, mean0, RhoCholSd, n);
  //Rprintf("lp = %f\n", lp);
  /* for(ii = 0; ii < n; ii++) { */
  /*   colI = n*ii; */
  /*   tmpSum = 0.0; */
  /*   for(jj = 0; jj < ii; jj++) { */
  /*     tmpSum += RhoCholSd[colI + jj] * qNorm[jj]; */
  /*   } */
  /*   qNorm[ii] = (qNorm[ii] - tmpSum)/RhoCholSd[colI + ii]; */
  /* } */
  /* tmpSum = 0.0; */
  /* for(ii = 0; ii < n; ii++) { */
  /*   tmpSum += qNorm[ii] * qNorm[ii]; */
  /* } */
  /* tmpSum *= 0.5; */
  /* for(ii = 0; ii < n; ii++) { */
  /*   tmpSum += log(RhoCholSd[n*ii + ii]); */
  /* } */
  /* lp -= tmpSum; */
  return(lp);
}

#endif
