#ifndef mvnUtils_h
#define mvnUtils_h 1

#include <Rcpp.h>
using namespace Rcpp;

// x = sd * z + mean.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
inline void xmvn(double *x, double *z, double *mean, double *cholSd, int n) {
  int ii, jj, colI;
  for(ii = 0; ii < n; ii++) {
    colI = n*ii;
    x[ii] = 0;
    for(jj = 0; jj <= ii; jj++) x[ii] += cholSd[colI + jj] * z[jj];
    x[ii] += mean[ii];
  }
  return;
}

// z = sd^{-1} * (x - mean).  only calculates first nMax values of z.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
inline void zmvn(double *z, double *x, double *mean, double *cholSd, int n, int nMax) {
  int ii, jj, colI;
  double tmpSum;
  for(ii = 0; ii < nMax; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii = 0; ii < nMax; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  return;
}

// log-normal density evaluation.  z[] is required as temporary storage of residuals.
// i.e., z = sd^{-1} * (x - mean)
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
// TODO: include pi factor
inline double lmvn(double *x, double *z, double *mean, double *cholSd, int n) {
  double tmpSum2 = 0.0;
  double tmpSum3 = 0.0;
  double resi, tmpSum, val;
  int ii, colI, jj;
  // forward substitution
  colI = 0;
  for(ii = 0; ii < n; ii++) {
    resi = x[ii] - mean[ii];
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    val = (resi - tmpSum) / cholSd[colI + ii];
    tmpSum3 += log(cholSd[colI + ii]);
    z[ii] = val;
    tmpSum2 += (val * val);
    colI += n;
  }
  tmpSum2 *= 0.5;
  tmpSum2 += tmpSum3;
  return(-tmpSum2);
}

/*
// iid normal draws
void rnormiid(double *z, int n) {
  for(int ii=0; ii<n; ii++) {
    z[ii] = norm_rand();
  }
  return;
}
*/

#endif
