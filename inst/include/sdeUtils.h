#ifndef sdeUtils_h
#define sdeUtils_h 1

// utilities for sde's

#include "LinAlgUtils.h"
#include "mvnUtils.h"
#include "sdeModel.h"

inline void mvEuler(double *mean, double *sd,
		    double *x, double dT, double sqrtDT,
		    double *theta, sdeModel *sde) {
  // mean = x + drift(x,t,theta)*dT
  sde->sdeDr(mean, x, theta);
  v_mult(mean, dT, sdeModel::nDims);
  for(int ii = 0; ii < sdeModel::nDims; ii++) {
    mean[ii] += x[ii];
  }
  // sd = diff(x,t,theta)*sqrt(dT)
  sde->sdeDf(sd, x, theta);
  if(!sdeModel::sdDiff) {
    chol_decomp(sd, sd, sdeModel::nDims);
  }
  U_mult(sd, sqrtDT, sdeModel::nDims);
  return;
}


#endif
