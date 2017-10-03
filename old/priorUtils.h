#ifndef priorUtils_h
#define priorUtils_h 1

#include<Rcpp.h>
using namespace Rcpp;

// parse prior arguments from R list
class priorArgs {
 public:
  int nArgs;
  double **phi;
  int *nEachArg;
  priorArgs(List phiIn);
  ~priorArgs() {
    delete [] nEachArg;
    delete [] phi;
  }
}

inline priorArgs::priorArgs(List phiIn) {
  nArgs = phiIn.length(); // at least 1, since hyper at least list(NULL)
  phi = new double*[nArgs];
  nEachArg = new int[nArgs];
  for(int ii=0; ii<nArgs; ii++) {
    if(Rf_isNull(phiIn[ii])) {
      nEachArg[ii] = 0;
    } else {
      nEachArg[ii] = as<NumericVector>(phiIn[ii]).length();
      phi[ii] = REAL(phiIn[ii]);
    }
  }
}

#endif
