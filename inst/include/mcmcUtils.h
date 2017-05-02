#ifndef mcmcUtils_h
#define mcmcUtils_h 1

#include <Rcpp.h>
using namespace Rcpp;

// x is either:
// T/F: in which case it's a yes or no
// integer, in which case it's every so many iterations (i.e., iter is a multiple of x)
// fraction between 0 an 1, in which case it's randomly updated that fraction of the time
inline bool updateComponent(double x, int iter) {
  bool update = false;
  if(x > 0.0) {
    if(x < 1.0) {
      update = unif_rand() <= x;
    }
    else {
      update = (iter % (int) x) == 0;
    }
  }
  return update;
}


#endif
