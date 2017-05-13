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

// metropolis-within-gibbs adaptive sampler

// inputs:
// max_val: maximum value
// dec_rate: adaption decay rate (default is .5)
// n: current step of algorithm
// acc: current acceptance rate

// output:
// delta(n)

/*
double mwg_adapt(double max_val, double dec_rate,
		 int n, double acc, );

class Tune_aMWG {
  double max_val;
  double dec_rate;
  const double targ_acc = .44;
 public:
  mwg_adapt(double mv, double dr) {
    max_val = mv;
    dec_rate = dr;
  }
  mwg_adapt() {
    max_val = .01;
    dec_rate = .5;
  }
  double adapt(int n, double acc, double sigma) {
    double delta = pow(n, dec_rate);
    double ls_new = log(sigma);
    if(max_val < delta) delta = max_val;
    if(acc < targ_acc) {
      ls_new -= delta;
    }
    else {
      ls_new += delta;
    }
    return(exp(ls_new));
  }
}
*/

// starts with log(rwJumpSd) = 0
// at each step moves up or down by delta(n) according to whether
// acceptance rate is +/- target value .44, where

// delta(n) = min(max_val, n^dec_rate)



#endif
