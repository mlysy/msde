#ifndef sdeAdapt_h
#define sdeAdapt_h 1

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
//#include <vector>
#include <RcppArmadillo.h>
#include <smctc.h>
#include <rngUtils.h>
#include <sdeUtils.h>

// adaptor class
// want to keep track of particle in fMove/fInitialise
// also want to input new theta at the beginning of every pf calculation.
template <class sMod>
class sdeAdapt: public smc::adaptMethods<sdeParticle<sMod>, sdeAlgParams<sMod> >
{
 private:
  bool updateTheta;
  double *theta;
  static const int nParams = sMod::nParams;
 public:
  sdeAdapt() {}
  sdeAdapt(double *theta_in) {
    theta = new double[nParams];
    set_theta(theta_in);
    updateTheta = false;
  }
  void set_theta(double *theta_in) {
    for(int ii=0; ii<nParams; ii++) {
      theta[ii] = theta_in[ii];
    }
  }
  void updateForMove(sdeAlgParams<sMod> & pf_calcs, const smc::population<sdeParticle<sMod> > & pop) {
    // set particle count to zero; each fMove increases it
    // cheap way of keeping track of particle in fMove/fInitialize
    pf_calcs.reset_counter();
    return;
  }
  ~sdeAdapt() {};
};

#endif
