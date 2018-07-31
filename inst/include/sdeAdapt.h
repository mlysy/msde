#ifndef sdeAdapt_h
#define sdeAdapt_h 1

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
#include <RcppArmadillo.h>
#include <smctc.h>
#include <rngUtils.h>
#include <sdeUtils.h>

// adaptor class
// the point of this call is to keep track of the
// particle number in fMove/fInitialise
// want to keep track of particle in fMove/fInitialise
template <class sMod>
class sdeAdapt: public smc::adaptMethods<sdeParticle<sMod>, sdeAlgParams<sMod> >
{
 public:
  // default constructor/destructor
  sdeAdapt() {}
  ~sdeAdapt() {};
  // set particle count to zero; each fMove increases it
  // cheap way of keeping track of particle in fMove/fInitialize
  void updateForMove(sdeAlgParams<sMod> & algParams, const smc::population<sdeParticle<sMod> > & pop) {
    algParams.reset_counter();
    return;
  }
};

#endif
