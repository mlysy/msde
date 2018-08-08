#ifndef sdeAlgPtr_h
#define sdeAlgPtr_h 1

//[[Rcpp::depends("msde")]]
#include <sdeAlgParams.h>

/// Wrapper class for sdeAlgParams object
///
/// given a pointer (algPtr) of type *sdeAlgParams, provides wrappers to the
/// public methods of sdeAlgParams, i.e., sdeAlgPtr.foo() = algPtr->foo().
/// Thus, algPtr can be changed externally without touching sdeAlgPtr.
template <class sMod>
class sdeAlgPtr {
 private:
  sdeAlgParams<sMod> *pParams;
 public:
  /// constructor
  sdeAlgPtr(sdeAlgParams<sMod> *pparams) {pParams = pparams;}
  /// default constructor
  sdeAlgPtr() {}
  /// copy constructor
  sdeAlgPtr & operator=(const sdeAlgPtr & other) {
    if(this != &other) {
      pParams = other.pParams;
    }
    return *this;
  }
  // wrappers to public methods
  double init(double *yNew) {
    return pParams->init(yNew);
  }
  double update(double *yNew, double *yOld, long lTime,
		int iPart, int iCore) {
    return pParams->update(yNew, yOld, lTime, iPart, iCore);
  }
  int get_counter() const {
    return pParams->get_counter();
  }
  void increase_counter() {
    pParams->increase_counter();
    return;
  }
  void reset_counter() {
    pParams->reset_counter();
    return;
  }
};

#endif
