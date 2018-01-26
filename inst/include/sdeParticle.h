/*
// particle filter for SDEs.

New design: create sampler inside class containing memory.

class sdeParticle: as before, the particle object.

class sdeFilter:
  - all memory is now private!
  - private SMC members: fInitialize, fMove, save_state
  - public member: Sampler, and set it up properly.  so no more algparams
  - as all private members of sdeFilter are global variables from the perspective of fInitialize, fMove, save_state, etc.

*/

#ifndef sdeParticle_h
#define sdeParticle_h

// -----------------------------------------------------------

// the particle itself
template <class sMod>
class sdeParticle {
 private:
  static const int nDims = sMod::nDims;
 public:
  double *yObs;
  sdeParticle() {
    yObs = new double[nDims];
  }
  int get_nDims() {
    return nDims;
  }
  void get_yObs(double *yOut) {
    for(int ii=0; ii<nDims; ii++) {
      yOut[ii] = yObs[ii];
    }
    return;
  }
  void set_yObs(double *yIn) {
    for(int ii=0; ii<nDims; ii++) {
      yObs[ii] = yIn[ii];
    }
    return;
  }
  ~sdeParticle() {
    //Rprintf("sdeParticle destructor called.\n");
    delete [] yObs;
  }
  sdeParticle & operator=(const sdeParticle & other) {
    if(this != &other) {
      set_yObs(other.yObs);
    }
    return *this;
  }
  // deep copy required to extract particle from Sampler.
  sdeParticle(const sdeParticle & from) {
    //Rprintf("sdeParticle copy constructor called.\n");
    yObs = new double[nDims];
    set_yObs(from.yObs);
  }
};

#endif
