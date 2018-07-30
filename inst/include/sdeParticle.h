#ifndef sdeParticle_h
#define sdeParticle_h 1

// particle class to pass to smc::sampler
// basically a thin wrapper to the double array of the complete data SDE vector at the given timepoint.

// the particle itself
template <class sMod>
class sdeParticle {
 private:
  static const int nDims = sMod::nDims; // number of dimensions
 public:
  double *yObs; // value of particle (both obs and missing data)
  // getter/setter for particle
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
  // constructor: just allocates memory
  sdeParticle() {
    yObs = new double[nDims];
  }
  // int get_nDims() {
  //   return nDims;
  // }
  // destructor: deletes dynamic memory
  ~sdeParticle() {
    //Rprintf("sdeParticle destructor called.\n");
    delete [] yObs;
  }
  // assignment overloading
  sdeParticle & operator=(const sdeParticle & other) {
    if(this != &other) {
      set_yObs(other.yObs);
    }
    return *this;
  }
  // deep copy: required to extract particle from Sampler.
  sdeParticle(const sdeParticle & from) {
    //Rprintf("sdeParticle copy constructor called.\n");
    yObs = new double[nDims];
    set_yObs(from.yObs);
  }
};


#endif
