/*
// --- SMC for multivariate SDEs -----------------------------------------------

Limitations:
- no smoothing, only filtering.
- no parallel processing.

Implementation:
* sdePF: 
  - class with public methods:
     - set_theta (and eventually set_x for initial data),
     - eval, which should return latest weight and particle.
  - storage for data, calculations, etc. are in nested class setAlgParams.

* sdeAlgParams: 
  - everything we need to calculate particle weights, i.e.,
  - data matrix, dT/sqrtDT, and incomplete observation index.
  - propMean/Sd, sde objects, etc.
  - option to precompute random draws (much more storage, but useful for debugging)
  - proposal distribution Euler-Maruyama {t = curr-1} -> {t = curr}, conditioned
    on observed components at {t = curr}.
  - invalid proposals automatically given a weight of zero.  undefined behavior when all particles yield invalid proposals.

* sdeAdapt:
  - currently used to manage particle-specific memory.

* sdeParticle:
  - lightweight class to represent a single particle at a given time, 
    as required by smctc.

* fMove, fInitialise:
  - need to be non-member functions.
  - referred to from sdePF constructor.

*/

#ifndef sdeAlgParams_h
#define sdeAlgParams_h 1

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
#include <RcppArmadillo.h>
#include <smctc.h>
#include <rngUtils.h>
#include <sdeUtils.h>

// observed data, i.e., data matrix, dT/sqrtDT, nvarObs
template <class sMod>
class sdeAlgParams {
 private:
  bool drawZ; // whether or not normal draws are pre-computed
  static const int nCores = 1; // for parallel processing.  change this later.
  int nDims2;
  int nPart; // number of particles
  int iPart; // which particle are we on
  // internal constructor (shared by overloaded constructors)
  void initialize(int np, int nc, double *dt, int *yIndex,
		  double *yInit, double *thetaInit);
  // internal destructor (for delete and deep-copy overwrite)
  void clear();
 public:
  int nDims, nParams, nComp;
  int *nObsComp;
  double *yObs, *yProp; // observed and proposal data arrays
  double *theta; // parameter value
  double *dT, *sqrtDT;
  double *propMean, *propSd, *propZ; // for storing normal calculations
  sMod *sde;
  // constructors: with and without precomputed random draws
  sdeAlgParams(int np, int nc, double *dt, int *yIndex,
	    double *yInit, double *thetaInit);
  sdeAlgParams(int np, int nc, double *dt, int *yIndex,
	    double *yInit, double *thetaInit, double *zInit);
  sdeAlgParams(); // default constructor: needed to pass algParams to smc::sampler
  ~sdeAlgParams(); // destructor
  // deep copy: required to pass in to smc::sampler
  sdeAlgParams & operator=(const sdeAlgParams & other);
  // internal version of fInitialise.
  // that is, initializes particle.yNew with first "row" of yObs
  // also sets logweight = 0.
  double init(double *yNew);
  // internal version of fMove.
  // that is, updates particleOld.yOld -> particleNew.yNew
  // also updates logweight.
  double update(double *yNew, double *yOld, long lTime, int iPart, int iCore);
  // counter required for particles to use different blocks of memory.
  int get_counter() const {return iPart;}
  void increase_counter() {iPart++; return;}
  void reset_counter() {iPart = 0;}
};


// internal version of fInitialise.
// that is, initializes particle.yNew with first "row" of yObs
// also sets logweight = 0.
template <class sMod>
inline double sdeAlgParams<sMod>::init(double *yNew) {
  for(int ii=0; ii<nDims; ii++) {
    yNew[ii] = yObs[ii];
  }
  return 0;
}


// internal version of fMove.
// that is, updates particleOld.yOld -> particleNew.yNew
// also updates logweight.
// assumed that yOld contains the right observed data where it ought to.
template <class sMod>
inline double sdeAlgParams<sMod>::update(double *yNew, double *yOld,
				      long lTime, int iPart, int iCore) {
  int ii;
  double *mean, *sd, *Z, *yTmp, lw;
  int nObs = nObsComp[lTime];
  mean = &propMean[iCore*nDims];
  sd = &propSd[iCore*nDims2];
  yTmp = &yProp[iCore*nDims];
  for(ii=0; ii<nObs; ii++) {
    yTmp[ii] = yObs[lTime*nDims + ii]; // fill in observed data
  }
  // normal draws
  if(drawZ) {
    Z = &propZ[iCore*nDims];
    for(ii=nObs; ii<nDims; ii++) {
      Z[iCore*nDims + ii] = sdeRNG::rnorm();
    }
  } else {
    Z = &propZ[(lTime-1)*nDims*nPart + iPart*nDims];
  }
  // mean and variance calculation
  mvEuler(mean, sd, yOld, dT[lTime-1], sqrtDT[lTime-1],
	  theta, &sde[iCore]);
  // partial observations
  if(nObs > 0 && nObs < nDims) {
    zmvn<sMod>(Z, yTmp, mean, sd, nObs);
  }
  // proposal
  if(nObs < nDims) {
    xmvn<sMod>(yTmp, Z, mean, sd, nObs, nDims);
  }
  // assignment
  for(ii=0; ii<nDims; ii++) {
    yNew[ii] = yTmp[ii];
  }
  // log-weight
  if(sde[iCore].isValidData(yTmp, theta)) {
    lw = lmvn<sMod>(yTmp, Z, mean, sd, nObs);
  } else {
    lw = -1.0/0.0;
  }
  return lw;
}

// internal constructor (shared by overloaded constructors)
// basically just allocates memory and copies data referenced by pointer inputs
template <class sMod>
inline void sdeAlgParams<sMod>::initialize(int np, int nc, double *dt,
					   int *yIndex, double *yInit,
					   double *thetaInit) {
  int ii,jj;
  iPart = 0;
  nComp = nc;
  nPart = np;
  //nCores = ncores;
  nDims = sMod::nDims;
  nDims2 = sMod::diagDiff ? nDims : nDims*nDims;
  nParams = sMod::nParams;
  // create storage space
  yObs = new double[nComp*nDims]; // all data stored here
  dT = new double[nComp-1];
  sqrtDT = new double[nComp-1];
  nObsComp = new int[nComp];
  // everything else just for calcs, so only need nCores of each
  yProp = new double[nCores*nDims];
  // propZ gets allocated elsewhere, depending on whether draws are precomputed.
  sde = new sMod[nCores];
  propMean = new double[nCores*nDims];
  propSd = new double[nCores*nDims2];
  theta = new double[nParams];
  // initialize 
  for(ii=0; ii<nComp; ii++) {
    nObsComp[ii] = yIndex[ii];
    for(jj=0; jj<nDims; jj++) {
      yObs[ii*nDims + jj] = yInit[ii*nDims + jj];
    }
    if(ii < nComp-1) {
      dT[ii] = dt[ii];
      sqrtDT[ii] = sqrt(dT[ii]);
    }
  }
  for(ii=0; ii<nParams; ii++) {
    theta[ii] = thetaInit[ii];
  }
  return;
}

// constructor
// this version draws its own normals at every step.
template <class sMod>
inline sdeAlgParams<sMod>::sdeAlgParams(int np, int nc, double *dt, int *yIndex,
				  double *yInit, double *thetaInit) {
  initialize(np, nc, dt, yIndex, yInit, thetaInit);
  drawZ = true;
  propZ = new double[nCores*nDims];
}

// constructor
// this version takes a pointer to normal draws and does not draw them ever.
// since some parts of the normal matrix does get updated at every step, normal
// draws are copied in at construction.
// this is mainly for debugging but could be useful elsewhere.
template <class sMod>
inline sdeAlgParams<sMod>::sdeAlgParams(int np, int nc, double *dt, int *yIndex,
					double *yInit, double *thetaInit,
					double *zInit) {
  initialize(np, nc, dt, yIndex, yInit, thetaInit);
  drawZ = false;
  propZ = new double[nPart*(nComp-1)*nDims];
  for(int ii=0; ii<nPart*(nComp-1)*nDims; ii++) {
    propZ[ii] = zInit[ii];
  }
}

// internal destructor (for delete and deep-copy overwrite)
template <class sMod>
inline void sdeAlgParams<sMod>::clear() {
  delete [] nObsComp;
  delete [] yObs;
  delete [] yProp;
  delete [] theta;
  delete [] dT;
  delete [] sqrtDT;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] sde;
  return;
}

// default constructor
// needed to pass algParams to smc::sampler
template <class sMod>
inline sdeAlgParams<sMod>::sdeAlgParams() {
  nObsComp = new int[1];
  yObs = new double[1];
  yProp = new double[1];
  theta = new double[1];
  dT = new double[1];
  sqrtDT = new double[1];
  propMean = new double[1];
  propSd = new double[1];
  propZ = new double[1];
  sde = new sMod[1];
}

// deep copy assignment
// required to pass in algParams to smc::sampler.
// (which first creates algParams with "default ctor", then deep-copies
// to replace it when setAlgParams is called)
template <class sMod>
inline sdeAlgParams<sMod> & sdeAlgParams<sMod>::operator=(const sdeAlgParams<sMod> & other) {
  if(this != &other) {
    // deallocate memory
    clear();
    // allocate new memory
    initialize(other.nPart, other.nComp, other.dT, other.nObsComp,
	       other.yObs, other.theta);
    // finalize allocation
    if(other.drawZ) {
      drawZ = true;
      propZ = new double[nCores*nDims];
      for(int ii=0; ii<nCores*nDims; ii++) {
	propZ[ii] = other.propZ[ii];
      }
    } else {
      drawZ = false;
      propZ = new double[nPart*(nComp-1)*nDims];
      for(int ii=0; ii<nPart*(nComp-1)*nDims; ii++) {
	propZ[ii] = other.propZ[ii];
      }
    }
  }
  return *this;
}

// destructor
template <class sMod>
inline sdeAlgParams<sMod>::~sdeAlgParams() {
  //Rprintf("sdeAlgParams destructor called.\n");
  clear();
}

#endif
