/*
// --- SMC for multivariate SDEs -----------------------------------------------

Limitations:
- undefined behavior for invalid data (isValidData = FALSE).  This is because
- no filtering, only smoothing.
- no parallel processing.

Implementation:
* sdeParticle: 
  - comprises of the complete data at that point.
    storage for calculations, etc. are elsewhere.

* sdeFilter: 
  - everything we need to calculate particle weights, i.e.,
  - data matrix, dT/sqrtDT, and incomplete observation index.
  - propMean/Sd, sde objects, etc.

* random draws:
  - normal draws should be sitting around somewhere as well, with the option to 
    redraw these.  
  - this should be done with the adaptor, but perhaps a small modification to 
    code base would have particle index in moveset construction...

* updates:
  - get current complete data point.
  - draw from conditional distribution of previous data point and this one for
    unobserved components.
  - compute loglikelihood.
  - update particle.

*/

#ifndef sdeSMC_h
#define sdeSMC_h 1

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
//#include <vector>
#include <RcppArmadillo.h>
#include <smctc.h>
#include <rngUtils.h>
#include <sdeUtils.h>

// observed data, i.e., data matrix, dT/sqrtDT, nvarObs
template <class sMod>
class sdeFilter {
 private:
  bool drawZ; // whether or not normal draws are pre-computed
  static const int nCores = 1; // for parallel processing.  change this later.
  int nDims2;
  int nPart; // number of particles
  int iPart; // which particle are we on
  // internal constructor (shared by overloaded constructors)
  void initialize(int np, int nc, double *dt, int *yIndex,
		  double *yInit, double *thetaInit);
  void clear(); // internal destructor (for delete and deep-copy overwrite)
 public:
  int nDims, nParams, nComp;
  int *nObsComp;
  double *yObs, *yProp; // observed and proposal data arrays
  double *theta; // parameter value
  double *dT, *sqrtDT;
  double *propMean, *propSd, *propZ; // for storing normal calculations
  sMod *sde;
  sdeFilter(int np, int nc, double *dt, int *yIndex,
	    double *yInit, double *thetaInit);
  sdeFilter(int np, int nc, double *dt, int *yIndex,
	    double *yInit, double *thetaInit, double *zInit);
  sdeFilter(); // default constructor
  sdeFilter & operator=(const sdeFilter & other); // deep copy
  double init(double *yNew);
  double update(double *yNew, double *yOld, long lTime, int iPart, int iCore);
  int get_counter() const {return iPart;}
  void increase_counter() {iPart++; return;}
  void reset_counter() {iPart = 0;}
  ~sdeFilter();
};

// assumed that yOld contains the right observed data where it ought to.
template <class sMod>
inline double sdeFilter<sMod>::update(double *yNew, double *yOld,
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

template <class sMod>
inline double sdeFilter<sMod>::init(double *yNew) {
  for(int ii=0; ii<nDims; ii++) {
    yNew[ii] = yObs[ii];
  }
  return 0;
}

template <class sMod>
inline void sdeFilter<sMod>::initialize(int np, int nc, double *dt,
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

// this version draws its own normals at every step.
template <class sMod>
inline sdeFilter<sMod>::sdeFilter(int np, int nc, double *dt, int *yIndex,
				  double *yInit, double *thetaInit) {
  initialize(np, nc, dt, yIndex, yInit, thetaInit);
  drawZ = true;
  propZ = new double[nCores*nDims];
}

// this version takes a pointer to normal draws and does not draw them ever.
// since some parts of the normal matrix does get updated at every step, normal
// draws are copied in at construction.
// this is mainly for debugging but could be useful elsewhere.
template <class sMod>
inline sdeFilter<sMod>::sdeFilter(int np, int nc, double *dt, int *yIndex,
				  double *yInit, double *thetaInit,
				  double *zInit) {
  initialize(np, nc, dt, yIndex, yInit, thetaInit);
  drawZ = false;
  propZ = new double[nPart*(nComp-1)*nDims];
  for(int ii=0; ii<nPart*(nComp-1)*nDims; ii++) {
    propZ[ii] = zInit[ii];
  }
}

template <class sMod>
inline void sdeFilter<sMod>::clear() {
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
template <class sMod>
inline sdeFilter<sMod>::sdeFilter() {
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
// required to pass in AlgParams to Sampler.
// (which first creates AlgParams with "default ctor", then deep-copies
// to replace it when setAlgParams is called)
template <class sMod>
inline sdeFilter<sMod> &
sdeFilter<sMod>::operator=(const sdeFilter<sMod> & other) {
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

template <class sMod>
inline sdeFilter<sMod>::~sdeFilter() {
  //Rprintf("sdeFilter destructor called.\n");
  clear();
}

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

// adaptor class
// want to keep track of particle in fMove/fInitialise
// also want to input new theta at the beginning of every pf calculation.
template <class sMod>
class sdeAdapt: public smc::adaptMethods<sdeParticle<sMod>, sdeFilter<sMod> >
{
 private:
  bool updateTheta;
  double *theta;
  int nParams = sMod::nParams;
 public:
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
  void updateForMove(sdeFilter<sMod> & pf_calcs, const smc::population<sdeParticle<sMod> > & pop) {
    // set particle count to zero; each fMove increases it
    // cheap way of keeping track of particle in fMove/fInitialize
    pf_calcs.reset_counter();
    return;
  }
  ~sdeAdapt() {};
};

#endif
