#ifndef sdePF_h
#define sdePF_h 1

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppSMC")]]
//[[Rcpp::depends("msde")]]
#include <RcppArmadillo.h>
#include <smctc.h>
#include <rngUtils.h>
#include <sdeUtils.h>
#include <sdeParticle.h>
#include <sdeAlgPtr.h>
#include <sdeAdapt.h>


/// Initialise smc::sampler
template <class sMod>
void fInitialise(sdeParticle<sMod>& value, double& logweight,
     sdeAlgPtr<sMod> & algParams) {
  //Rprintf("made it to fInitialise.\n");
  // set particle
  logweight = algParams.init(value.yObs);
  algParams.increase_counter();
  return;
}

/// Move smc::sampler
template <class sMod>
void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
     sdeAlgPtr<sMod> & algParams) {
  int iCore = 0;
  //Rprintf("counter = %i\n", algParams.get_counter());
  logweight += algParams.update(value.yObs, value.yObs, lTime,
             algParams.get_counter(), iCore);
  algParams.increase_counter();
  return;
}


/// Particle filter for SDE models
template <class sMod>
class sdePF {
 private:
  /// particle object required for extracting particles from smc::sampler
  sdeParticle<sMod> tParticle;
  /// Adapt object required to access different memory per particle
  smc::adaptMethods<sdeParticle<sMod>, sdeAlgPtr<sMod> > *Adapt;
  /// pointer to Sampler object
  smc::sampler<sdeParticle<sMod>, sdeAlgPtr<sMod> > *Sampler;
  /// pointer to Moveset object
  smc::moveset<sdeParticle<sMod>, sdeAlgPtr<sMod> > *Moveset;
  /// pointer to AlgParams object
  sdeAlgParams<sMod> *Params;
  /// common initialization of Sampler/Moveset
  void ctor_init(int npart, ResampleType::Enum resample, double dThreshold);
 public:
  /// constructor: run-time random draws
  sdePF(int npart, int ncomp, double *dt, int *yIndex,
	double *yInit, double *thetaInit,
	ResampleType::Enum resample, double dThreshold);
  /// constructor: precomputed random draws
  sdePF(int npart, int ncomp, double *dt, int *yIndex,
	double *yInit, double *thetaInit, double *zInit,
	ResampleType::Enum resample, double dThreshold);
  /// destructor
  ~sdePF();
  /// setter for theta
  void set_theta(double *thetaIn);
  /// setter for first observation
  void set_yInit(double *yIn);
  /// getter for all particles and corresponding loglikelihooods at current step
  void eval(double *yOut, double *lwgt);
  /// wrapper for Sampler->Initialise
  void Initialise() {
    Sampler->Initialise();
    return;
  }
  /// wrapper for Sampler->Iterate 
  void Iterate() {
    Sampler->Iterate();
    return;
  }
};

// allocate memory for Adapt, Sampler, Moveset
// and pass these on to Sampler via Sampler->SetXYZ
template <class sMod>
inline void sdePF<sMod>::ctor_init(int npart,
				   ResampleType::Enum resample,
				   double dThreshold) {
  /// allocate memory
  Adapt = new sdeAdapt<sMod>;
  Sampler = new smc::sampler<sdeParticle<sMod>, sdeAlgPtr<sMod> >((long)npart, HistoryType::NONE);
  Moveset = new smc::moveset<sdeParticle<sMod>, sdeAlgPtr<sMod> >(fInitialise<sMod>, fMove<sMod>, NULL);
  /// pass to Sampler
  Sampler->SetResampleParams(resample, dThreshold);
  Sampler->SetMoveSet(*Moveset);
  Sampler->SetAdaptMethods(Adapt);
  return;
}

/// Call ctor_init, allocate memory for Params, then pass on to Sampler
/// \param npart Number of particles
/// \param ncomp Number of observations
/// \param dt Interobservation times
/// \param yIndex Number of observed components per timepoint
/// \param yInit Matrix of observed data.  Entries corresponding to unobserved components are ignored.
/// \param thetaInit Initial parameter value.
/// \param resample Type of resampling.
/// \param dThreshold ESS threshold for resampling.
template <class sMod>
inline sdePF<sMod>::sdePF(int npart, int ncomp, double *dt, int *yIndex,
			  double *yInit, double *thetaInit,
			  ResampleType::Enum resample,
			  double dThreshold) {
  ctor_init(npart, resample, dThreshold);
  /// pointer to AlgParams
  Params = new sdeAlgParams<sMod>(npart, ncomp, dt, yIndex, yInit, thetaInit);
  Sampler->SetAlgParam(sdeAlgPtr<sMod>(Params));
}

/// Call ctor_init, allocate memory for Params, then pass on to Sampler
/// \param npart Number of particles
/// \param ncomp Number of observations
/// \param dt Interobservation times
/// \param yIndex Number of observed components per timepoint
/// \param yInit Matrix of observed data.  Entries corresponding to unobserved components are ignored.
/// \param thetaInit Initial parameter value.
/// \param zInit Array of precomputed normal draws.  Ordering is nDims x nPart x nComp.
/// \param resample Type of resampling.
/// \param dThreshold ESS threshold for resampling.
template <class sMod>
inline sdePF<sMod>::sdePF(int npart, int ncomp, double *dt, int *yIndex,
			  double *yInit, double *thetaInit, double *zInit,
			  ResampleType::Enum resample,
			  double dThreshold) {
  ctor_init(npart, resample, dThreshold);
  /// pointer to AlgParams
  Params = new sdeAlgParams<sMod>(npart, ncomp, dt, yIndex,
				  yInit, thetaInit, zInit);
  Sampler->SetAlgParam(sdeAlgPtr<sMod>(Params));
}

/// destructor
template <class sMod>
inline sdePF<sMod>::~sdePF() {
  delete Adapt;
  delete Sampler;
  delete Moveset;
  delete Params;
}

/// Logweights in Sampler are in fact normalized and divided by nPart,
/// so must undo this to get the unnormalized weights.
/// \param yOut Array of particle values at current step (lTime).  Ordering is nDims x nPart.
/// \param lwgt Array of corresponding PF likelihoods (unnormalized).
template <class sMod>
inline void sdePF<sMod>::eval(double *yOut, double *lwgt) {
  double logNC = Sampler->GetLogNCPath() + log(Sampler->GetNumber());
  for(long iPart=0; iPart<Sampler->GetNumber(); iPart++) {
    tParticle = Sampler->GetParticleValueN(iPart);
    tParticle.get_yObs(&yOut[iPart*sMod::nDims]);
    lwgt[iPart] = Sampler->GetParticleLogWeightN(iPart) + logNC;
  }
  return;
}

/// \param thetaIn Parameter vector.
template <class sMod>
inline void sdePF<sMod>::set_theta(double *thetaIn) {
  for(int ii=0; ii<sMod::nParams; ii++) {
    Params->theta[ii] = thetaIn[ii];
  }
  return;
}

template <class sMod>
inline void sdePF<sMod>::set_yInit(double *yIn) {
  for(int ii=0; ii<sMod::nDims; ii++) {
    Params->yObs[ii] = yIn[ii];
  }
  return;
}

#endif
