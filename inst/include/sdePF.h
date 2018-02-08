#ifndef sdePF_h
#define sdePF_h 

// This header file contains the definition of particleEval template and 
// class template sdePF with member functions fInitialise, fMove, save_state)
//
// TODO:
// change the header into a particle filter template class for sde
// member functions: fInitialise, fMove, save_state, particleEval

//[[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <smctc.h>
// using namespace Rcpp;
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
// if we define Rcpp::NumericMatrix as Matrix, an error 
// "reference to 'Matrix' is ambiguous" would be encountered
typedef Rcpp::NumericMatrix NumericMatrix;
typedef Rcpp::List List;
#include "sdeInterface.h"
#include "sdeSMC.h"

// template <class sMod> 
// class sdePF {
//   public:
//       // define fInitialise, fMove, save_state as pure virtual functions
//       virtual void fInitialise(sdeParticle<sMod>& value, double& logweight,
//                     sdeFilter<sMod> & pf_calcs) = 0;
//       virtual void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
//                     sdeFilter<sMod> & pf_calcs) = 0;
//       virtual void save_state(double *yOut, double *lwgt,
//                     smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> > & Sampler,
//                     sdeParticle<sMod> & pTmp) = 0;
//       virtual ~sdePF() = 0;
// };

// definition of fInitialise template
template <class sMod>
void fInitialise(sdeParticle<sMod>& value, double& logweight,
     sdeFilter<sMod> & pf_calcs) {
  //Rprintf("made it to fInitialise.\n");
  // set particle
  logweight = pf_calcs.init(value.yObs);
  pf_calcs.increase_counter();
  return;
}

// definition of fMove template
template <class sMod>
void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
     sdeFilter<sMod> & pf_calcs) {
  int iCore = 0;
  //Rprintf("counter = %i\n", pf_calcs.get_counter());
  logweight += pf_calcs.update(value.yObs, value.yObs, lTime,
             pf_calcs.get_counter(), iCore);
  pf_calcs.increase_counter();
  return;
}

// definition of save_state template
template <class sMod>
void save_state(double *yOut, double *lwgt,
    smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> > & Sampler,
    sdeParticle<sMod> & pTmp) {
  for(long iPart=0; iPart<Sampler.GetNumber(); iPart++) {
    pTmp = Sampler.GetParticleValueN(iPart);
    pTmp.get_yObs(&yOut[iPart*sMod::nDims]);
    lwgt[iPart] = Sampler.GetParticleLogWeightN(iPart);
  }
  return;
}

// definition of sdeRobj<sMode, sPi>::particleEval
// particleEval is a member function of sdeRobj & sdeCobj
// there is no prior specification for now, we keep sPi only for consistency
template <class sMod, class sPi>
  inline List sdeRobj<sMod, sPi>::particleEval(Numeric initParams,
  	NumericMatrix initData, Numeric dT, Integer nDimsPerObs,
    int nPart, int resample, double dThreshold) {
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  // without NormalDraws we then need users to give nPart
  //int nPart = NormalDraws.nrow()/nDims;
  int nStride = nDims*nPart;
  // for debugging purposes, output whole history
  NumericMatrix dataOut(nDims*nPart, nComp);
  NumericMatrix LogWeightOut(nPart, nComp);
  // pointers to Rcpp memory, i.e., double* representation
  double *yOut = REAL(dataOut); 
  double *lwgt = REAL(LogWeightOut);
  // Sampler can only deep-copy particles out from its storage
  sdeParticle<sMod> pTmp;
  // this is for eventual parallel implementation.
  smc::adaptMethods<sdeParticle<sMod>, sdeFilter<sMod> > *Adapt;
  Adapt = new sdeAdapt<sMod>;

  // determine resample mode
  ResampleType::Enum resampleMode;
  switch (resample) {
    // MULTINOMIAL
    case 0:
      resampleMode = ResampleType::MULTINOMIAL;
      break;
    // RESIDUAL
    case 1:
      resampleMode = ResampleType::RESIDUAL;
      break;
    // STRATIFIED
    case 2:
      resampleMode = ResampleType::STRATIFIED;
      break;
    // SYSTEMATIC
    case 3:
      resampleMode = ResampleType::SYSTEMATIC;
      break;
  }

  //Rprintf("before SMC.\n");
  // SMC
  try {
    smc::sampler<sdeParticle<sMod>, sdeFilter<sMod> >
      Sampler((long)nPart, HistoryType::NONE);
    smc::moveset<sdeParticle<sMod>, sdeFilter<sMod> >
      Moveset(fInitialise<sMod>, fMove<sMod>, NULL);
    //Rprintf("Sampler and Moveset created.\n");

    // AlgParam needs to be deep-copied into Sampler
    // Rprintf("right before SetAlgParams.\n");
    //
    // The sdeFilter constructor can overload without zInit
    // sdeFilter(int np, int nc, 
    //  double *dt, int *yIndex,
    //  double *yInit, double *thetaInit)
    //
    // Also note nPart & nComp are just C++ int objects
    // Don't use INTEGER() to wrap them up
    Sampler.SetAlgParam(sdeFilter<sMod>(nPart, nComp,
          REAL(dT), INTEGER(nDimsPerObs),
          REAL(initData), REAL(initParams)));
    //Rprintf("algParams passed in.\n");
    Sampler.SetResampleParams(resampleMode, dThreshold);
    Sampler.SetMoveSet(Moveset);
    Sampler.Initialise();
    // extract particle from Sampler
    save_state<sMod>(yOut, lwgt, Sampler, pTmp); 
    Sampler.SetAdaptMethods(Adapt); // TBD for parallel processing
    //Rprintf("Sampler initialized.\n");
    // long lTime;
    for(int ii=1; ii<nComp; ii++) {
      //Rprintf("lTime = %i\n", ii);
      Sampler.Iterate();
      // save_state<sMod>(&yOut[ii*(nPart*nDims)],&lwgt[ii*nPart], Sampler, pTmp);
    }
    // only save the last observation (yOut & lgwt)
    save_state<sMod>(&yOut[(nComp-1)*(nPart*nDims)],&lwgt[(nComp-1)*nPart], Sampler, pTmp);
    //delete Sampler;
    delete Adapt;
    // output only the last observation
    // nComp - 1 is the last column since the index starts from 0
    return List::create(Rcpp::Named("X") = dataOut( _ , nComp - 1),
      Rcpp::Named("lwgt") = LogWeightOut( _ , nComp - 1));
    // return List::create(Rcpp::Named("X") = dataOut,
    //   Rcpp::Named("lwgt") = LogWeightOut);
  }
  catch(smc::exception e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue;
}

#endif