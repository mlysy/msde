#ifndef sdeSMC_h
#define sdeSMC_h 

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
#include "sdePF.h"
// #include "sdeParticle.h"
// #include "sdeAlgPtr.h"
// #include "sdeAdapt.h"

// template <class sMod> 
// class sdePF {
//   public:
//       // define fInitialise, fMove, save_state as pure virtual functions
//       virtual void fInitialise(sdeParticle<sMod>& value, double& logweight,
//                     sdeAlgParams<sMod> & algParams) = 0;
//       virtual void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
//                     sdeAlgParams<sMod> & algParams) = 0;
//       virtual void save_state(double *yOut, double *lwgt,
//                     smc::sampler<sdeParticle<sMod>, sdeAlgParams<sMod> > & Sampler,
//                     sdeParticle<sMod> & pTmp) = 0;
//       virtual ~sdePF() = 0;
// };

// // definition of fInitialise template
// template <class sMod>
// void fInitialise(sdeParticle<sMod>& value, double& logweight,
//      sdeAlgPtr<sMod> & algParams) {
//   //Rprintf("made it to fInitialise.\n");
//   // set particle
//   logweight = algParams.init(value.yObs);
//   algParams.increase_counter();
//   return;
// }

// // definition of fMove template
// template <class sMod>
// void fMove(long lTime, sdeParticle<sMod>& value, double& logweight,
//      sdeAlgPtr<sMod> & algParams) {
//   int iCore = 0;
//   //Rprintf("counter = %i\n", algParams.get_counter());
//   logweight += algParams.update(value.yObs, value.yObs, lTime,
//              algParams.get_counter(), iCore);
//   algParams.increase_counter();
//   return;
// }

// // definition of save_state template
// template <class sMod>
// void save_state(double *yOut, double *lwgt,
//     smc::sampler<sdeParticle<sMod>, sdeAlgPtr<sMod> > & Sampler,
//     sdeParticle<sMod> & pTmp) {
//   // logweights are in fact normalized and divided by nPart
//   // undo this to just get the raw unnormalized weights.
//   double logNC = Sampler.GetLogNCPath() + log(Sampler.GetNumber());
//   for(long iPart=0; iPart<Sampler.GetNumber(); iPart++) {
//     pTmp = Sampler.GetParticleValueN(iPart);
//     pTmp.get_yObs(&yOut[iPart*sMod::nDims]);
//     lwgt[iPart] = Sampler.GetParticleLogWeightN(iPart) + logNC;
//   }
//   return;
// }

// definition of sdeRobj<sMode, sPi>::particleEval
// particleEval is a member function of sdeRobj & sdeCobj
// there is no prior specification for now, we keep sPi only for consistency
template <class sMod, class sPi>
  inline List sdeRobj<sMod, sPi>::particleEval(Numeric initParams,
					       NumericMatrix initData,
					       Numeric dT, Integer nDimsPerObs,
					       int nPart, int resample,
					       double dThreshold,
					       NumericMatrix NormalDraws,
					       bool hasNormalDraws,
					       bool historyOut) {
  int nDims = sMod::nDims;
  int nComp = initData.ncol();
  int nCompOut = historyOut ? nComp : 1;
  // without NormalDraws we then need users to give nPart
  //int nPart = NormalDraws.nrow()/nDims;
  // int nStride = nDims*nPart;
  // for debugging purposes, output whole history
  NumericMatrix dataOut(nDims*nPart, nCompOut);
  NumericMatrix LogWeightOut(nPart, nCompOut);
  // pointers to Rcpp memory, i.e., double* representation
  double *yOut = REAL(dataOut); 
  double *lwgt = REAL(LogWeightOut);
  sdePF<sMod> *PFilter;
  // // Sampler can only deep-copy particles out from its storage
  // sdeParticle<sMod> pTmp;
  // // this is for eventual parallel implementation.
  // // also required when NormalDraws are provided
  // smc::adaptMethods<sdeParticle<sMod>, sdeAlgPtr<sMod> > *Adapt;
  // Adapt = new sdeAdapt<sMod>;
  // // pointer to algParams
  // sdeAlgParams<sMod> *Params;

  // determine resample mode
  ResampleType::Enum resampleMode = static_cast<ResampleType::Enum>(resample);
  // switch (resample) {
  //   // MULTINOMIAL
  //   case 0:
  //     resampleMode = ResampleType::MULTINOMIAL;
  //     break;
  //   // RESIDUAL
  //   case 1:
  //     resampleMode = ResampleType::RESIDUAL;
  //     break;
  //   // STRATIFIED
  //   case 2:
  //     resampleMode = ResampleType::STRATIFIED;
  //     break;
  //   // SYSTEMATIC
  //   case 3:
  //     resampleMode = ResampleType::SYSTEMATIC;
  //     break;
  // }

  //Rprintf("before SMC.\n");
  // SMC
  try {
    // smc::sampler<sdeParticle<sMod>, sdeAlgPtr<sMod> >
    //   Sampler((long)nPart, HistoryType::NONE);
    // smc::moveset<sdeParticle<sMod>, sdeAlgPtr<sMod> >
    //   Moveset(fInitialise<sMod>, fMove<sMod>, NULL);
    // //Rprintf("Sampler and Moveset created.\n");
    // Sampler.SetResampleParams(resampleMode, dThreshold);
    // Sampler.SetMoveSet(Moveset);
    // Sampler.SetAdaptMethods(Adapt);
    // AlgParam needs to be deep-copied into Sampler
    // Rprintf("right before SetAlgParams.\n");
    //
    // The sdeAlgParams constructor can overload without zInit
    // sdeAlgParams(int np, int nc, 
    //  double *dt, int *yIndex,
    //  double *yInit, double *thetaInit)
    //
    // Also note nPart & nComp are just C++ int objects
    // Don't use INTEGER() to wrap them up
    if(hasNormalDraws) {
      PFilter = new sdePF<sMod>(nPart, nComp,
				REAL(dT), INTEGER(nDimsPerObs),
				REAL(initData), REAL(initParams),
				REAL(NormalDraws),
				resampleMode, dThreshold);
      // Params = new sdeAlgParams<sMod>(nPart, nComp,
      // 				      REAL(dT), INTEGER(nDimsPerObs),
      // 				      REAL(initData), REAL(initParams),
      // 				      REAL(NormalDraws));
    } else {
      PFilter = new sdePF<sMod>(nPart, nComp,
				REAL(dT), INTEGER(nDimsPerObs),
				REAL(initData), REAL(initParams),
				resampleMode, dThreshold);
      // Params = new sdeAlgParams<sMod>(nPart, nComp,
      // 				      REAL(dT), INTEGER(nDimsPerObs),
      // 				      REAL(initData), REAL(initParams));
    }
    // Sampler.SetAlgParam(sdeAlgPtr<sMod>(Params));
    //Rprintf("algParams passed in.\n");
    // Sampler.Initialise();
    PFilter->Initialise();
    if(historyOut) {
      // extract particle from Sampler
      // save_state<sMod>(yOut, lwgt, Sampler, pTmp);
      PFilter->eval(yOut, lwgt);
    }
    //Rprintf("Sampler initialized.\n");
    // long lTime;
    for(int ii=1; ii<nComp; ii++) {
      //Rprintf("lTime = %i\n", ii);
      // Sampler.Iterate();
      PFilter->Iterate();
      if(historyOut) {
        // save_state<sMod>(&yOut[ii*(nPart*nDims)],&lwgt[ii*nPart], Sampler, pTmp);
	PFilter->eval(&yOut[ii*(nPart*nDims)],&lwgt[ii*nPart]);
      }
    }
    if(!historyOut) {
      // only save the last observation (yOut & lgwt)
      // save_state<sMod>(yOut,lwgt, Sampler, pTmp);
      PFilter->eval(yOut, lwgt);
    }
    /// deallocate memory
    delete PFilter;
    // delete Sampler;
    // delete Adapt;
    // delete Params;
    // output
    return List::create(Rcpp::Named("data") = dataOut,
			Rcpp::Named("lwgt") = LogWeightOut);
  }
  catch(smc::exception e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue;
}

#endif
