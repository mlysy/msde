#ifndef sdeData_h
#define sdeData_h 1

/// Hold data for an SDE model.
///
/// Creates space for the data itself, timepoints,
/// and storage space for mean/variance and normal likelihood evaluations.
///
/// \note A lot of public members here should in fact be protected/private...

template <class sMod>
class sdeData {
 private:
  /// allocate dynamic memory (used by constructors)
  void initialize(int nc, int nmv, int nsde);
 protected:
  int nDims2; ///< number of variance dimensions
 public:
  int nDims; ///< number of dimensions (get from sMod)
  int nParams; ///< number of parameters (get from sMod)
  int nComp; ///< number of (complete data) observations
  double *dT; ///< array of interobservation times (length nComp-1)
  double *sqrtDT; ///< array of sqrt(dT)
  double *XComp; ///< complete data
  int *nObsComp; ///< number of observed dimensions per observation
  double *propMean; ///< storage for sde means
  double *propSd; ///< storage for sde standard deviations
  sMod *sde; ///< storage for drift/diffusion calculations
  double *propZ; ///< storage for normal likelihoods
  /// constructor
  sdeData(int nc, double *dt, int nmv, int nz, int nsde);
  /// destructor
  ~sdeData();
};

template <class sMod>
inline sdeData<sMod>::sdeData(int nc, ///< number of (complete data observations)
			      double *dt, ///< array of interobservation times
			      int nmv, ///< number of mean/variance units to store
			      int nz, ///< number of additional z vectors to store
			      int nsde ///< number of sdeModel objects to store
			      ) {
  nComp = nc;
  nDims = sMod::nDims;
  nDims2 = sMod::diagDiff ? nDims : nDims*nDims;
  nParams = sMod::nParams;
  // create storage space
  dT = new double[nComp];
  sqrtDT = new double[nComp];
  XComp = new double[nComp*nDims];
  propMean = new double[nmv*nDims];
  propSd = new double[nmv*nDims2];
  propZ = new double[nz*nDims];
  sde = new sMod[nsde];
  // assign times
  for(int ii=0; ii<nComp-1; ii++) {
    dT[ii] = dt[ii];
    sqrtDT[ii] = sqrt(dT[ii]);
  }
}

template <class sMod>
inline sdeData<sMod>::~sdeData() {
  delete [] XComp;
  delete [] sde;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] dT;
  delete [] sqrtDT;
}

#endif
