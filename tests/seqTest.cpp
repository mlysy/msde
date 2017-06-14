#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
int seqTest(IntegerVector parIndex, int nDims) {
  int ii, jj, II, JJ;
  int nComp = parIndex.length();
  int *nObsComp = INTEGER(parIndex);
  // missing data
  int nMiss0 = nDims-nObsComp[0]; // unobserved states in first observation
  int nMissN = nDims-nObsComp[nComp-1]; // unobserved states in last observation
  // identify missing _intermediate_ data indices,
  // i.e. at least one component to update and not first or last observation
  int nMiss = 0;
  for(ii = 1; ii < nComp-1; ii++) {
    nMiss += (nObsComp[ii] < nDims);
  }
  int *missInd = new int[nMiss + (nMiss == 0)];
  jj = 0;
  for(ii = 1; ii < nComp-1; ii++) {
    if(nObsComp[ii] < nDims) {
      missInd[jj++] = ii;
    }
  }
  // Rprintf("Before split:\n");
  // for(ii=0; ii<nMiss; ii++) {
  //   Rprintf("missInd[%i] = %i\n", ii, missInd[ii]);
  // }
  // print double sequence
  for(JJ=0; JJ<2; JJ++) {
    Rprintf("JJ = %i\n", JJ);
    for(II=JJ; II<nMiss; II = II+2) {
      Rprintf("missInd[%i] = %i\n", II, missInd[II]);
    }
  }
  delete [] missInd;
  return 0;
}
