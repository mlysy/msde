require(Rcpp)

code <- "
int CoreCount(int n, int nCores) {
  int nBlock = n/nCores+1;
  for(int ii=0; ii<nCores; ii++) {
    int JJstart = ii*nBlock;
    int JJend = (ii+1)*nBlock;
    JJend = JJend < n ? JJend : n;
    Rprintf(\"core = %i, JJstart = %i, JJend = %i\\n\", ii, JJstart, JJend);
    for(int jj=JJstart; jj<JJend; jj++) {
      Rprintf(\"jj = %i\\n\", jj);
    }
  }
  return(0);
}
"

cppFunction(code = code)

CoreCount(8,3)

code <- "
#include <Rcpp.h>
using namespace Rcpp;
double* getPointer(double *x, int i) {
  return &x[i];
}
void foo(double *x, double y) {
  x[0] = y;
  return;
}
//[[Rcpp::export]]
int TestPointer(int n, int i) {
  int ii;
  NumericVector x(n);
  for(ii=0; ii<n; ii++) {
    x(ii) = 1.0*ii;
  }
  foo(getPointer(REAL(x), i), 77.0);
  for(ii=0; ii<n; ii++) {
    Rprintf(\"x[%i] = %f\\n\", ii, x[ii]);
  }
  return 0;
}
"

sourceCpp(code = code, verbose = TRUE, rebuild = TRUE)

TestPointer(5, 3)


#--- scratch ---------------------------------------------------------------

code <- "
sdeModel::nDims
"

evalCpp(code = code, depends = "msdeHeaders",
        includes = "#include \"sdeModel.h\"", verbose = TRUE)
