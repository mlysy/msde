
#include <Rcpp.h>
using namespace Rcpp;
double* getPointer(double *x, int i) {
  return &x[i];
}
void foo(double *x, double y) {
  x[0] = y;
  return;
}
////[[Rcpp::export]]
int TestPointer(int n, int i) {
  int ii;
  NumericVector x(n);
  for(ii=0; ii<n; ii++) {
    x(ii) = 1.0*ii;
  }
  foo(getPointer(REAL(x), i), 77.0);
  for(ii=0; ii<n; ii++) {
    Rprintf("x[%i] = %f\n", ii, x[ii]);
  }
  return 0;
}
