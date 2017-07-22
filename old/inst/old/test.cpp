#include <Rcpp.h>
using namespace Rcpp;

/*
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
    Rprintf("x[%i] = %f\n", ii, x[ii]);
  }
  return 0;
}

*/

class mwgAdapt {
 private:
  double *adaptMax, *adaptRate;
  bool *doAdapt;
  int nRV;
 public:
  mwgAdapt(double *amax, double *arate, int *adapt, int n) {
    nRV = n;
    adaptMax = new double[nRV];
    adaptRate = new double[nRV];
    doAdapt = new bool[nRV];
    for(int ii=0; ii<nRV; ii++) {
      adaptMax[ii] = amax[ii];
      adaptRate[ii] = arate[ii];
      doAdapt[ii] = (adapt[ii] != 0); // R conversion
    }
  }
  ~mwgAdapt() {
    delete [] adaptMax;
    delete [] adaptRate;
    delete [] doAdapt;
  }
  void adapt(double *mwgSd, int *nAccept, int nIter) {
    double acc;
    double lsig;
    double delta;
    const double targAcc = 0.44;
    for(int ii=0; ii<nRV; ii++) {
      if(doAdapt[ii]) {
	acc = (double) nAccept[ii] / (double) nIter;
	delta = pow((double) nIter, -adaptRate[ii]);
	if(delta > adaptMax[ii]) delta = adaptMax[ii];
	lsig = log(mwgSd[ii]);
	lsig += acc < targAcc ? -delta : delta;
	mwgSd[ii] = exp(lsig);
      }
    }
    return;
  }
};

/*
void mwgAdapt(double *mwgSd, double *adaptMax, double *adaptRate,
	      bool *doAdapt, int *nAccept, int nIter, int nRV) {
  double acc;
  double lsig;
  double delta;
  const double targAcc = 0.44;
  for(int ii=0; ii<nRV; ii++) {
    if(doAdapt[ii]) {
      acc = (double) nAccept[ii] / (double) nIter;
      delta = pow((double) nIter, -adaptRate[ii]);
      if(delta > adaptMax[ii]) delta = adaptMax[ii];
      //Rprintf("ii = %i, acc = %f, delta = %f\n", ii, acc, delta);
      lsig = log(mwgSd[ii]);
      lsig += acc < targAcc ? -delta : delta;
      mwgSd[ii] = exp(lsig);
      //Rprintf("delta[%i] = %f, mwgSd[%i] = %f\n", ii, delta, ii, mwgSd[ii]);
    }
  }
  return;
}
*/

//[[Rcpp::export]]
NumericVector mwgAdaptTest(NumericVector mwgSd,
			   NumericVector adaptMax, NumericVector adaptRate,
			   LogicalVector doAdapt, IntegerVector nAccept,
			   int nIter) {
  int nRV = mwgSd.length();
  NumericVector mwgSdOut = Rcpp::clone(mwgSd);
  mwgAdapt tuneMCMC(REAL(adaptMax), REAL(adaptRate), LOGICAL(doAdapt), nRV);
  tuneMCMC.adapt(REAL(mwgSdOut), INTEGER(nAccept), nIter);
  // bool *doAdapt_ = new bool[nRV];
  // for(int ii=0; ii<nRV; ii++) {
  //   doAdapt_[ii] = (bool) doAdapt[ii];
  // }
  // mwgAdapt(REAL(mwgSdOut), REAL(adaptMax), REAL(adaptRate),
  // 	   doAdapt_, INTEGER(nAccept), nIter, nRV);
  // delete [] doAdapt_;
  return mwgSdOut;
}

//[[Rcpp::export]]
int LogicalTest(LogicalVector x) {
  int n = x.length();
  bool *y = reinterpret_cast<bool*>(LOGICAL(x));
  bool *z = new bool[n];
  for(int ii=0; ii<n; ii++) {
    z[ii] = (bool)x[ii];
    Rprintf("x[%i] = %i, y[%i] = %i, z[%i] = %i\n",
	    ii, x[ii], ii, y[ii], ii, z[ii]);
  }
  delete [] z;
  return 0;
}

//[[Rcpp::export]]
int IntegerTest(IntegerVector x) {
  int n = x.length();
  int *y = INTEGER(x);
  for(int ii=0; ii<n; ii++) {
    Rprintf("x[%i] = %i, y[%i] = %i\n", ii, x[ii], ii, y[ii]);
  }
  return 0;
}

//[[Rcpp::export]]
int RealTest(NumericVector x) {
  int n = x.length();
  double *y = REAL(x);
  for(int ii=0; ii<n; ii++) {
    Rprintf("x[%i] = %f, y[%i] = %f\n", ii, x[ii], ii, y[ii]);
  }
  return 0;
}
