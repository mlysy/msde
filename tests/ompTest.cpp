#include <Rcpp.h>
using namespace Rcpp;
#include <omp.h>

class foo {
  int nCores;
public:
  foo(int ncores) {
    nCores = ncores;
  }
  void print(int n) {
    omp_set_num_threads(nCores);
    Rprintf("omp_get_num_threads = %i\n", omp_get_num_threads());
    #pragma omp parallel for
    for(int ii=0; ii<n; ii++) {
      Rprintf("printing %i on thread %i of %i\n",
	      ii, omp_get_thread_num(), omp_get_num_threads());
    }
    return;
  }
};

//[[Rcpp::export]]
int ompTest(int nCores) {
  // omp_set_num_threads(nCores);
  // Rprintf("omp_get_num_threads = %i\n", omp_get_num_threads());
  // #pragma omp parallel for
  // for(int ii=0; ii<10; ii++) {
  //   Rprintf("printing %i on thread %i of %i\n",
  // 	    ii, omp_get_thread_num(), omp_get_num_threads());
  // }
  foo bar(nCores);
  bar.print(10);
  return 0;
}
