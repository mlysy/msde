#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends("msde")]]
#include <Interface.h>
#include "Derived.h"

//[[Rcpp::export("construct")]]
SEXP construct() {
  Cobj *foo = new Robj<myCobj>;
  XPtr<Cobj> fooptr(foo, true);
  return fooptr;
}
