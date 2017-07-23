// generic R wrappers to exported members of instantiated Cobj

#include <Rcpp.h>
using namespace Rcpp;
#include "Interface.h"

//[[Rcpp::export(".Cobj_eval")]]
double Cobj_eval(SEXP Cobj_ptr, double y) {
  XPtr<Cobj> foo(Cobj_ptr);
  return foo->eval(y);
}

//[[Rcpp::export(".Cobj_setx")]]
void Cobj_setx(SEXP Cobj_ptr, double x) {
  XPtr<Cobj> foo(Cobj_ptr);
  foo->set_x(x);
  return;
}

//[[Rcpp::export(".Cobj_getx")]]
double Cobj_getx(SEXP Cobj_ptr) {
  XPtr<Cobj> foo(Cobj_ptr);
  return foo->get_x();
}
