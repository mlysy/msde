// R object constructor

#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::depends("msde")]]

#include "sdeModel.h"
#include "sdePrior.h"
#include <sdeInterface.h>

//[[Rcpp::export(".sdeMakeModel")]]
SEXP sdeMakeModel() {
  sdeCobj *sde = new sdeRobj<sdeModel, sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}
