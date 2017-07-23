// export sample models

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeInterface.h"
#include "mvnPrior_pkg.h"
#include "hestModel_pkg.h"
#include "lotvolModel_pkg.h"
#include "biouModel_pkg.h"
#include "pgnetModel_pkg.h"

//[[Rcpp::export(".hest_MakeModel")]]
SEXP hestMakeModel() {
  sdeCobj *sde = new sdeRobj<hestModel, mvnPrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".biou_MakeModel")]]
SEXP biouMakeModel() {
  sdeCobj *sde = new sdeRobj<biouModel, mvnPrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".pgnet_MakeModel")]]
SEXP pgnetMakeModel() {
  sdeCobj *sde = new sdeRobj<pgnetModel, mvnPrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".lotvol_MakeModel")]]
SEXP lotvolMakeModel() {
  sdeCobj *sde = new sdeRobj<lotvolModel, mvnPrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}
