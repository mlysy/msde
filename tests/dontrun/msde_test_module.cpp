#include <Rcpp.h>
using namespace Rcpp ;
//[[Rcpp::depends(msde)]]
//[[Rcpp::depends(RcppProgress)]]
#include <sdeRobj.h>
#include "sdeModel.h"
#include "sdePrior.h"
typedef sdeRobj<hest::sdeModel, mvn::sdePrior>  sdeRobj_hest_sdeModel_mvn_sdePrior;

RCPP_MODULE(msde_hestModel) {


    class_<sdeRobj<hest::sdeModel, mvn::sdePrior> >("msde_hestModel")

    .constructor()


    .method("nDims", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::get_nDims)
    .method("nParams", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::get_nParams)
    .method("isData", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::isData)
    .method("isParams", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::isParams)
    .method("Drift", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::Drift)
    .method("Diff", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::Diff)
    .method("Loglik", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::LogLik)
    .method("Prior", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::Prior)
    .method("Sim", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::Sim)
    .method("Post", &sdeRobj<hest::sdeModel, mvn::sdePrior> ::Post)
    ;
}
