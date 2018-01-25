#' particle filter prototype for eou model
#'
#' @export
eou.pf <- function(init, npart, Z) {
    ndims <- 2 # for eou model
    nparams <- 5
    # initial values
    if(class(init) != "sde.init") {
        # all argument checking done at construction time
        stop("init must be an sde.init object.")
    }
    # dt <- init$dt.m
    # par.index <- init$nvar.obs.m
    # init.data <- init$data
    # init.params <- init$params
    ncomp <- nrow(init$data)
    # normal draws
    if(missing(Z)) Z <- matrix(rnorm(npart*ndims*(ncomp-1)), ncomp-1, npart*ndims)
    ans <- .pf_eval(initParams = init$params, initData = t(init$data),
                    dT = init$dt.m, nDimsPerObs = init$nvar.obs.m,
                    NormalDraws = t(Z))
    ans$X <- t(ans$X)
    ans$lwgt <- t(ans$lwgt)
    return(ans)
}

# todo:
# 1.  augment this to sde.pf, i.e, any sdeModel.
#     - make sdePF <= particleEval  a virtual function of sdeCobj/sdeRobj
#     - every call to sde.pf allocates/deallocates full memory.
#     - sdeSMC.cpp should become (possible multiple) header files.
# 2.  add smc debug to each model using msde-test_debug.R
#
# 1) I need to change sdeInterface.h
# 2) translate particleEval to a function template and include it in a header file (similar as sdeSim.h/sdePost.h)
