#' particle filter prototype for any sde model
#'
#' @export
sde.pf <- function(model, init, npart, Z) {
  # model constants
  if (class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initial values
  if (class(init) != "sde.init") {
    # all argument checking done at construction time
    stop("init must be an sde.init object.")
  }
  init.params <- init$params
  init.data <- t(init$data)
  dt <- init$dt.m
  par.index <- init$nvar.obs.m

  ncomp <- nrow(init$data)
  # normal draws
  if (missing(Z))
    Z <- matrix(rnorm(npart * ndims * (ncomp - 1)), ncomp - 1, npart * ndims)
    Z <- t(Z)

  # run particle filter
  ans <- .pf_eval(sdeptr = model$ptr, initParams = as.double(init$params), 
            initData = as.matrix(init.data), dT = as.double(dt),
            nDimsPerObs = as.integer(par.index), NormalDraws = as.matrix(Z))
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
# my steps:
# 1) change sdeInterface.h accordingly
# 2) translate particleEval to a function template and include it in a header file (similar as sdeSim.h/sdePost.h)
# 3) modify the function parameter list of particleEval, just as sdeSim.h & sdePost.h
# 4) add particleEval in sdeExports.cpp
# 5) change eou.pf.R into sde.pf.R and update the function body (refer to sde.sim.R & sde.post.R)
