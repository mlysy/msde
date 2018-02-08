#' particle filter prototype for any sde model
#' @param model An \code{sde.model} object.
#' @param init Initialization of SDE.
#' @param npart Number of particles.
#' @param resample An integer indicates which resampling scheme the particle filter should use, 0:multinomial; 1:residual; 2:stratified; 3:systematic. Default method is 0: Multinomial.
#' @param threshold A real number between 0 and 1 to indicate the threshold for resampling. By default it is disabled, set to be -0.1
#' @details ...
#'
#' @export
sde.pf <- function(model, init, npart, resample = 0, threshold = -0.1) {
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

  # run the PF without pre-specified Z
  ans <- .pf_eval(sdeptr = model$ptr, initParams = as.double(init$params), 
            initData = as.matrix(init.data), dT = as.double(dt),
            nDimsPerObs = as.integer(par.index), nPart = npart,
            resample = resample, dThreshold = threshold)
  
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
# 2) translate particleEval to a function template and include it in a header file sdePF.h (similar as sdeSim.h/sdePost.h)
# 3) add particleEval in sdeExports.cpp
# 4) change eou.pf.R into sde.pf.R and update the function body (similar as sde.sim.R & sde.post.R)
