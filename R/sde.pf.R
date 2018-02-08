#' Particle filter for SDE model
#' @param model An \code{sde.model} object.
#' @param init Initialization of SDE.
#' @param npart Number of particles.
#' @param resample An integer indicates which resampling scheme the particle filter should use, 0:multinomial; 1:residual; 2:stratified; 3:systematic. Default method is 0: Multinomial.
#' @param threshold A real number between 0 and 1 to indicate the threshold for resampling. By default it is disabled, set to be -0.1
#' @details ...
#' 
#' @return A list with elements:
#' \describe{
#'   \item{\code{Yup}}{A vector of last observation data}
#'   \item{\code{lgwt}}{A vector of normalized log weights corresponding to the last observation.}
#' }
#' 
#' @examples
#' # load pre-compiled model
#' model <- sde.examples("eou")
#'
#' # initial parameters
#' theta0 <- c(alpha = .1, gamma = 1, eta = .3, sigma = .2, rho = -.63)
#' # number of observations
#' nObs <- 100 
#' # number of dimensions
#' nDims <- model$ndims
#' # time between observations (1 year has about 252 trading days)
#' dt <- 1/252 
#' # internal observation time
#' dt.sim <- dt/10
#' # initial SDE values
#' Y0 <- c(X = rnorm(1), V = rnorm(1))
#' # simulate SDE data
#' esim <- sde.sim(model, x0 = Y0, theta = theta0,
#'                 nobs = nObs, # nObs steps forward
#'                 dt = dt, dt.sim = dt/10)
#' # initialization
#' minit <- sde.init(model, x = esim$data, dt = dt, theta = theta0,
#'                  nvar.obs = sample(nDims, nObs, replace = TRUE), m = 1)
#' 
#' # number of particles
#' nPart <- 50
#' # particle filter
#' pf <- sde.pf(model = model, init = minit, npart = nPart, resample = 1, threshold = 0.2)
#' # output the last observation and normalized log-weights
#' X <- pf$X
#' lwgt <- pf$lwgt
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
