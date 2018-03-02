#' Particle filter for SDE model.
#'
#' @param model An \code{sde.model} object.
#' @param init Initialization of SDE.
#' @param npart Number of particles.
#' @param resample The type of particle resampling scheme to use. These are: multi(nomial), resid(ual), strat(ified), sys(tematic).
#' @param threshold A scalar less than 1 to indicate the threshold for resampling. A negative number disables resampling.
#' @param Z Optional array of dimensions \code{(ncomp - 1) x ndims x npart} providing the standard normal draws for the filter to use.  This is most useful for debugging, in conjunction with setting \code{threshold} to a negative value.
#' @param history Logical; whether or not the entire history of the particle filter should be output, or only draws for the last observation.
#' @details ...
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{data}}{If \code{history = FALSE}, a \code{npart x ndims} matrix of particles for the last observation.  If \code{history = TRUE}, an array of dimension \code{npart x ndims x nComp}.}
#'   \item{\code{lgwt}}{If \code{history = FALSE}, a \code{npart}-length vector of normalized log weights for the last observation.  Otherwise a matrix of dimension \code{npart x nComp}.}
#' }
#'
#' @examples
#' # load pre-compiled model
#' emod <- sde.examples("eou")
#'
#' # initial parameters
#' theta0 <- c(alpha = .1, gamma = 1, eta = .3, sigma = .2, rho = -.63)
#' # number of observations
#' nObs <- 100
#' # number of dimensions
#' nDims <- emod$ndims
#' # time between observations (1 year has about 252 trading days)
#' dt <- 1/252
#' # internal observation time
#' dt.sim <- dt/10
#' # initial SDE values
#' Y0 <- c(X = rnorm(1), V = rnorm(1))
#' # simulate SDE data
#' esim <- sde.sim(emod, x0 = Y0, theta = theta0,
#'                 nobs = nObs, # nObs steps forward
#'                 dt = dt, dt.sim = dt/10)
#' # initialization
#' minit <- sde.init(emod, x = esim$data, dt = dt, theta = theta0,
#'                  nvar.obs = sample(nDims, nObs, replace = TRUE), m = 1)
#'
#' # number of particles
#' nPart <- 37
#' Z <- array(rnorm(nPart*nDims*(nObs-1)), c(nObs-1, nDims, nPart))
#' # particle filter (without pre-specified Z)
#' pf <- sde.pf(emod, init = minit, npart = nPart,
#'              resample = "multi", threshold = -1,
#'              Z = Z, history = TRUE)
#' # output the last observation and normalized log-weights
#' data <- pf$data
#' lwgt <- pf$lwgt
#'
#' @export
sde.pf <- function(model, init, npart,
                   resample = c("multi", "resid", "strat", "sys"),
                   threshold = 0.5, Z, history = FALSE) {
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
  # code resample into integer value
  resample <- match.arg(resample)
  resample <- switch(resample,
                     multi = 0L, resid = 1L, strat = 2L, sys = 3L)

  # hasZ == TRUE when Z is provided; otherwise hasZ == FALSE
  hasZ <- !missing(Z)
  if(hasZ) {
    # check the dimensions of Z input
    if(!all(dim(Z) == c(ncomp-1, ndims, npart))) {
 +    stop("Z must be an array of dimensions (ncomp-1) x ndims x npart.")
 +  }
    # transform the input 3-d array Z into a matrix
    dim(Z) <- c((ncomp-1), ndims*npart)
    Z <- t(Z)
  } else {
    Z <- as.matrix(0)
  }
  # run the particle filter
  ans <- .pf_eval(sdeptr = model$ptr, initParams = as.double(init$params),
                  initData = as.matrix(init.data), dT = as.double(dt),
                  nDimsPerObs = as.integer(par.index), nPart = npart,
                  resample = as.integer(resample), dThreshold = threshold,
                  NormalDraws = Z, hasNormalDraws = hasZ,
                  historyOut = history)

  # output 
  # ans$data should also be a 3-d array if history == TRUE; otherwise, a matrix
  if(history == TRUE) {
    ans$data <- array(c(ans$data), dim = c(ndims, npart, ncomp))
    ans$data <- aperm(ans$data, perm = c(2, 1, 3)) # notice the option is "perm" not "dim"
    dimnames(ans$data)[[2]] <- data.names
  } else {
    ans$data <- matrix(ans$data, nrow = ndims, ncol = npart)
    ans$data <- t(ans$data)
  }
  ans$lwgt <- t(ans$lwgt)

  out <- list(data = ans$data, lwgt = ans$lwgt)
  return(out)
}
