#' SDE Loglikelihood Function
#'
#' @param model An \code{sde.model} object.
#' @param x A matrix or 3-d array of data with \code{dim(x)[1]} observations and \code{dim(x)[2] == ndims}.
#' @param dt A scalar or vector of length \code{dim(x)[1]-1} of time intervals between observations.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @param ncores If \code{model} is compiled with \code{OpenMP}, the number of cores to use for parallel processing.  Otherwise, uses \code{ncores = 1} and gives a warning.
#' @return A vector of likelihood evaluations.  If input contains invalid data or parameters an error is thrown.
#' @example
#' # we assume the model's header file is created and it is call, hestModel.h
#' library(msdeHeaders)
#' 
#' # We first initialize the names for data and parameters
#' param.names <- c('alpha','gamma','beta','sigma','rho')
#' data.names <- c("X","Z")
#'
#' # Initialize the model
#' modfile <- "hestModel.h"
#' hmod <- sde.make.model(ModelFile = modfile,
#'                        param.names = param.names,
#'                        data.names = data.names)
#'
#' # Randomly pick the starting point for the process and parameters
#' nDim <- hmod$ndims
#' nparams <- hmod$nparams
#' X0 <- c(X = log(1000), Z = 0.1)
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' dT <- 1/252
#' nobs.sim <- 2000
#' burn <- 500
#' hest.sim <- sde.sim(model = hmod,
#'                     x0 = X0,
#'                     theta = theta,
#'                     dt = dT,
#'                     dt.sim = dT,
#'                     nobs = nobs.sim,
#'                     burn = burn)
#'
#' # ll contains the log-likelihood of Heston model
#' ll <- sde.loglik(model = hmod,
#'                  x = hest.sim$data,
#'                  dt = hest.sim$dt,
#'                  theta = theta)
#'
#' # end of sde.loglik example.
#' @export
sde.loglik <- function(model, x, dt, theta, ncores = 1, debug = FALSE) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  ## # model constants
  ## ndims <- model$ndims
  ## data.names <- model$data.names
  ## nparams <- model$nparams
  ## param.names <- model$param.names
  ## # initialize
  ## if(is.matrix(x)) {
  ##   ncomp <- nrow(x)
  ##   x <- array(t(x), dim = c(ndims,ncomp,1))
  ## } else {
  ##   ncomp <- dim(x)[2]
  ##   x <- aperm(x, perm = 3:1)
  ## }
  ## # fixme: what if ndims = 1?
  ## # probably good to full error check on all of x and theta
  ## # fixme: validate theta's before or in C++ code
  ## if(!is.matrix(theta)) {
  ##   theta <- matrix(theta, ncol = 1)
  ## } else {
  ##   theta <- t(theta)
  ## }
  # FIXME: what is preferred format for x?
  if(debug) browser()
  x <- .format.data(x, model$data.names, type = "array")
  theta <- .format.params(theta, model$param.names)
  # problem dimensions
  ncomp <- dim(x)[2]
  if(ncomp <= 2) {
    stop("likelihood calculation requires at least two observations.")
  }
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) {
    stop("x and dt have incompatible dimensions.")
  }
  nreps <- c(dim(x)[3], ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  ## if(dim(x)[3] == 1) x <- array(x, dim = c(ndims,ncomp,nreps))
  ## if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  ## if(!dim(x)[3] == ncol(theta)) {
  ##   stop("x and theta must have the same number of samples.")
  ## }
  ## if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  ## if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  # multicore functionality
  if(ncores < 1) stop("ncores must be a positive integer.")
  if(!model$omp && ncores > 1) {
    warning("model not compiled with openMP: ncores set to 1.")
    ncores <- 1
  }
  # validate
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!single.theta) {
    theta.ind <- rep(1:nreps, each = ncomp)
  } else {
    theta.ind <- 1
  }
  if(!all(.is.valid.data(model, x, theta[,theta.ind],
                         single.x, single.theta, nreps*ncomp))) {
    stop("x contains invalid sde data.")
  }
  ## if(!all(model$is.params(thetaIn = as.double(theta),
  ##                         nReps = as.integer(ifelse(single.theta, 1
  ##                                                 , nreps))))) {
  ##   stop("theta contains invalid sde parameters.")
  ## }
  ## if(!single.theta) {
  ##   theta.ind <- rep(1:nreps, each = ncomp)
  ## } else {
  ##   theta.ind <- 1
  ## }
  ## if(!all(model$is.data(xIn = as.double(x),
  ##                       thetaIn = as.double(theta[,theta.ind]),
  ##                       singleX = as.logical(single.x),
  ##                       singleTheta = as.logical(single.theta),
  ##                       nReps = as.integer(nreps)))) {
  ##   stop("x contains invalid sde data.")
  ## }
  # compute
  model$loglik(xIn = as.double(x), dTIn = as.double(dt),
               thetaIn = as.double(theta),
               nComp = as.integer(ncomp),
               nReps = as.integer(nreps),
               singleX = as.logical(single.x),
               singleTheta = as.logical(single.theta),
               nCores = as.integer(ncores))
}
#' @example
#' # we assume the model's header file is created and it is call, hestModel.h
#' library(msdeHeaders)
#' 
#' # We first initialize the names for data and parameters
#' param.names <- c('alpha','gamma','beta','sigma','rho')
#' data.names <- c("X","Z")
#'
#' # Initialize the model
#' modfile <- "hestModel.h"
#' hmod <- sde.make.model(ModelFile = modfile,
#'                        param.names = param.names,
#'                        data.names = data.names)
#'
#' # randomly pick the starting point for the process and parameters
#' nDim <- hmod$ndims
#' nparams <- hmod$nparams
#' X0 <- c(X = log(1000), Z = 0.1)
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' dT <- 1/252
#' nobs.sim <- 2000
#' burn <- 500
#' hest.sim <- sde.sim(model = hmod,
#'                     x0 = X0,
#'                     theta = theta,
#'                     dt = dT,
#'                     dt.sim = dT,
#'                     nobs = nobs.sim,
#'                     burn = burn)
#'
#' # intialize argument for sde.init
#' nobs.post <- 128 # observations for posterior
#' nvar.obs <- c(nDim, sample(c(0:nDim),nobs.post-1,replace=TRUE))
#' m <- sample(c(1:5),1) # this can be chosen at random
#'
