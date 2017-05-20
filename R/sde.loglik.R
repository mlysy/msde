#' SDE Loglikelihood Function
#'
#' @export
sde.loglik <- function(model, x, dt, theta, ncores = 1, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(debug) browser()
  if(is.matrix(x)) {
    ncomp <- nrow(x)
    x <- array(t(x), dim = c(ndims,ncomp,1))
  } else {
    ncomp <- dim(x)[2]
    x <- aperm(x, perm = 3:1)
  }
  if(!is.matrix(theta)) {
    theta <- matrix(theta, ncol = 1)
  } else {
    theta <- t(theta)
  }
  nreps <- max(dim(x)[3], ncol(theta))
  if(dim(x)[3] == 1) x <- array(x, dim = c(ndims,ncomp,nreps))
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(!dim(x)[3] == ncol(theta)) {
    stop("x and theta must have the same number of samples.")
  }
  if(length(dt) == 1) dt <- rep(dt, ncomp-1)
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  # multicore functionality
  if(ncores < 1) stop("ncores must be a positive integer.")
  if(!model$omp && ncores > 1) {
    warning("model not compiled with openMP: ncores set to 1.")
    ncores <- 1
  }
  # compute
  ans <- model$loglik(xIn = as.double(x), dTIn = as.double(dt),
                      thetaIn = as.double(theta),
                      nComp = as.integer(ncomp),
                      nReps = as.integer(nreps),
                      nCores = as.integer(ncores))
  ans
}
