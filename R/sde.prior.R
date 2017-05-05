#' SDE Prior Function
#'
#' @note TODO: pass \code{fixed.params} and \code{nmiss0} to \code{C++} code.
#' @export
sde.prior <- function(model, theta, x, phi, fixed.params, nmiss0,
                      debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
  if(debug) browser()
  if(!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
  } else {
    x <- t(x)
  }
  if(!is.matrix(theta)) {
    theta <- matrix(theta, ncol = 1)
  } else {
    theta <- t(theta)
  }
  nreps <- max(ncol(x), ncol(theta))
  if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(!ncol(x) == ncol(theta)) {
    stop("x and theta must have the same number of samples.")
  }
  # format hyperparameters
  phi <- model$prior.spec(phi, nparams, ndims, fixed.params, nmiss0)
  # C++ format check (is phi a list with vector-double elements)
  if(!is.valid.hyper(phi)) {
    stop("Unintended behavior.  Please contact package maintainer.")
  }
  # compute
  ans <- model$logprior(thetaIn = as.double(theta), xIn = as.double(x),
                        nReps = as.integer(nreps), phi = phi)
  ans
}
