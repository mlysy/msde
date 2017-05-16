#' SDE Drift Function
#'
#' @param model An \code{sde.model} object.
#' @param x A vector or matrix of data with \code{ndims} columns.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @return A matrix with \code{ndims} columns containing the drift funtion evaluated at \code{x} and \code{theta}.
#' @export
sde.drift <- function(model, x, theta, debug = FALSE) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # model constants
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize
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
  if(debug) browser()
  nreps <- max(ncol(x), ncol(theta))
  if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  if(!all(c(ncol(x), ncol(theta)) == nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  # compute
  ans <- model$drift(xIn = as.double(x), thetaIn = as.double(theta),
                     nReps = as.integer(nreps))
  dr <- matrix(ans, nreps, ndims, byrow = TRUE)
  dr
}
