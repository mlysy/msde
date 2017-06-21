#' SDE Diffusion Function
#'
#' @param model An \code{sde.model} object.
#' @param x A vector or matrix of data with \code{ndims} columns.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @return A matrix with \code{ndims^2} columns containing the diffusion funtion evaluated at \code{x} and \code{theta}.  Each row corresponds to the upper triangular cholesky factor of the diffusion matrix.  If either input contains invalid sde data or parameters an error is thrown.
#' @example
#' library(msdeHeaders)
#' modfile <- "hestModel.h"
#' 
#' # Initialize Heston model in C++ using sde.make.model
#' param.names <- c("alpha","gamma","beta","sigma","rho")
#' data.names <- c("X","Z")
#' hmod <- sde.make.model(ModelFile = modfile,
#'                        param.names = param.names,
#'                        data.names = data.names)
#'
#' # Initialize simulated data
#' X0 <- c(X = log(1000), Z = 0.1)
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' dT = 1/252
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
#' # sde.diff will calculate the diffusion function given data
#' dr <- sde.diff(model = hmod,
#'                x = hest.sim$data[1:100,],
#'                theta = theta)
#' # end of sde.diff example
#' @export
sde.diff <- function(model, x, theta, debug = FALSE) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  ## # model constants
  ## ndims <- model$ndims
  ## data.names <- model$data.names
  ## nparams <- model$nparams
  ## param.names <- model$param.names
  ## # initialize
  ## if(!is.matrix(x)) {
  ##   x <- matrix(x, ncol = 1)
  ## } else {
  ##   x <- t(x)
  ## }
  ## if(!is.matrix(theta)) {
  ##   theta <- matrix(theta, ncol = 1)
  ## } else {
  ##   theta <- t(theta)
  ## }
  # initialize
  ndims <- model$ndims
  x <- .format.data(x, model$data.names, type = "matrix")
  theta <- .format.params(theta, model$param.names)
  # check singles and compatible x and theta
  nreps <- c(ncol(x), ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  ## if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  ## if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  ## if(!all(c(ncol(x), ncol(theta)) == nreps)) {
  ##   stop("x and theta have incompatible dimensions.")
  ## }
  if(debug) browser()
  # validate
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, x, theta, single.x, single.theta, nreps))) {
    stop("x contains invalid sde data.")
  }
  ## if(!all(model$is.params(thetaIn = as.double(theta),
  ##                         nReps = as.integer(ifelse(single.theta, 1
  ##                                                 , nreps))))) {
  ##   stop("theta contains invalid sde parameters.")
  ## }
  ## if(!all(model$is.data(xIn = as.double(x),
  ##                       thetaIn = as.double(theta),
  ##                       singleX = as.logical(single.x),
  ##                       singleTheta = as.logical(single.theta),
  ##                       nReps = as.integer(nreps)))) {
  ##   stop("x contains invalid sde data.")
  ## }
  ## val <- model$is.params(thetaIn = as.double(theta),
  ##                        nReps = as.integer(ifelse(single.theta, 1, nreps)))
  ## if(single.theta) val <- rep(val, nreps)
  ## val <- val & model$is.data(xIn = as.double(x),
  ##                            thetaIn = as.double(theta),
  ##                            singleX = as.logical(single.x),
  ##                            singleTheta = as.logical(single.theta),
  ##                            nReps = as.integer(nreps))
  # compute
  ## df <- matrix(NA, ndims^2, nreps)
  ## if(any(val)) {
  ##   ans <- model$diff(xIn = as.double(x[,val]),
  ##                     thetaIn = as.double(theta[,val]),
  ##                     singleX = as.logical(single.x),
  ##                     singleTheta = as.logical(single.theta),
  ##                     nReps = as.integer(sum(val)))
  ##   df[,val] <- ans
  ## }
  ## # put zeros into the off-triangular elements
  ## df[val,lower.tri(diag(ndims))] <- 0
  df <- model$diff(xIn = as.double(x),
                   thetaIn = as.double(theta),
                   singleX = as.logical(single.x),
                   singleTheta = as.logical(single.theta),
                   nReps = as.integer(nreps))
  df <- matrix(df, nrow = nreps, ncol = ndims^2, byrow = TRUE)
  # put zeros into the off-triangular elements
  df[,lower.tri(diag(ndims))] <- 0
  df
}
