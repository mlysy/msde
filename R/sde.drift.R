#' SDE Drift Function
#'
#' @param model An \code{sde.model} object.
#' @param x A vector or matrix of data with \code{ndims} columns.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @return A matrix with \code{ndims} columns containing the drift funtion evaluated at \code{x} and \code{theta}.  If either input contains invalid sde data or parameters an error is thrown.
#' @examples
#' \donttest{
#' # compile model
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # single input
#' x0 <- c(X = log(1000), Z = 0.1)
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' sde.drift(model = hmod, x = x0, theta = theta)
#'
#' # multiple inputs
#' nreps <- 10
#' Theta <- apply(t(replicate(nreps,theta)),2,jitter)
#' X0 <- apply(t(replicate(nreps,x0)),2,jitter)
#' sde.drift(model = hmod, x = X0, theta = Theta)
#' }
#' @export
sde.drift <- function(model, x, theta) {
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
  ## if(debug) browser()
  ## nreps <- max(ncol(x), ncol(theta))
  ## if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  ## if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  ## if(!all(c(ncol(x), ncol(theta)) == nreps)) {
  ##   stop("x and theta have incompatible dimensions.")
  ## }
  # initialize
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
#  if(debug) browser()
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
  ## dr <- matrix(NA, ndims, nreps)
  ## if(any(val)) {
  ##   ans <- model$drift(xIn = as.double(x[,val]),
  ##                      thetaIn = as.double(theta[,val]),
  ##                      singleX = as.logical(single.x),
  ##                      singleTheta = as.logical(single.theta),
  ##                      nReps = as.integer(sum(val)))
  ##   dr[,val] <- ans
  ## }
  dr <- model$drift(xIn = as.double(x),
                    thetaIn = as.double(theta),
                    singleX = as.logical(single.x),
                    singleTheta = as.logical(single.theta),
                    nReps = as.integer(nreps))
  matrix(dr, nrow = nreps, ncol = model$ndims, byrow = TRUE)
}
