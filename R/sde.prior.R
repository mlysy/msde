#' SDE Prior Function
#'
#' @param model An \code{sde.model} object.
#' @param theta A vector or matrix of parameters with \code{nparams} columns.
#' @param x A vector or matrix of data with \code{ndims} columns.
#' @param hyper The hyperparameters of the SDE prior.  See Details.
#' @return A vector of log-prior densities evaluated at the inputs.
#' @details The prior is constructed at the \code{C++} level by defining a function (i.e., public member) \code{double logPrior(double *theta, double *x)} within the \code{sdePrior} class (see Examples).  At the \code{R} level, the \code{prior.spec} argument of \code{sde.make.model} is a function with arguments \code{hyper}, \code{param.names}, \code{data.names} used to convert \code{hyper} into a list of \code{NULL} or double-vectors which get passed on to the \code{C++} code.  This function can also be used to throw \code{R}-level errors to protect the \code{C++} code from invalid inputs, as is done for the default prior in \code{\link{mvn.prior.spec}}.
#' @examples
#' \donttest{
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # setting prior for 3 parameters
#' rv.names <- c("alpha","gamma","rho")
#' mu <- rnorm(3)
#' Sigma <- crossprod(matrix(rnorm(9),3,3))
#' names(mu) <- rv.names
#' colnames(Sigma) <- rv.names
#' rownames(Sigma) <- rv.names
#' prior.args <- list(mu = mu, Sigma = Sigma)
#'
#' # Simulate data
#' nreps-10
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' x0 <- c(X = log(1000), Z = 0.1)
#' Theta <- apply(t(replicate(nReps,theta)),2,jitter)
#' X0 <- apply(t(replicate(nReps,x)),2,jitter)
#'
#' sde.prior(model = hmod, x = X0, theta = Theta, hyper = prior.args)
#' }
#' @export
sde.prior <- function(model, theta, x, hyper) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # initialize
  param.names <- model$param.names
  data.names <- model$data.names
  x <- .format.data(x, data.names, type = "matrix")
  theta <- .format.params(theta, param.names)
  # check singles and compatible x and theta
  nreps <- c(ncol(x), ncol(theta))
  single.x <- nreps[1] == 1
  single.theta <- nreps[2] == 1
  if(!is.valid.nreps(nreps)) {
    stop("x and theta have incompatible dimensions.")
  }
  nreps <- max(nreps)
  ## # model constants
  ## ndims <- model$ndims
  ## data.names <- model$data.names
  ## nparams <- model$nparams
  ## param.names <- model$param.names
  ## # initialize
  ## if(debug) browser()
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
  ## nreps <- max(ncol(x), ncol(theta))
  ## if(ncol(x) == 1) x <- matrix(x, ndims, nreps)
  ## if(ncol(theta) == 1) theta <- matrix(theta, nparams, nreps)
  ## if(!ncol(x) == ncol(theta)) {
  ##   stop("x and theta must have the same number of samples.")
  ## }
  # format hyperparameters
  phi <- model$prior.spec(hyper = hyper,
                          param.names = param.names, data.names = data.names)
  # C++ format check (is phi a list with vector-double elements)
  if(!is.valid.hyper(phi)) {
    stop("model$prior.spec must convert hyper to a list with NULL or double-vector elements.")
  }
  # compute
  ans <- model$logprior(thetaIn = as.double(theta), xIn = as.double(x),
                        singleTheta = as.logical(single.theta),
                        singleX = as.logical(single.x),
                        nReps = as.integer(nreps), phiIn = phi)
  ans
}
