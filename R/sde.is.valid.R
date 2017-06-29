#' Data validator.
#' @param model An \code{sde.model} object.
#' @param x A matrix of data.
#' @param theta A length \code{nparams} vector of parameter values.
#' @examples
#' \donttest{
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#'
#' # example of not valid data
#' x0 <- c(X = log(1000), Z = -0.1)
#' sde.valid.data(model = hmod, x = x0, theta = theta)
#'
#' # example of valid data
#' x0 <- c(X = log(1000), Z = 0.1)
#' sde.valid.data(model = hmod, x = x0, theta = theta)
#' }
#' @export
sde.valid.data <- function(model, x, theta) {
  if(class(model) != "sde.model")
    stop("model must be an sde.model object.")
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
  model$is.data(xIn = as.double(x),
                thetaIn = as.double(theta),
                singleX = as.logical(single.x),
                singleTheta = as.logical(single.theta),
                nReps = as.integer(nreps))
}

#' Parameter validator.
#' @param model An \code{sde.model} object.
#' @param theta A length \code{nparams} vector of parameter values.
#' @examples
#' \donttest{
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#'
#' # example of not valid param
#' theta <- c(alpha = 0.1, gamma = -4, beta = 0.8, sigma = 0.6, rho = -0.8)
#' sde.valid.params(model = hmod, theta = theta)
#'
#' # example of valid param
#' theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
#' sde.valid.params(model = hmod, theta = theta)
#' }
#' @export
sde.valid.params <- function(model, theta) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # initialize
  theta <- .format.params(theta, model$param.names)
  # check singles and compatible x and theta
  nreps <- ncol(theta)
  model$is.params(thetaIn = as.double(theta),
                  nReps = as.integer(nreps))
}

#--- internal versions: no argument checking/formatting ------------------------

.is.valid.data <- function(model, x, theta, single.x, single.theta, nreps) {
  model$is.data(xIn = as.double(x),
                thetaIn = as.double(theta),
                singleX = as.logical(single.x),
                singleTheta = as.logical(single.theta),
                nReps = as.integer(nreps))
}

.is.valid.params <- function(model, theta, single.theta, nreps) {
  rep(model$is.params(thetaIn = as.double(theta),
                      nReps = as.integer(ifelse(single.theta, 1, nreps))),
      times = ifelse(single.theta, nreps, 1))
}
