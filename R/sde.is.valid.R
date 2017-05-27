#' Data validator
#'
#' @export
sde.valid.data <- function(model, x, theta) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
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

#' Parameter validator
#'
#' @export
sde.valid.params <- function(model, theta) {
  if(class(model) != "sde.model")
    stop("Expecting object of class sde.model.  Use sde.make.model to create.")
  # initialize
  theta <- .format.params(theta, model$param.names)
  # check singles and compatible x and theta
  nreps <- ncol(theta)
  model$is.params(xIn = as.double(x),
                  thetaIn = as.double(theta),
                  nReps = as.integer(nreps))
}
