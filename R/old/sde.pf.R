#' Particle filter for SDE model
#' which is currently restricted for eou model and level-1 Euler approximation, later we will make it to be more generalized.
#'
#' Simulate a particle filter based on discretized Euler-Maruyama approximation of a given SDE model
#' @param model An \code{sde.model} object.
#' @param theta A vector or matrix of size \code{nreps x nparams} of SDE parameters.
#' @param Y0 A vector or a matrix of size \code{nreps x ndims} of the SDE values at time 0.
#' @param nObs The number of SDE observations per trajectory to generate.
#' @param dt Scalar interobservation time.
#' @param dt.sim Scalar interobservation time for simulation.  That is, interally the interobservation time is \code{dt.sim} but only one out of every \code{dt/dt.sim} simulation steps is kept in the output.
#' @param nPart Number of particles.
#' @param verbose Whether or not to display information on the simulation.
#' @details ...
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{Yup}}{A matrix of all the historical data}
#'   \item{\code{lgwt}}{A matrix of all the historical normalized log weights.}
#' }
#' @examples
#' # load pre-compiled model
#' emod <- sde.examples("eou")
#'
#' # initialize the sde model
#' theta <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63) # true parameter values
#' nObs <- 100 # number of observations
#' dt <- 1/252 # time between observations (1 year has about 252 trading days)
#' dt.sim <- dt/10
#' Y0 <- c(X = rnorm(1), V = rnorm(1)) # initial SDE values
#' nPart <- 50
#' 
#' output <- sde.pf(emod, theta, Y0, nObs, dt, dt.sim, nPart, verbose = FALSE)
#'
#' @export
sde.pf <- function(model, theta, Y0, nObs, dt, dt.sim, nPart, verbose = FALSE) {
  # model checking
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # initial values checking
  if(is.null(theta)) {
  	message("You didn't give SDE parameters. Default settings will be used.")
    theta <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63)
    cat("The default parameters are:\n", theta)
  }
  if(is.null(Y0)) {
  	message("You didn't give SDE initial values at time 0. Default settings will be used.")
    Y0 <- c(X = rnorm(1), V = rnorm(1))
    cat("The initial SDE values are generated from std normal distributions:\n", Y0)
  }
  if(is.null(nObs)) {
  	message("You didn't give the number of SDE observations. Default value will be used.")
  	nObs <- 100
  	cat("The default number of SDE observations is ", nObs)
  }
  if(is.null(dt)) {
  	message("You didn't give the value of dt. Default value 1/252 will be used.")
  	dt <- 1/252
  }
  if(is.null(dt.sim)) {
  	message("You didn't give the value of dt.sim. Default value dt/10 will be used.")
  	dt.sim <- dt/10
  }
  if(is.null(nPart)) {
  	message("You didn't give the number of particles. Default value 50 will be used.")
  	nPart <- 50
  }

  # initialize the model
  nDims <- model$ndims # number of dimensions
  # simulate some data
  sim <- sde.sim(model, x0 = Y0, theta = theta, nobs = nObs, dt = dt, dt.sim = dt.sim)
  # extract the simulated SDE values (X, V), Yt is an nObs x nDims matrix
  Yt <- sim$data 

  # normal draws
  Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)
  init <- sde.init(model, x = Yt, dt = dt, theta = theta,
  				nvar.obs = 1, # number of observed variables per timepoint/row in data Yt
  				#nvar.obs = sample(nDims, nObs, replace = TRUE),
  				m = 1) # assume no artificial missing points, will be generalized in the future

  # particle filters (in C++ by using SMCTC)
  ans <- .pf_eval(initParams = init$params, initData = t(Yt),
           dT = init$dt.m, nDimsPerObs = init$nvar.obs.m,
           NormalDraws = t(Z))

  out <- list(Yup = t(ans$X), lwgt = t(ans$lwgt))

  return(out)
}