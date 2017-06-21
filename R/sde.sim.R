#' Simulation of Multivariate SDE Trajectories
#'
#' @param model An \code{sde.model} object.
#' @param x0 A vector or a matrix of size \code{nreps x ndims} of the SDE values at time 0.
#' @param theta A vector or matrix of size \code{nreps x nparams} of SDE parameters.
#' @param dt Scalar interobservation time.
#' @param dt.sim Scalar, interobservation time for simulation.  That is, interally the interobservation time is \code{dt.sim} but only one out of every \code{dt/dt.sim} simulation steps is kept in the output.
#' @param nobs The number of SDE observations per trajectory to generate.
#' @param burn Scalar burnin value.  Either an integer giving the number of burn-in steps, or a value between 0 and 1 giving the fraction of burn-in relative to \code{nobs}.
#' @param nreps The number of trajectories to generate.
#' @param max.bad.draws The maximum number of times that invalid forward steps are proposed.  See Details.
#' @param verbose Whether or not to display information on the simulation.
#' @details The simulation algorithm is a Markov process with \eqn{Y_0 = x_0} and
#' \deqn{
#' Y_{t+1} \sim N(Y_t + dr(Y_t, \theta) dt_{sim}, df(Y_t, \theta)df(Y_t, \theta)' dt_{sim}).
#' }
#' At each step, a while-loop is used until a valid SDE draw is produced.  The simulation algorithm terminates after \code{nreps} trajectories or once a total of \code{max.bad.draws} are reached.
#' @return A list with elements:
#' \describe{
#'   \item{data}{An array of size \code{nobs x ndims x nreps} containing the simulated SDE trajectories.}
#'   \item{params}{The vector or matrix of parameter values used to generate the data.}
#'   \item{dt,dt.sim}{The actual and internal interobservation times.}
#'   \item{nbad}{The total number of bad draws.}
#' }
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
#' # end of drift example
#' @export
sde.sim <- function(model, x0, theta, dt, dt.sim,
                    nobs, burn = 0, nreps = 1,
                    max.bad.draws = 5e3, verbose = TRUE,
                    debug = FALSE) {
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  # format data and parameters
  ndims <- model$ndims
  x0 <- .format.data(x0, model$data.names, type = "matrix")
  theta <- .format.params(theta, model$param.names)
  # check validity
  N <- c(ncol(x0), ncol(theta))
  single.x <- N[1] == 1
  single.theta <- N[2] == 1
  if(!is.valid.nreps(c(N, nreps))) {
    stop("Incompatible dimensions of x0, theta, and nreps.")
  }
  if(!all(.is.valid.params(model, theta, single.theta, nreps))) {
    stop("theta contains invalid sde parameters.")
  }
  if(!all(.is.valid.data(model, x0, theta, single.x, single.theta, nreps))) {
    stop("x0 contains invalid sde data.")
  }
  # time
  if(dt.sim <= dt) {
    rr <- ceiling(dt/dt.sim)
    dT <- dt/rr
  } else {
    rr <- 1
    dT <- dt
  }
  if(burn < 1) burn <- nobs*burn
  burn <- floor(burn)
  if(verbose) {
    message("Number of normal draws required: ",
            round((nobs+burn-1)*rr*nreps, 2))
    message("Running simulation...")
  }
  if(debug) browser()
  tm <- chrono()
  ans <- model$sim(nDataOut = as.integer(nobs*ndims*nreps),
                   N = as.integer(nobs),
                   burn = as.integer(burn),
                   reps = as.integer(nreps),
                   r = as.integer(rr),
                   dT = as.double(dT),
                   MAXBAD = as.integer(max.bad.draws),
                   initData = as.double(x0),
                   params = as.double(theta),
                   singleX = as.logical(single.x),
                   singleTheta = as.logical(single.theta))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("dataOut", "nBadDraws")
  if(verbose) message("Bad Draws = ", ans$nBadDraws)
  data <- aperm(array(ans$dataOut, dim = c(ndims, nobs, nreps)),
                perm = c(2,1,3))
  dimnames(data) <- list(NULL, model$data.names, NULL)
  out <- list(data = drop(data), params = drop(t(theta)),
              dt = dt, dt.sim = dT,
              nbad = ans$nBadDraws)
  out
}
