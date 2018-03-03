#' Particle MCMC sampling for SDE
#' 
#' @param model An \code{sde.model} object constructed with \code{\link{sde.make.model}}.
#' @param init An \code{sde.init} object constructed with \code{\link{sde.init}}.
#' @param theta0 Initial parameter values.
#' @param fixed.params Logical vector of length \code{nparams} indicating which parameters are to be held fixed in the MCMC sampler.
#' @param nsamples Number of particle MCMC iterations.
#' @param npart Number of particles, fixed integer value.
#' @param dT Scalar interobservation time.
#' @param resample The type of particle resampling scheme to use. These are: multi(nomial), resid(ual), strat(ified), sys(tematic).
#' @param threshold A scalar less than 1 to indicate the threshold for resampling. A negative number disables resampling.
#' @param mwg.sd Standard deviation jump size for Metropolis-within-Gibbs on parameters and missing components of first SDE observation.
#'
#' @return A list of the following elements:
#' \describe{
#'   \item{\code{params}}{An \code{nsamples x nparams} matrix of posterior parameter draws.}
#'   \item{\code{accept}}{A named list of acceptance rates for the various components of the MCMC sampler.}
#'   \item{\code{init}}{The \code{sde.init} object which initialized the sampler.}
#'   \item{\code{mwg.sd}}{A named vector of Metropolis-within-Gibbs standard devations used at the last posterior iteration.}
#' }
#'
#'@export
sde.pmcmc <- function(model, init, theta0, fixed.params, 
                      nsamples, npart, dT,
                      resample = "multi", threshold = 0.5, mwg.sd = 1) {
  # iteration = 0, initialize the sde model for PMCMC
  # get initial log-weights & log marginal likelihood via particle filtering
  pf <- sde.pf(model, init, npart, resample, threshold, history = FALSE)
  lwgt <- pf$lwgt
  logYt <- .logY(lwgt, npart)
  
  # posterior sampling
  M <- nsamples
  nParams <- length(theta0)
  acc <- rep(NA, M) # accept or reject indicator vector
  thetaMatrix <- matrix(NA, M+1, nParams)
  colnames(thetaMatrix) <- model$param.names
  thetaMatrix[1, ] <- theta0
  # start iteration of PMCMC
  for(ii in 2:(M+1)) {
    theta_old <- thetaMatrix[ii-1, ]
    # random walk proposal for theta
    # only change the unknown (non-fixed) Lambda1
    sd <- rep(NA, nParams)
    sd[1:nParams][fixed.params] <- 0
    sd[1:nParams][!fixed.params] <- mwg.sd
    theta_prop <- rnorm(nParams, mean = theta_old, sd = sd)
    
    names(theta_prop) <- model$param.names
    
    # run the particle filter based on theta_prop
    # assume no artificial missing points placed between observations, m = 1 by default
    tmp_init <- sde.init(model, x = YObs, dt = dT, theta = theta_prop, nvar.obs = 1)
    tmp_pf <- sde.pf(model, init = tmp_init, npart = npart,
                     resample, threshold, history = FALSE)
    # calculate log p_theta_prop(y_T) 
    lwgt <- tmp_pf$lwgt
    logYt_prop <- .logY(lwgt, npart)
    
    # calculate the log acceptance ratio
    # remember we use log density 
    # we have assumed prior of theta to be 1 for simplicity
    logRatio <- logYt_prop - logYt + 
      sum(dnorm(theta_old[!fixed.params], mean = theta_prop[!fixed.params], sd = sd[!fixed.params], log = TRUE)) -
      sum(dnorm(theta_prop[!fixed.params], mean = theta_old[!fixed.params], sd = sd[!fixed.params], log = TRUE))
    
    if (logRatio > log(runif(1))) {
      # accept the proposal
      thetaMatrix[ii, ] <- theta_prop
      logYt <- logYt_prop
      acc[ii-1] <- 1
    } else {
      thetaMatrix[ii, ] <- theta_old
      #logYT <- logYT
      acc[ii-1] <- 0
    }
  }
  accRate <- sum(acc)/M
  out <- list(params = thetaMatrix, accept = accRate, init = init, mwg.sd = mwg.sd)
  return(out)
}

# ---------- helper functions -------------------------------------------------------
# Approximation of log marginal likelihood log p(theta | y_T)
# x: Vector of log weights (not incremental weights) at a given time
# nPart: Number of particles
# thataPrior: Prior of theta, by default it is 1
# return: log marginal approximated density log p_theta(theta | y_T)
.logY <- function(x, nPart, thetaPrior = 1) {
  # check dimension compatibility
  if(length(x) != nPart)
    stop("Check the input vector, its length doesn't match the number of particles")
  c   <- max(x)
  ans <- log(sum(exp(x - c))/nPart) + c # to avoid possible overflow
  ans <- ans * thetaPrior
  return(ans)
}
