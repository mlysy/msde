# Particle MCMC sampling for SDE (currently not exported)
# 
# @param model An \code{sde.model} object constructed with \code{\link{sde.make.model}}.
# @param init An \code{sde.init} object constructed with \code{\link{sde.init}}.
# @param theta0 Initial parameter values.
# @param fixed.params Logical vector of length \code{nparams} indicating which parameters are to be held fixed in the MCMC sampler.
# @param hyper The hyperparameters of the SDE prior. For details, please see \code{\link{sde.prior}}.
# @param nsamples Number of particle MCMC iterations.
# @param npart Number of particles, fixed integer value.
# @param dT Scalar interobservation time.
# @param resample The type of particle resampling scheme to use. These are: multi(nomial), resid(ual), strat(ified), sys(tematic).
# @param threshold A scalar less than 1 to indicate the threshold for resampling. A negative number disables resampling.
# @param rw.sd Standard deviation jump size for random walk proposal step size on parameters and missing components of first SDE observation.
# @param delta Tuning parameter for adaptive rw.sd. It has the same length as rw.sd
#
# @return A list of the following elements:
# \describe{
#   \item{\code{params}}{An \code{nsamples x nparams} matrix of posterior parameter draws.}
#   \item{\code{accept}}{A named list of acceptance rates for the various components of the MCMC sampler.}
#   \item{\code{init}}{The \code{sde.init} object which initialized the sampler.}
#   \item{\code{rw.sd}}{A named vector of RW standard devations used at the last posterior iteration.}
#   \item{\code{delta}}{A numeric value or a vector of tuning parameter(s) for rw.sd}
# }
#
sde.pmcmc <- function(model, init, theta0, fixed.params, hyper, 
                      nsamples, npart, dT,
                      resample = "multi", threshold = 0.5, rw.sd = 1, delta = .618) {
  YObs <- init$data
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
    # adaptive sd, increase sd if the last draw was accepted; otherwise, decrease sd
    # note ii starts from 2, so the current iteration is the (ii-1)th iteration
    if(ii-1 == 1) {
      rw.sd <- rw.sd
    } else if(ii-1 > 1) {
      if(acc[(ii-1)-1] == TRUE) {
        rw.sd <- exp(log(rw.sd) + delta/ii)
      } else {
        rw.sd <- exp(log(rw.sd) - delta/ii)
      }
    }
    # only change the unknown (non-fixed) Lambda1
    sd <- rep(NA, nParams)
    sd[1:nParams][fixed.params] <- 0
    sd[1:nParams][!fixed.params] <- rw.sd
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
    # calculate prior densities
    logprior_old <- sde.prior(model = model, x = YObs[1,], theta = theta_old, hyper = hyper)
    logprior_prop <- sde.prior(model = model, x = YObs[1,], theta = theta_prop, hyper = hyper)
    # calculate the log acceptance ratio
    # remember we use log density 
    # we have assumed prior of theta to be 1 for simplicity
    logRatio <- logYt_prop - logYt + 
      logprior_prop - logprior_old +
      sum(dnorm(theta_old[!fixed.params], mean = theta_prop[!fixed.params], sd = sd[!fixed.params], log = TRUE)) -
      sum(dnorm(theta_prop[!fixed.params], mean = theta_old[!fixed.params], sd = sd[!fixed.params], log = TRUE))
    if (logRatio > log(runif(1))) {
      # accept the proposal
      thetaMatrix[ii, ] <- theta_prop
      logYt <- logYt_prop
      acc[ii-1] <- TRUE
    } else {
      thetaMatrix[ii, ] <- theta_old
      #logYT <- logYT
      acc[ii-1] <- FALSE
    }
  }
  accRate <- sum(acc)/M
  out <- list(params = thetaMatrix, accept = accRate, init = init, rw.sd = rw.sd)
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
