#' Particle MCMC sampler for the SDE posterior.
#'
#' Implementation of an independent-coordinate multivariate random-walk Metropolis sampler on the marginal posterior distribution of parameters and missing components of the first SDE observation (see Details).
#'
#' @param model An \code{sde.model} object constructed with \code{\link{sde.make.model}}.
#' @param init An \code{sde.init} object constructed with \code{\link{sde.init}}.
#' @param hyper The hyperparameters of the SDE prior. See \code{\link{sde.prior}}.
#' @param nsamples Number of particle MCMC iterations.
#' @param burn Integer number of burn-in samples, or fraction of \code{nsamples} to prepend as burn-in.
#' @param rw.sd Standard deviation jump size for random walk proposal step size on parameters and missing components of first SDE observation. So the component for theta is before that for missing obs. rw.sd = (sd.theta.component, sd.x.component)
#' @param fixed.params Logical vector of length \code{nparams} indicating which parameters are to be held fixed in the MCMC sampler.
#' @param last.miss.out Logical, whether to return the missing sde components of the last observation.
#' @param npart Number of particles.
#' @param resample The type of particle resampling scheme to use. These are: multi(nomial), resid(ual), strat(ified), sys(tematic).
#' @param threshold A scalar less than 1 to indicate the threshold for resampling. A negative number disables resampling.
#' @return A list of the following elements:
#' \describe{
#'   \item{\code{params}}{An \code{nsamples x nparams} matrix of posterior parameter draws.}
#'   \item{\code{init}}{The \code{sde.init} object which initialized the sampler.}
#'   \item{\code{rw.sd}}{A named vector of RW standard devations used at the last posterior iteration.}
#'   \item{\code{accept}}{A named list of acceptance rates for the various components of the MCMC sampler.}
#' }
sde.pmcmc <- function(model, init, hyper,
                      nsamples, burn, rw.sd = 1, fixed.params,
                      last.miss.out = FALSE, adapt = FALSE,
                      npart, resample, threshold) {
  # model constants
  if(class(model) != "sde.model") {
    stop("model must be an sde.model object.")
  }
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initial values
  if(class(init) != "sde.init") {
    # all argument checking done at construction time
    stop("init must be an sde.init object.")
  }
  if(missing(fixed.params)) fixed.params <- rep(FALSE, nparams)
  # missing data in 1st observation
  nmiss0 <- ndims - init$nvar.obs.m[1]
  # missing data in last observation
  nmissN <- ndims - init$nvar.obs.m[nrow(init$data)] # the number of missing data
  imissN <- ndims-nmissN+1:nmissN # the indices of missing data
  # burn-in
  if(missing(burn)) burn <- max(.1, 1e3)
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # random walk standard deviations tuning params
  tune.par <- .set.jump(rw.sd, adapt, param.names, data.names)
  # output
  Theta <- matrix(NA, nsamples, nparams) # parameters
  colnames(Theta) <- param.names
  X0 <- matrix(NA, nsamples, ndims) # data in first observation
  colnames(X0) <- data.names
  if(last.miss.out) {
    # data in last observation
    XN <- matrix(NA, nsamples, nmissN)
    colnames(XN) <- data.names[imissN]
  }
  # initialize MCMC sampler
  theta.curr <- init$params # it is not init$theta
  x.curr <- init$data[1,]
  # particle filter estimate of log marginal posterior
  lp.marg <- function(lwgt, theta, x) {
    npart <- length(lwgt)
    mx <- max(lwgt)
    lp <- log(sum(exp(lwgt - mx))/npart) + mx # to avoid possible overflow
    lp + sde.prior(model, theta, x, hyper) # add prior
  }
  # initial log-posterior
  pf <- sde.pf(model, init, npart, resample, threshold, history = FALSE)
  lacc.curr <- lp.marg(pf$lwgt, theta.curr, x.curr)
  accept <- 0 # number of accepted draws
  # MCMC ------------------------------------------------
  for(ii in (-burn+1):nsamples) {
    # many of these functions will have to be written
    # I hope the meaning is clear
    # generate a proposal for theta and x in first obs
    # probably easiest for result to have length nparams + ndims
    # with some elements of theta/x never changing if they are fixed
    niter <- ii - (-burn)
    if(tune.par$adapt == TRUE) {
      rw.sd <- adapt.rw(niter, tune.par$sig, tune.par$beta)
    } else {
      rw.sd <- tune.par$sig
    }
    .rw.prop <- function(theta.curr, x.curr, rw.sd, fixed.params, nmiss) {
      # set sd to be 0 for fixed params and observed obs
      # rw.sd = (sd.theta.component, sd.x.component) which is consistent with sde.post
      rw.sd[1:nparams][fixed.params] <- 0
      rw.sd[nparams+1:ndims][1:ndims <= init$nvar.obs.m[1]] <- 0
      # draw full proposal
      rnorm(nparams + ndims, mean = c(theta.curr, x.curr), sd = rw.sd)
    }
    full.prop <- .rw.prop(theta.curr, x.curr, rw.sd, fixed.params, nmiss0)
    theta.prop <- full.prop[1:nparams]
    x.prop <- full.prop[nparams+1:ndims] # or return (x, theta), whichever is more conveninent
    # calculate lacc.prop
    init$params <- theta.prop
    init$data[1,] <- x.prop
    if(sde.valid.params(model, theta.prop) &&
       sde.valid.data(model, x.prop, theta.prop)) {
      # only calculate target density if proposal is valid
      pf <- sde.pf(model, init, npart, resample, threshold, history = FALSE)
      lacc.prop <- lp.marg(pf$lwgt, theta.prop, x.prop)
      lacc <- lacc.prop - lacc.curr
      if(lacc > 0 || runif(1) < exp(lacc)) {
        # proposal is accepted
        x.curr <- x.prop
        theta.curr <- theta.prop
        lacc.curr <- lacc.prop
        accept <- accept + 1
      }
    }
    # storage
    if(ii >= 1) {
      Theta[ii,] <- theta.curr
      X0[ii,] <- x.curr
    }
  }
  return(list(Theta = Theta, X0 = X0, acc = accept/(burn+nsamples)))
  # -----------------------------------------------------
  # # Previous code, kept for reference
  # YObs <- init$data
  # # iteration = 0, initialize the sde model for PMCMC
  # # get initial log-weights & log marginal likelihood via particle filtering
  # pf <- sde.pf(model, init, npart, resample, threshold, history = FALSE)
  # lwgt <- pf$lwgt
  # logYt <- .logY(lwgt, npart)
  # # posterior sampling
  # M <- nsamples
  # nParams <- length(theta0)
  # acc <- rep(NA, M) # accept or reject indicator vector
  # thetaMatrix <- matrix(NA, M+1, nParams)
  # colnames(thetaMatrix) <- model$param.names
  # thetaMatrix[1, ] <- theta0
  # # start iteration of PMCMC
  # for(ii in 2:(M+1)) {
  #   theta_old <- thetaMatrix[ii-1, ]
  #   # random walk proposal for theta
  #   # adaptive sd, increase sd if the last draw was accepted; otherwise, decrease sd
  #   # note ii starts from 2, so the current iteration is the (ii-1)th iteration
  #   if(ii-1 == 1) {
  #     rw.sd <- rw.sd
  #   } else if(ii-1 > 1) {
  #     if(acc[(ii-1)-1] == TRUE) {
  #       rw.sd <- exp(log(rw.sd) + delta/ii)
  #     } else {
  #       rw.sd <- exp(log(rw.sd) - delta/ii)
  #     }
  #   }
  #   # only change the unknown (non-fixed) Lambda1
  #   sd <- rep(NA, nParams)
  #   sd[1:nParams][fixed.params] <- 0
  #   sd[1:nParams][!fixed.params] <- rw.sd
  #   theta_prop <- rnorm(nParams, mean = theta_old, sd = sd)
  #   names(theta_prop) <- model$param.names
  #   # run the particle filter based on theta_prop
  #   # assume no artificial missing points placed between observations, m = 1 by default
  #   tmp_init <- sde.init(model, x = YObs, dt = dT, theta = theta_prop, nvar.obs = 1)
  #   tmp_pf <- sde.pf(model, init = tmp_init, npart = npart,
  #                    resample, threshold, history = FALSE)
  #   # calculate log p_theta_prop(y_T)
  #   lwgt <- tmp_pf$lwgt
  #   logYt_prop <- .logY(lwgt, npart)
  #   # calculate prior densities
  #   logprior_old <- sde.prior(model = model, x = YObs[1,], theta = theta_old, hyper = hyper)
  #   logprior_prop <- sde.prior(model = model, x = YObs[1,], theta = theta_prop, hyper = hyper)
  #   # calculate the log acceptance ratio
  #   # remember we use log density
  #   # we have assumed prior of theta to be 1 for simplicity
  #   logRatio <- logYt_prop - logYt +
  #     logprior_prop - logprior_old +
  #     sum(dnorm(theta_old[!fixed.params], mean = theta_prop[!fixed.params], sd = sd[!fixed.params], log = TRUE)) -
  #     sum(dnorm(theta_prop[!fixed.params], mean = theta_old[!fixed.params], sd = sd[!fixed.params], log = TRUE))
  #   if (logRatio > log(runif(1))) {
  #     # accept the proposal
  #     thetaMatrix[ii, ] <- theta_prop
  #     logYt <- logYt_prop
  #     acc[ii-1] <- TRUE
  #   } else {
  #     thetaMatrix[ii, ] <- theta_old
  #     #logYT <- logYT
  #     acc[ii-1] <- FALSE
  #   }
  # }
  # accRate <- sum(acc)/M
  # out <- list(params = thetaMatrix, accept = accRate,
  #             init = init, rw.sd = rw.sd)
  # return(out)
}

#--- helper functions ----------------------------------------------------------

# Approximation of log marginal likelihood log p(theta | y_T)
# x: Vector of log weights (not incremental weights) at a given time
# nPart: Number of particles
# thataPrior: Prior of theta, by default it is 1
# return: log marginal approximated density log p_theta(theta | y_T)
# .logY <- function(x, nPart, thetaPrior = 1) {
#   # check dimension compatibility
#   if(length(x) != nPart)
#     stop("Check the input vector, its length doesn't match the number of particles")
#   c   <- max(x)
#   ans <- log(sum(exp(x - c))/nPart) + c # to avoid possible overflow
#   ans <- ans * thetaPrior
#   return(ans)
# }

# sig: empirical sd of each variable
adapt.rw <- function(iter, sig, bb = .05) {
  dd <- length(sig)
  if((iter <= 2*dd) || (runif(1) < bb)) {
    rw.sd <- rep(.1/sqrt(dd), dd)
  } else {
    rw.sd <- 2.38/sqrt(dd) * sig
  }
  rw.sd
}
# jump sizes
.set.jump <- function(rw.sd, adapt, param.names, data.names) {
  nparams <- length(param.names)
  ndims <- length(data.names)
  .format.arg <- function(x) {
    if(length(x) == 1) {
      x <- rep(x, nparams+ndims)
      names(x) <- c(param.names, data.names)
    } else if(length(x) == nparams+ndims && is.null(names(x))) {
      names(x) <- c(param.names, data.names)
    }
    x
  }
  # not sure what id is used for
  .set.arg <- function(x, id) {
    y <- rep(0, nparams+ndims)
    names(y) <- c(param.names, data.names)
    y[names(x)] <- x
    y
  }
  rw.sd <- .format.arg(rw.sd)
  if(!is.valid.vars(names(rw.sd), c(param.names, data.names))) {
    stop("names(rw.sd) must be a unique subset of param.names and data.names.")
  }
  rw.sd <- .set.arg(rw.sd)
  if(any(rw.sd < 0)) stop("rw.sd must be non-negative.")
  # adaptive MCMC
  if (adapt) {
    bb <- .05
    adapt <- TRUE
  } else {
    bb <- 0
    adapt <- FALSE
  }
  list(sig = as.double(rw.sd), beta = as.double(bb), adapt = as.logical(adapt))
}

# -------------------------------------------------------------------------------------
# the following helper functions can be removed when we put sde.pmcmc.R under R folder
# -------------------------------------------------------------------------------------
# check if vars are in var.names, non-NULL, no duplicates
is.valid.vars <- function(vars, var.names) {
  valid <- !is.null(vars)
  valid <- valid && !anyDuplicated(vars)
  valid <- valid && all(vars %in% var.names)
  valid
}

# check if hyperparameters are a list of NULL or vector-doubles
is.valid.hyper <- function(phi) {
  valid <- is.list(phi)
  valid <- valid && all(sapply(phi, function(x) {
    is.null(x) || (is.double(x) && is.vector(x))
  }))
  valid
}

# check a vector of lengths and make sure that they
# are all the same or equal to 1
is.valid.nreps <- function(nreps) {
  nreps <- unique(nreps)
  nreps <- c(1, nreps[nreps > 1])
  length(nreps) <= 2
}

# format data for sde.drift, sde.diff, sde.loglik
# output is a matrix or 3-d array with dimensions 1 and 2 permuted
.format.data <- function(x, data.names, type = c("matrix", "array"),
                         strict = FALSE, debug = FALSE) {
  ndims <- length(data.names)
  type <- match.arg(type)
  if(debug) browser()
  # check is.numeric, ndims, and data.names
  if(!is.numeric(x) || !all(is.finite(x))) {
    stop("x must be numeric with finite values.")
  }
  if(type == "matrix") {
    if(is.vector(x)) {
      if(strict) {
        stop("x must be a matrix.")
      } else {
        x <- t(x)
      }
    }
    if(!is.matrix(x)) {
      stop("x must be a vector or matrix.")
    }
    if(ncol(x) != ndims) {
      stop("dimensions of x are incompatible with ndims.")
    }
    if(!is.null(colnames(x)) && !identical(colnames(x), data.names)) {
      stop("names of x do not match data.names.")
    }
    x <- t(x)
  } else {
    if(is.matrix(x)) {
      if(strict) {
        stop("x must be an array.")
      } else {
        dn <- dimnames(x)
        x <- array(x, dim = c(dim(x), 1))
        if(!is.null(dn)) dimnames(x) <- c(dn, list(NULL))
      }
    }
    if(length(dim(x)) != 3) {
      stop("x must be a matrix or 3-d array.")
    }
    if(dim(x)[2] != ndims) {
      stop("dimensions of x are incompatible with ndims.")
    }
    if(!is.null(dimnames(x)) && !identical(dimnames(x)[[2]], data.names)) {
      stop("names of x do not match data.names.")
    }
    x <- aperm(x, c(2,1,3))

  }
  if(is.null(dimnames(x))) dimnames(x) <- rep(list(NULL), length(dim(x)))
  dimnames(x)[[1]] <- data.names
  x
}

# format parameters for sde.drift, sde.diff, sde.loglik
# output is a matrix with dimensions 1 and 2 permuted
.format.params <- function(theta, param.names) {
  nparams <- length(param.names)
  if(is.vector(theta)) theta <- t(theta)
  # check is.numeric, nparams, param.names
  if(!is.numeric(theta) || !all(is.finite(theta))) {
    stop("theta must be a numeric with finite values.")
  }
  if(!is.matrix(theta)) {
    stop("theta must be a vector or matrix.")
  }
  if(ncol(theta) != nparams) {
    stop("dimensions of theta are incompatible with nparams.")
  }
  if(!is.null(colnames(theta)) && !identical(colnames(theta), param.names)) {
    stop("names of theta do not match param.names.")
  }
  if(is.null(dimnames(theta))) {
    dimnames(theta) <- rep(list(NULL), length(dim(theta)))
  }
  dimnames(theta)[[2]] <- param.names
  t(theta)
}
