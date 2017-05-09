#' Gibbs sampler for Multivariate SDEs
#'
#' @param model An \code{sde.model} object constructed with \code{\link{sde.make.model}}.
#' @param init.data An \code{sde.data} object constructed with \code{\link{sde.data}}.
#' @param init.params Initial values for all \code{nparams} parameters.
#' @param hyper.params Hyper parameters for the model, to be read by \code{sde.mmodel$prior.spec}.
#' @param nsamples Number of iterations.
#' @param burn Integer number of burnin samples, or fraction of \code{nsamples} to prepend as burnin.
#' @param mwg.sd standard deviation jump size for Metropolis-within-Gibbs steps (see Details).
#' @param adapt Adaptive Metropolis-within-Gibbs sampling (see Details).
#' @param loglik.out Logical, whether to return the loglikelihood at each step.
#' @param last.miss.out Logical, whether to return the missing sde components of the last observation.
#' @param update.data Logical, whether to update the missing data.
#' @param data.out Determines the subset of data to be returned (see Details).
#' @param update.params Logical, whether to update the model parameters.
#' @param verbose Logical, whether to periodically output MCMC status.
#' @param details The MWG jump sizes can be specified as a scalar, a vector or length \code{nparams + ndims}, or a named vector containing the elements defined by \code{nmiss0} and \code{fixed.params}.  The default jump sizes for each MWG random variable are \code{.25 * |initial_value|}.
#'
#'
#' \code{adapt = TRUE} implements an adaptive MCMC by Rosenthal and Roberts (2005).  At step \eqn{n} of the MCMC, the jump size of each MWG random variable is increased or decreased by \eqn{\delta(n)}, depending on whether the cumulative acceptance rate is above or below the optimal value of 0.44.  If \eqn{\sigma_n} is the size of the jump at step \eqn{n}, then the next jump size is determined by
#' \deqn{
#' \log(\sigma_{n+1}) = \log(\sigma_n) \pm \delta(n), \qquad \delta(n) = \min(.01, n^{-1/2}).
#' }
#' When \code{adapt} is not logical, it is a list with elements \code{max} and \code{rate}, such that \code{delta(n) = min(max, 1/n^rate)}.  These elements can be scalars or vectors in the same manner as \code{mwg.sd}.
#'
#' For SDE models with thousands of latent variables, \code{data.out} can be used to thin the MCMC missing data output.  An integer vector or scalar returns specific or evenly-spaced posterior samples from the \code{ncomp x ndims} complete data matrix.  A list
#' @export
sde.post <- function(model, init.data, init.params, fixed.params, hyper.params,
                     nsamples, burn, mwg.sd = NULL, adapt = TRUE,
                     loglik.out = FALSE, last.miss.out = FALSE,
                     update.data = TRUE, data.out,
                     update.params = TRUE, verbose = TRUE, debug = FALSE) {
## sde.post <- function(model, init,
##                      init.data, init.params, par.index, dt, nsamples, burn,
##                      data.out.ind, fixed.params,
##                      prior, rw.jump.sd = NULL, adapt = TRUE,
##                      update.data = TRUE, update.params = TRUE,
##                      loglik.out = FALSE, last.miss.out = FALSE,
##                      verbose = TRUE, debug = FALSE) {
  # model constants
  if(class(model) != "sde.model") {
    stop("model must be of class sde.model.  Use sde.make.model to create.")
  }
  ndims <- model$ndims
  data.names <- model$data.names
  nparams <- model$nparams
  param.names <- model$param.names
  # initialize x and dt
  if(class(init.data) != "sde.data") {
    stop("init.data must be of class sde.data.  Use sde.init to create.")
  }
  dt <- init.data$dt
  par.index <- init.data$par.index
  init.data <- init.data$data
  # initialize parameters
  init.params <- init.data$params
  if(missing(fixed.params)) fixed.params <- rep(nparams, FALSE)
  ## if(!missing(init)) {
  ##   if(!is.null(init$data)) init.data <- init$data
  ##   if(!is.null(init$dt)) dt <- init$dt
  ##   if(!is.null(init$par.index)) par.index <- init$par.index
  ##   if(!is.null(init$params)) init.params <- init$params
  ## }
  # parse inputs
  ncomp <- nrow(init.data)
  nmiss0 <- ndims - par.index[1]
  nparams2 <- nparams+nmiss0
  if(length(dt) == 1) dt <- rep(dt, ncomp-1) # time
  .check.init(model, init.data, dt, init.params)
  if(missing(burn)) burn <- max(.1, 1e3)
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # storage
  # data can be subset both by time and by iteration
  if(all(par.index == ndims)) update.data <- FALSE
  data.out <- .set.data.out(data.out, nsamples, ncomp, update.data)
  nsamples.out <- length(data.out$row)
  ncomp.out <- length(data.out$col)
  ndata.out <- ifelse(update.data, nsamples.out*ndims*ncomp.out, 1)
  nparams.out <- ifelse(update.params, nsamples*nparams, 1)
  nmissN <- ndims - par.index[ncomp] # last missing output
  if(nmissN == 0) last.miss.out <- FALSE
  nlast.miss.out <- ifelse(last.miss.out, nsamples*nmissN, 1)
  ## if(update.data) {
  ##   ndata.out <- nsamples.out*ndims*ncomp.out
  ## } else {
  ##   ndata.out <- 1
  ##   data.out$row <- 1:nsamples
  ## }
  # prior specification
  prior <- hyper.params
  if(debug) browser()
  # format hyperparameters
  prior <- model$prior.spec(prior, nparams, ndims, fixed.params, nmiss0)
  # C++ format check (is phi a list with vector-double elements)
  .check.hyper(prior)
  ## if(!is.valid.hyper(prior)) {
  ##   stop("Unintended behavior.  Please contact package maintainer.")
  ## }
  # random walk jump size
  if(is.null(mwg.sd)) mwg.sd <- .25 * abs(c(init.params, init.data[1,]))
  tune.par <- .set.jump(model, mwg.sd, fixed.params, nmiss0, adapt)
  ## if(is.null(rw.jump.sd)) {
  ##   rw.jump.sd <- abs(init.params)/4
  ##   if(nmiss0 > 0) rw.jump.sd <- c(rw.jump.sd, abs(init.data[par.index[1]+(1:nmiss0)])/4)
  ## }
  ## if(length(rw.jump.sd) < nparams + nmiss0)
  ##   stop("Incorrectly specified random walk jump sd's (need at least nparams + nmiss0).")
  if(verbose) {
    message("Output size:")
    if(update.params) message("params = ", round(nparams.out, 2))
    if(update.data) message("data = ", round(ndata.out, 2))
    message("Running posterior sampler...")
  }
  if(debug) browser()
  tm <- chrono()
  # compute
  ans <- model$post(initParams = as.double(init.params),
                    initData = as.double(t(init.data)),
                    dT = as.double(dt),
                    nDimsPerObs = as.integer(par.index),
                    fixedParams = as.logical(fixed.params),
                    nSamples = as.integer(nsamples),
                    burn = as.integer(burn),
                    nParamsOut = as.integer(nparams.out),
                    nDataOut = as.integer(ndata.out),
                    dataOutRow = as.integer(data.out$row-1),
                    dataOutCol = as.integer(data.out$col-1),
                    updateParams = as.double(update.params),
                    updateData = as.double(update.data),
                    ## priorType = which(prior$type == .PriorTypes),
                    ## priorParams = prior$args,
                    priorArgs = prior,
                    rwJumpSd = as.double(tune.par$sd),
                    updateLogLik = as.integer(loglik.out),
                    nLogLikOut = as.integer(ifelse(loglik.out, nsamples, 1)),
                    updateLastMiss = as.integer(last.miss.out),
                    nLastMissOut = as.integer(nlast.miss.out))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("paramsOut", "dataOut", "paramAccept", "gibbsAccept",
                  "logLikOut", "lastMissOut", "lastIter")
  # acceptance rates
  if(debug) browser()
  accept <- .set.accept(ans$gibbsAccept, ans$paramAccept,
                        nsamples, burn, par.index, nparams, ndims,
                        param.names, data.names, fixed.params,
                        update.data, update.params, verbose)
  ## accept <- NULL
  ## if(update.data) {
  ##   accept <- c(accept, list(data = ans$gibbsAccept/(nsamples+burn)))
  ##   if(verbose) {
  ##     message("Gibbs accept: ", signif(mean(accept$data[par.index < ndims])*100,3), "%")
  ##   }
  ## }
  ## if(update.params) {
  ##   accept <- c(accept, list(params = ans$paramAccept[1:nparams]/(nsamples+burn)))
  ##   if(verbose) {
  ##     for(ii in 1:nparams)
  ##       message(param.names[ii], " accept: ", signif(accept$params[ii]*100,3), "%")
  ##   }
  ## }
  ## if(update.data && (nmiss0 > 0)) {
  ##   accept <- c(accept,
  ##               list(miss0 = ans$paramAccept[nparams+(1:nmiss0)]/(nsamples+burn)))
  ##   for(ii in 1:nmiss0) {
  ##     message(data.names[par.index[1] + ii], "0 accept: ",
  ##             signif(accept$miss0[ii]*100,3), "%")
  ##   }
  ## }
  out <- list()
  if(update.params) {
    theta.out <- matrix(ans$paramsOut,
                        ncol = nparams, byrow = TRUE,
                        dimnames = list(NULL, param.names))
    out <- c(out, list(params = theta.out))
  } else out <- c(out, list(params = init.params))
  if(update.data) {
    x.out <- array(ans$dataOut,
                   dim = c(ndims, ncomp.out, nsamples.out),
                   dimnames = list(data.names, NULL, NULL))
    x.out <- aperm(x.out, perm = 3:1)
    out <- c(out, list(data = x.out))
  } else out <- c(out, list(data = init.data))
  if(loglik.out) out <- c(out, list(loglik = ans$logLikOut))
  out <- c(out, list(dt = dt, par.index = par.index, data.out = data.out,
                     init.data = init.data, init.params = init.params,
                     rw.jump.sd = rw.jump.sd, hyper.params = hyper.params))
  last.iter <- list(params = ans$lastIter[1:nparams],
                    data = matrix(ans$lastIter[nparams + 1:(ncomp*ndims)],
                                  ncomp, ndims, byrow = TRUE))
  names(last.iter$params) <- param.names
  colnames(last.iter$data) <- data.names
  out <- c(out, list(last.iter = last.iter))
  if(last.miss.out) {
    last.miss <- matrix(ans$lastMissOut, nsamples, nmissN, byrow = TRUE)
    colnames(last.miss) <- data.names[par.index[ncomp]+(1:nmissN)]
    out <- c(out, list(last.miss = last.miss))
  }
  out <- c(out, list(accept = accept))
  out
}

# which MCMC iterations and which time points are returned
.set.data.out <- function(data.out, nsamples, ncomp, update.data) {
  if(missing(data.out)) data.out <- 2e3
  if(!is.list(data.out)) {
    # keep all observations
    data.out.row <- data.out
    data.out.col <- 1:ncomp
  } else {
    # which samples and time points
    if(!identical(sort(names(data.out)), sort(c("row", "col")))) {
      stop("data.out must be scalar, vector, or list with elements row and col.")
    }
    data.out.row <- data.out$row
    data.out.col <- data.out$col
  }
  if(length(data.out.row) == 1) {
    # evenly space returned samples
    data.out.row <- unique(floor(seq(1, nsamples, len = data.out.row)))
  }
  # return the data once
  if(!update.data) {
    data.out.row <- 1:nsamples
    data.out.col <- 1:ncomp
  }
  if(is.logical(data.out.row)) data.out.row <- which(data.out.row)
  if(is.logical(data.out.col)) data.out.col <- which(data.out.col)
  list(row = data.out.row, col = data.out.col)
}

# jump sizes
.set.jump <- function(model, sigma0, fixed.params, nmiss0, adapt) {
  nparams <- model$nparams
  ndims <- model$ndims
  .format.arg <- function(x) {
    if(length(x) == 1) {
      x <- rep(x, nparams+ndims)
    } else if(length(x) == nparams+ndims && is.null(names(x))) {
      names(x) <- c(model$param.names, model$data.names)
    }
    x
  }
  sigma0 <- .format.arg(sigma0)
  ## if(length(sigma0) == 1) {
  ##   sigma0 <- rep(sigma0, nparams+ndims)
  ## } else if(length(sigma0) == nparams+ndims && is.null(names(sigma0))) {
  ##   names(sigma0) <- c(model$param.names, model$dim.names)
  ## }
  # check input
  var.id <- sde.active.vars(model = model, vars = names(sigma0),
                            fixed.params = fixed.params, nmiss = nmiss)
  sigma0[!var.id] <- 0
  # adaptive MCMC
  if(is.logical(adapt)) {
    if(adapt) {
      amax <- .01
      arate <- .5
    } else {
      amax <- 0
      arate <- 0
    }
  } else {
    adapt <- TRUE
    amax <- adapt$max
    arate <- adapt$rate
    if(is.null(amax) || is.null(arate)) {
      stop("adapt must be logical or a list with elements max and rate.")
    }
  }
  # check args
  amax <- .format.arg(amax)
  var.id <- sde.active.vars(model = model, vars = names(amax),
                            fixed.params = fixed.params, nmiss = nmiss)
  amax[!var.ind] <- 0
  arate <- .format.arg(arate)
  var.id <- sde.active.vars(model = model, vars = names(arate),
                            fixed.params = fixed.params, nmiss = nmiss)
  arate[!var.ind] <- 0
  list(sd = sigma0, max = amax, rate = arate)
}

# name and dim check for data and params, length check for dt
.check.init <- function(model, init.data, dt, init.params) {
  if(model$ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != model$data.names))
      stop("Incorrect data.names.")
  }
  if(model$nparams != length(init.params))
    stop("init.params does not have the right length.")
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != model$param.names))
      stop("Incorrect param.names.")
  }
  if(length(dt) != nrow(init.data)-1) {
    stop("init.data and dt have incompatible sizes.")
  }
  TRUE
}

# check that phi is a list of double vectors
.check.hyper <- function(phi) {
  if(!is.list(phi)) stop("phi must be a list.")
  is.valid.phi <- sapply(phi, function(x) {
    is.null(x) || (is.double(x) & is.vector(x))
  })
  is.valid.phi <- all(is.valid.phi)
  if(!is.valid.phi) {
    stop("Each element of phi must be NULL or a double vector.")
  }
  is.valid.phi
}

# parse acceptance rates
.set.accept <- function(bb.accept, vnl.accept, nsamples,
                        par.index, nparams, ndims,
                        param.names, data.names, fixed.params,
                        update.data, update.params, verbose) {
  accept <- NULL
  if(update.data) {
    accept <- c(accept, list(data = bb.accept/nsamples))
    if(verbose) {
      message("Gibbs accept: ",
              signif(mean(accept$data[par.index < ndims])*100,3), "%")
    }
  }
  if(update.params) {
    accept <- c(accept, list(params = vnl.accept[1:nparams]/nsamples))
    if(verbose) {
      for(ii in 1:nparams)
        message(param.names[ii], " accept: ",
                signif(accept$params[ii]*100,3), "%")
    }
  }
  nmiss0 <- ndims-par.index[1]
  if(update.data && (nmiss0 > 0)) {
    accept <- c(accept,
                list(miss0 = vnl.accept[nparams+(1:nmiss0)]/nsamples))
    for(ii in 1:nmiss0) {
      message(data.names[par.index[1] + ii], "0 accept: ",
              signif(accept$miss0[ii]*100,3), "%")
    }
  }
  accept
}
