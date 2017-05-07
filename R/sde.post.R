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
#' @param details The MWG jump sizes can be specified as a scalar or vector, the latter with length either \code{(theta0, x0)} or as prescribed by \code{init.data} and \code{fixed.params}.  \code{adapt = TRUE} implements an adaptive MCMC by Rosenthal and Roberts (2005).  At step \eqn{n} of the MCMC, the jump size of each MWG random variable is increased or decreased by \eqn{\delta(n)}, depending on whether the cumulative acceptance rate is above or below the optimal value of 0.44.  If \eqn{\sigma_n} is the size of the jump at step \eqn{n}, then the next jump size is determined by
#' \deqn{
#' \log(\sigma_{n+1}) = \log(\sigma_n) \pm \delta(n), \qquad \delta(n) = \min(.01, n^{-1/2}).
#' }
#' When \code{adapt} is not logical, it is a list with elements \code{max} and \code{rate}, such that \code{delta(n) = min(max, n^rate)}.  These elements can be scalars or vectors in the same manner as \code{mwg.sd}.
#' @export
sde.post <- function(model, init.data, init.params, fixed.params, hyper.params,
                     nsamples, burn, rw.jump.sd = NULL, adapt = TRUE,
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
  if(missing(burn)) burn <- .1
  if(burn < 1) burn <- nsamples*burn
  burn <- floor(burn)
  # storage
  # data can be subset both by time and by iteration
  data.out <- .set.data.out(data.out, nsamples, ncomp)
  nsamples.out <- length(data.out$row)
  ncomp.out <- length(data.out$col)
  if(all(par.index == ndims)) update.data <- FALSE
  if(update.data) {
    ndata.out <- nsamples.out*ndims*ncomp.out
  } else {
    ndata.out <- 1
    data.out$row <- 1:nsamples
  }
  nparams.out <- ifelse(update.params, nsamples*nparams, 1)
  # prior specification
  #prior2 <- prior
  prior <- hyper.params
  if(debug) browser()
  # format hyperparameters
  prior <- model$prior.spec(prior, nparams, ndims, fixed.params, nmiss0)
  # C++ format check (is phi a list with vector-double elements)
  if(!is.valid.hyper(prior)) {
    stop("Unintended behavior.  Please contact package maintainer.")
  }
#--- # random walk jump size ---------------------------------------------------

  if(missing(fixed.params)) fixed.params <- rep(FALSE, nparams)

  #

  if(is.null(rw.jump.sd)) {
    rw.jump.sd <- abs(init.params)/4
    if(nmiss0 > 0) rw.jump.sd <- c(rw.jump.sd, abs(init.data[par.index[1]+(1:nmiss0)])/4)
  }
  if(length(rw.jump.sd) < nparams + nmiss0)
    stop("Incorrectly specified random walk jump sd's (need at least nparams + nmiss0).")
  # last missing output
  nmissN <- ndims - par.index[ncomp]
  if(nmissN == 0) last.miss.out <- FALSE
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
                    rwJumpSd = as.double(rw.jump.sd),
                    updateLogLik = as.integer(loglik.out),
                    nLogLikOut = as.integer(ifelse(loglik.out, nsamples, 1)),
                    updateLastMiss = as.integer(last.miss.out),
                    nLastMissOut = as.integer(ifelse(last.miss.out, nsamples*nmissN, 1)))
  tm <- chrono(tm, display = verbose)
  names(ans) <- c("paramsOut", "dataOut", "paramAccept", "gibbsAccept",
                  "logLikOut", "lastMissOut", "lastIter")
  # acceptance rates
  if(debug) browser()
  accept <- NULL
  if(update.data) {
    accept <- c(accept, list(data = ans$gibbsAccept/(nsamples+burn)))
    if(verbose) {
      message("Gibbs accept: ", signif(mean(accept$data[par.index < ndims])*100,3), "%")
    }
  }
  if(update.params) {
    accept <- c(accept, list(params = ans$paramAccept[1:nparams]/(nsamples+burn)))
    if(verbose) {
      for(ii in 1:nparams)
        message(param.names[ii], " accept: ", signif(accept$params[ii]*100,3), "%")
    }
  }
  if(update.data && (nmiss0 > 0)) {
    accept <- c(accept,
                list(miss0 = ans$paramAccept[nparams+(1:nmiss0)]/(nsamples+burn)))
    for(ii in 1:nmiss0) {
      message(data.names[par.index[1] + ii], "0 accept: ",
              signif(accept$miss0[ii]*100,3), "%")
    }
  }
  out <- list()
  if(update.params) {
    out <- c(out, list(params = matrix(ans$paramsOut,
                                       ncol = nparams, byrow = TRUE,
                                       dimnames = list(NULL, param.names))))
  } else out <- c(out, list(params = init.params))
  if(update.data) {
    out <- c(out, list(data = aperm(array(ans$dataOut,
                                          dim = c(ndims, ncomp.out, nsamples.out),
                                          dimnames = list(data.names, NULL, NULL)), perm = 3:1)))
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

.set.jump <- function(rw.jump.sd, adapt) {
}

format.data.out <- function(data.out, nsamples, ncomp) {
  if(missing(data.out)) data.out <- 2e3
  if(!is.list(data.out)) {
    # keep all observations
    data.out.row <- data.out
    data.out.col <- 1:ncomp
  } else {
    # which samples and
    data.out.row <- data.out$row
    data.out.col <- data.out$col
  }
  if(length(data.out.row) == 1)
    data.out.row <- unique(floor(seq(1, nsamples, len = data.out.row)))
  if(is.logical(data.out.row)) data.out.row <- which(data.out.row)
  if(is.logical(data.out.col)) data.out.col <- which(data.out.col)
  list(row = data.out.row, col = data.out.col)
}

tune.par <- .set.jump(sigma0, theta0, fixed.params, x0, nmiss0, adapt) {
  if(is.null(sigma0)) {
    sigma0 <- abs(init.params)/4
    if(nmiss0 > 0) sigma0 <- c(sigma0, abs(init.data[par.index[1]+(1:nmiss0)])/4)
  }
  if(length(sigma0) < nparams + nmiss0)
    stop("Incorrectly specified random walk jump sd's (need at least nparams + nmiss0).")
}

.check.init(model, init.data, dt, init.params) {
  if(ndims != ncol(init.data))
    stop("init.data does not have the right number of components.")
  if(!is.null(colnames(init.data))) {
    if(any(colnames(init.data) != data.names))
      stop("Incorrect data.names.")
  }
  if(nparams != length(init.params)) stop("init.params does not have the right length.")
  if(!is.null(colnames(init.params))) {
    if(any(colnames(init.params) != param.names))
      stop("Incorrect param.names.")
  }
  # time
  if(length(dt) != ncomp-1) stop("Incorrectly specified dt.")
  invisible(TRUE)
}
