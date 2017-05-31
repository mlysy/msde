#' MCMC Initialization
#'
#' Specifies the observed SDE data, interobservation times, and initial parameter and missing data values for \code{\link{sde.post}}.
#' @param x An \code{nobs x ndims} matrix of data.
#' @param dt A scalar or length \code{nobs-1} vector of interobservations times.
#' @param theta A length \code{nparams} vector of parameter values.
#' @param nvar.obs A scalar or length \code{nobs} vector of integers between 0 and \code{ndims} denoting the number of observed SDE variables in each row of \code{data}.  See Details.
#' @param m Positive integer, such that \code{m-1} evenly-spaced missing data time points are placed between observations.  See Details.
#' @return An \code{sde.init} object.  A list with elements:
#' \describe{
#'   \item{data}{An \code{ncomp x ndims} matrix of complete data, where \code{ncomp = N_m = m * (nobs-1)+1}.}
#'   \item{dt.m}{The complete data interobservation time, \code{dt_m = dt/m}.}
#'   \item{nvar.obs.m}{The number of variables observed per row of \code{data}.  Note that \code{nvar.obs.m[(i-1)*m+1] == nvar.obs[ii]}, and that \code{nvar.obs.m[i-1] == 0} if \code{i} is not a multiple of \code{m}.}
#'   \item{params}{Parameter initial values.}
#' }
#' @export
sde.init <- function(x, dt, m = 1, nvar.obs, theta, debug = FALSE) {
  nobs <- nrow(x)
  ndims <- ncol(x)
  #mm <- m
  #m <- m-1
  #if(missing(m)) m <- 2^k-1
  ncomp <- (nobs-1)*mm+1
  init.data <- matrix(NA, ncomp, ndims)
  colnames(init.data) <- colnames(x)
  # interpolation to create missing data
  if(debug) browser()
  dt1 <- length(dt) == 1
  if(dt1) dt <- rep(dt, nobs-1)
  if(length(dt) != nobs-1) stop("x and dt have incompatible sizes.")
  dtnew <- rep(dt/mm, each = mm)
  told <- cumsum(c(0, dt))
  tnew <- cumsum(c(0, dtnew))
  for(ii in 1:ndims) {
    init.data[,ii] <- approx(x = told, y = x[,ii],
                             xout = tnew)$y
  }
  if(dt1) dtnew <- dtnew[1]
  #par.index <- nvar.obs
  #if(missing(par.index)) par.index <- ndims
  par.ind <- rep(0, ncomp)
  par.ind[seq(1, ncomp, len = nobs)] <- nvar.obs
  ans <- list(data = init.data, dt,m = dtnew,
              nvar.obs.m = par.ind, params = theta)
  #if(!missing(params)) ans$params <- params
  class(ans) <- "sde.init"
  ans
}
