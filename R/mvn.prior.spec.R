#--- argument checking for the mvn prior ---------------------------------------

# prior.args is a list with elements \code{mu} and \code{Sigma} giving the mean and variance of components specified by \code{names(mu) == rownames(Sigma) == colnames(Sigma)}.  The remaining components of \code{theta, x0} are given Lebesgue priors.

mvn.prior.spec <- function(prior.args, fixed.params, nmiss0,
                           param.names, data.names) {
  nparams <- length(param.names)
  ndims <- length(data.names)
  prior.names <- c("mu", "Sigma")
  if(!is.null(prior.args)) {
    # error checking
    if(is.null(names(prior.args)) ||
       !identical(sort(names(prior.args)), sort(prior.names))) {
      stop("prior.args must be NULL or a list with elements mu and Sigma.")
    }
    # check argument names
    ## mu.names <- names(prior.args$mu)
    ## Sigma.rnames <- rownames(prior.args$Sigma)
    ## Sigma.cnames <- colnames(prior.args$Sigma)
    ## if(!identical(mu.names, Sigma.rnames) ||
    ##    !identical(mu.names, Sigma.cnames)) {
    ##   stop("names(mu), rownames(Sigma), and colnames(Sigma) are not consistent.")
    ## }
    if(findInterval(length(prior.args$mu), c(1, nparams+ndims+1)) != 1) {
      stop("mu must have length between 1 and (nparams + ndims).")
    }
    if(length(prior.args$mu)^2 != length(prior.args$Sigma)) {
      stop("mu and Sigma have inconsistent dimensions.")
    }
    # format arguments
    prior.args <- prior.args[c("mu", "Sigma")]
    prior.args$Sigma <- chol(prior.args$Sigma)
    prior.args <- lapply(prior.args, as.double)
  } else {
    prior.args <- list(NULL)
  }
  prior.args
}
