#--- argument checking for the mvn prior ---------------------------------------

mvn.prior.spec <- function(prior.args, nParams, nDims, fixed.params, nMiss0) {
  prior.names <- c("mu", "Sigma")
  if(!is.null(prior.args)) {
    # error checking
    if(is.null(names(prior.args)) ||
       !identical(sort(names(prior.args)), sort(prior.names))) {
      stop("prior.args must be NULL or a list with elements mu and Sigma.")
    }
    if(findInterval(length(prior.args$mu), c(1, nParams+nDims+1)) != 1) {
      stop("mu must have length between 1 and (nParams + nDims).")
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
