# helper functions for pmcmc

# approximation of log conditional likelihood p_theta(y_n | y_n-1)
#'@param x Vector of log-weights at a given time
#'@param nPart Number of particles
#'@return log marginal approximated density log p_theta(y_n | y_n-1)
logCondY <- function(x, nPart) {
    # check dimension compatibility
    if(length(x) != nPart)
        stop("Check the input vector, its length doesn't match the number of particles")
    c   <- max(x)
    ans <- log(sum(exp(x - c))/nPart) + c
    return(ans)
}

# approximation of log marginal likelihood log p_theta(y_T)
#'@param lwgt A matrix of log-weights with dimension nObs x nPart
#'@param nPart The number of particles
#'@param T Time/Observation time T
#'@return log p_theta(y_T)
logMarginY <- function(lwgt, nPart, T) {
    # check dimension compatibility
    if(ncol(lwgt) != nPart || nrow(lwgt) != T)
        stop("incompatible dimensions!")
    
    logMarginVector <- rep(NA, T)
    logMarginVector <- apply(lwgt, 1, logCondY, nPart = nPart)
    # log p_theta(y_T)
    ans <- sum(logMarginVector)
    return(ans)
}
