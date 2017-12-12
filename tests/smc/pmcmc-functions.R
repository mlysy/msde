# pmcmc helper functions

# log approximation of marginal likelihood p_theta(y_n | y_n-1)
#'@param x Vector of log-weights at a given time
#'@param nPart Number of particles
#'@return log marginal approximated density log p_theta(y_n | y_n-1)
logMargin <- function(x, nPart) {
    # check dimension compatibility
    if(length(x) != nPart)
        stop("Check the input vector, its length doesn't match the number of particles")
    c   <- max(x)
    ans <- log(sum(exp(x - c))/nPart) + c
    return(ans)
}

#'@param lwgt A matrix of log-weights with dimension nObs x nPart
#'@param nPart The number of particles
#'@param nObs The number of observation time T
#'@return log p_theta(y_T)
logMarginT <- function(lwgt, nPart, nObs) {
    # check dimension compatibility
    if(ncol(lwgt) != nPart || nrow(lwgt) != nObs)
        stop("incompatible dimensions!")
    
    logMarginVector <- rep(NA, nObs)
    logMarginVector <- apply(lwgt, 1, logMargin, nPart = nPart)
    # log p_theta(y_T)
    ans <- sum(logMarginVector)
    return(ans)
}

# inverse-gamma density.
# for improper density does not compute normalizing constant.
dinvgam <- function(x, shape, scale, log = FALSE) {
    # unnormalized density
    lp <- -(shape + 1) * log(x) - scale/x
    # normalizing constant
    lp <- lp + ifelse(shape > 0 & scale > 0,
                      shape * log(scale) - lgamma(shape), 0)
    if(!log) lp <- exp(lp)
    return(lp)
}

# sample from inverse-gamma distribution
rinvgam <- function(n, shape, scale) {
    1/rgamma(n, shape = shape, rate = scale)
}

# update theta proposal
# to create more appropriate proposal for the parameter
# not finished yet
update_theta <- function(theta) {
    alpha <- theta[1]
    gamma <- theta[2]
    eta <- theta[3]
    sigma <- theta[4]
    rho <- theta[5]
    
    alpha_prop <- rnorm(1, mean = alpha, sd = .1)
    gamma_prop <- rinvgam(1, 0.1, 0.1)
    eta_prop <- rnorm(1, mean = eta, sd = .2)
}