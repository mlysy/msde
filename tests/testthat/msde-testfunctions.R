# maximum absolute and relative difference between arrays
max.diff <- function(x1, x2) {
  c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/max(abs(x1), 1e-8)))
}


# heston model drift and diffusion
hest.dr <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .125 * x[,2]^2, theta[,3]/x[,2] - .5*theta[,2]*x[,2])
}
hest.df <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cv <- .5*theta[,5]*theta[,4]*x[,2]
  ans <- cbind(.25 * x[,2]^2, cv, cv, theta[,4]^2)
  t(apply(ans, 1, function(x) chol(matrix(x,2,2))))
}

# generate heston data/parameters
hest.data <- function(nreps) {
  X0 <- c(X = log(1000), Z = 0.1)
  if(nreps > 1) X0 <- apply(t(replicate(nreps, X0)), 2, jitter)
  X0
}
hest.params <- function(nreps) {
  Theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
  if(nreps > 1) Theta <- apply(t(replicate(nreps, Theta)), 2, jitter)
  Theta
}

#--- R versions of xmvn, zmvn, lmvn --------------------------------------------

xmvn <- function(z, mean, cholsd) {
  n <- length(z)
  cholsd <- matrix(cholsd, n, n)
  c(mean + z %*% cholsd)
}

zmvn <- function(x, mean, cholsd) {
  n <- length(mean)
  cholsd <- matrix(cholsd, n, n)
  backsolve(r = cholsd, x = x - mean, transpose = TRUE)
}

lmvn <- function(x, mean, cholsd) {
  n <- length(mean)
  cholsd <- matrix(cholsd, n, n)
  z <- backsolve(r = cholsd, x = x - mean, transpose = TRUE)
  -.5 * sum(z^2) - sum(log(diag(cholsd)))
}

rmvn <- function(mean, cholsd) {
  n <- length(mean)
  cholsd <- matrix(cholsd, n, n)
  xmvn(rnorm(n), mean, cholsd)
}

# R likelihood
hest.ll <- function(x, theta, dt) {
  ncomp <- nrow(x)
  mu <- x[-ncomp,,drop=FALSE] + hest.dr(x[-ncomp,], theta) * dt
  cholsd <- hest.df(x[-ncomp,,drop=FALSE], theta) * sqrt(dt)
  sum(sapply(2:ncomp-1, function(ii) {
    lmvn(x[ii+1,], mu[ii,], cholsd[ii,])
  }))
}

# R simulation
hest.sim <- function(nobs, dt, rr, x0, theta) {
  X <- matrix(NA, nobs, length(x0))
  x <- x0
  for(ii in 1:nobs) {
    for(jj in 1:rr) {
      mu <- x + hest.dr(x, theta) * (dt/rr)
      csd <- hest.df(x, theta) * sqrt(dt/rr)
      x <- rmvn(mu, csd)
      while(x[2] <= 0) x <- rmvn(mu, csd)
    }
    X[ii,] <- x
  }
  X
}

#--- initialize single/multiple inputs -----------------------------------------

input.init <- function(nreps, sx, st) {
  has.ncomp <- length(nreps) > 1
  if(!has.ncomp) {
    ncomp <- 1
  } else {
    ncomp <- nreps[1]
    nreps <- nreps[2]
  }
  X <- hest.data(nreps*ncomp)
  Theta <- hest.params(nreps)
  if(nreps == 1) {
    Theta <- t(Theta)
  }
  if(nreps*ncomp == 1) {
    X <- t(X)
  }
  if(sx) {
    if(has.ncomp) {
      X <- X[1:ncomp,]
      X.R <- array(X, dim = c(ncomp, ncol(X), nreps))
    } else {
      X.R <- X[rep(1,nreps),,drop=FALSE]
      X <- X[1,]
    }
  } else {
    if(has.ncomp) {
      X <- array(X, dim = c(ncomp, ncol(X), nreps))
    }
    X.R <- X
  }
  if(st) {
    Theta.R <- Theta[rep(1,nreps),,drop=FALSE]
    Theta <- Theta[1,]
  } else Theta.R <- Theta
  list(X = X, X.R = X.R, Theta = Theta, Theta.R = Theta.R)
}
