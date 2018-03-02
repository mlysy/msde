# update particles based on a fixed set of normal draws.
source("msde-testfunctions.R")

# conditional distribution
# returns the normalized log-pdf
cmvn <- function(x2, x1, mean, cholsd) {
  n <- length(mean)
  Sigma <- crossprod(matrix(cholsd, n, n))
  n1 <- length(x1)
  n2 <- length(x2)
  if(n1 == 0) {
    # no conditioning
    mn <- mean
    Vn <- Sigma
  } else {
    m1 <- mean[1:n1]
    m2 <- mean[n1+1:n2]
    V11 <- Sigma[1:n1,1:n1,drop=FALSE]
    V12 <- Sigma[1:n1,n1+1:n2,drop=FALSE]
    V22 <- Sigma[n1+1:n2,n1+1:n2,drop=FALSE]
    # all inversion in one step
    iV <- crossprod(V12, solve(V11, cbind(V12, x1 - m1)))
    mn <- m2 + iV[,n2+1]
    Vn <- V22 - iV[,1:n2]
  }
  # double check the elements selection
  lmvn(x2, mn, chol(Vn))
}

# update
#'@param X The SDE data which is an nObs x nDims matrix
#'@param Z The normal draws, an nObs-1 by nPart*nDims matrix
smc.update <- function(X, Z, dt, nvar.obs, theta, dr, df) {
  nobs <- nrow(X)
  ndims <- ncol(X)
  x <- X[1,]
  if(length(dt) == 1) dt <- rep(dt, nobs-1)
  lw <- rep(NA, nobs)
  lw[1] <- 0
  for(ii in 2:nobs) {
    mu <- x + c(dr(x, theta)) * dt[ii-1]
    csd <- c(df(x, theta)) * sqrt(dt[ii-1])
    z <- zmvn(X[ii,], mu, csd)
    qn <- nvar.obs[ii]
    if(qn < ndims) {
      z[(qn+1):ndims] <- Z[ii-1,(qn+1):ndims]
    }
    x <- xmvn(z, mu, csd)
    X[ii,] <- x
    ## cSd <- matrix(csd, ndims, ndims)
    ## lw[ii] <- lmvn(x[1:qn], mu[1:qn], cSd[1:qn,1:qn,drop=FALSE])
    lw[ii] <- lmvn(x, mu, csd)
    if(qn < ndims) {
      lw[ii] <- lw[ii] - cmvn(x[(qn+1):ndims], x[1:qn], mu, csd)
    }
  }
  list(X = X, lwgt = lw)
}

pf.fun <- function(init, dr, df, Z, history = FALSE) {
  Yt <- init$data
  nObs <- nrow(Yt)
  nDims <- ncol(Yt)
  nPart <- ncol(Z)/nDims
  data <- matrix(NA, nObs, nDims*nPart)
  lwgt <- matrix(NA, nObs, nPart)
  for(ipart in 1:nPart) {
    ind <- (ipart-1)*nDims+(1:nDims) # column index for Xt, Vt corresponding to particle ipart
    tmp <- smc.update(X = Yt, Z = Z[,ind], # update each particle
                      dt = init$dt.m, nvar.obs = init$nvar.obs.m,
                      theta = init$params,
                      dr = dr, df = df)
    data[,ind] <- tmp$X
    lwgt[,ipart] <- tmp$lwgt
  }

  lwgtn <- apply(lwgt, 2, cumsum)
  # return the whole history if history == TRUE
  # or the last observation if history == FALSE
  if (history == FALSE) {
    out <- list(data = data[nObs, ], lwgt = lwgtn[nObs, ])
  } else {
    out <- list(data = data, lwgt = lwgtn)
  }
  return(out)
}
