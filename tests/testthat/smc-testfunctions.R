# update particles based on a fixed set of normal draws.
source("msde-testfunctions.R")

# conditional distribution p(
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

## .cmvn <- function(x, mu, Sigma, qn) {
##   ndims <- length(x)
##   rn <- ndims-qn
##   x <- x[1:qn] # selecting the observed Xn's
##   m2 <- mu[1:qn]
##   m1 <- mu[(qn+1):ndims]
##   # double check the elements selection
##   V22 <- Sigma[1:qn,1:qn,drop=FALSE] # qn x qn
##   V21 <- Sigma[1:qn,(qn+1):ndims,drop=FALSE] # qn x (d - qn)
##   V12 <- Sigma[(qn+1):ndims,1:qn,drop=FALSE] # (d - qn) x qn
##   V11 <- Sigma[(qn+1):ndims,(qn+1):ndims,drop=FALSE] # (d - qn) x (d - qn)
##   # all inversion in one step
##   iV <- crossprod(V21, .solveV(V22, cbind(V21, x - m2)))
##   mn <- m1 + iV[,rn+1]
##   Vn <- V11 - iV[,1:rn]
##   list(m = mn, V = Vn)
## }


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
  # normalize R log-weights
  # the outer apply(.., 1, func...) will implictly change the dimension of lwgt
  # since each row it extracts will be treated as
  # a column vector in func(x){...}
  ## lwgtn <- apply(apply(lwgt, 2, cumsum), 1, function(x) {
  ##   mx <- max(x)
  ##   x - (log(sum(exp(x - mx))) + mx) # for avoiding enumerical overflow
  ## })
  ## lwgtn <- t(lwgtn)
  lwgtn <- apply(lwgt, 2, cumsum)
  # return the whole history if history == TRUE
  # or the last observation if history == FALSE
  if (history == FALSE) {
    out <- list(data = data[nObs, ], lwgt = lwgtn[nObs, ])
  } else {
    out <- list(data = data, lwgt = lwgtn)
  }
  return(out)
  # data <- t(data)
  # data <- as.vector(data)
  # data <- array(data, dim = c(nDims, nPart, nObs))
  # data <- aperm(data, dim = c(2, 1, 3))
  # ans <- list(data = data, lwgt = lwgtn)
  # if(!history) {
  #   ans$data <- ans$data[,,1]
  #   ans$lwgt <- ans$lwgt[,1]
  # }
  # return(ans)
}
