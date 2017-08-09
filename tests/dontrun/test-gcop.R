#--- check the gaussian copula prior -------------------------------------------

require(msde)
require(GaussCop)
source("gcop.hyper.check.R")

# compile model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       PriorFile = "gcopPrior.h",
                       hyper.check = gcop.hyper.check,
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

# ok fit a gcop model
nRV <- ndims + nparams
mu <- rnorm(nRV)
Sigma <- crossprod(matrix(rnorm(nRV^2), nRV, nRV))
X <- t(t(chol(Sigma)) %*% matrix(rnorm(1e5*nRV), nRV, 1e5) + mu)
colnames(X) <- c(param.names, data.names)
gcop <- gcopFit(X = X, fitXD = "normal")

# check parser
tmp <- replicate(n = 100, {
  nrv <- sample(nRV, 1)
  rv.names <- sample(c(param.names, data.names), nrv, replace = FALSE)
  gcs <- gcopSub(gcop, rv.names)
  #gcop.hyper.check(gcs, param.names, data.names)
  # check density
  tx <- rgcop(1, gcop)
  #dgcop(X = rgcop(1, gcs), gcs, log = TRUE)
  c(dgcop(X = tx[,rv.names], gcs, decomp = FALSE, log = TRUE),
    sde.prior(hmod, theta = tx[1:5], x = tx[6:7], hyper = gcs))
})
range(apply(tmp, 2, diff))
