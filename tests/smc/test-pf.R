# SMC for multivariate SDEs.

# 1. assume level-1 Euler discretization is good enough.
# 2. no filtering, only smoothing.

# will need to create an object from the class template.

#--- quick tests ---------------------------------------------------------------

require(Rcpp)
require(RcppSMC)

sourceCpp(file = "sdeSMC.cpp")

pf_test() # returns nDims

#--- ok now see about particle updates (no loglik calculation yet...) ----------

require(msde)
source("smc-functions.R")

data.names <- c("X", "A")
param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
emod <- sde.make.model("eouModel.h",
                       data.names = data.names, param.names = param.names)

# simulate some data
theta0 <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63)
nObs <- 10
nDims <- emod$ndims
dT <- runif(1)
Y0 <- c(X = rnorm(1), A = rnorm(1))

esim <- sde.sim(emod, x0 = Y0, theta = theta0,
                dt = dT, dt.sim = dT/10, nobs = nObs)
Yt <- esim$data

# ok infill

# normal draws
nPart <- 50 # number of particles
Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)

einit <- sde.init(emod, x = Yt, dt = dT, theta = theta0,
                  nvar.obs = sample(nDims, nObs, replace = TRUE))

Yup <- matrix(NA, nObs, nDims*nPart)
lwgt <- matrix(NA, nObs, nPart)
for(ipart in 1:nPart) {
  ind <- (ipart-1)*nDims+(1:nDims)
  tmp <- smc.update(Yt, Z = Z[,ind],
                    dt = einit$dt.m, nvar.obs = einit$nvar.obs.m,
                    theta = einit$params,
                    dr = function(x, theta) sde.drift(emod, x, theta)[1,],
                    df = function(x, theta) sde.diff(emod, x, theta)[1,])
  Yup[,ind] <- tmp$X
  lwgt[,ipart] <- tmp$lwgt
}

# ok now in C++: without SMCTC
tmp <- pf_update(initParams = einit$params, initData = t(Yt),
                 dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
                 NormalDraws = t(Z))
Yup2 <- t(tmp$X)
lwgt2 <- t(tmp$lwgt)

range(Yup - Yup2)
range(lwgt - lwgt2)

# with SMCTC.
tmp2 <- pf_eval(initParams = einit$params, initData = t(Yt),
                dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
                NormalDraws = t(Z))

Yup2 <- t(tmp2$X)
lwgt2 <- apply(apply(lwgt, 2, cumsum), 1, function(x) {
  mx <- max(x)
  x - (log(sum(exp(x - mx))) + mx)
})

range(Yup - t(tmp2$X))
range(lwgt2 - tmp2$lwgt)

# todo:
# 1. wrap pf_eval as sde.pf(...) in R using msde interface (add to package)
# 2. write a testthat unit test formalizing the tests above.
