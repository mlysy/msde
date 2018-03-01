#--- memory allocation issue ---------------------------------------------------

require(Rcpp)
require(RcppSMC)
require(msde)
source("smc-functions.R")

sourceCpp(file = "sdeSMC.cpp")

# build sde model
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

# normal draws
nPart <- 4 # number of particles
Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)
einit <- sde.init(emod, x = Yt, dt = dT, theta = theta0,
                  nvar.obs = sample(nDims, nObs, replace = TRUE))
einit$nvar.obs.m

## tmp <- pf_eval_test(initParams = einit$params, initData = t(Yt),
##                dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
##                NormalDraws = t(Z))

tmp <- pf_eval(initParams = einit$params, initData = t(Yt),
               dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
               NormalDraws = t(Z))
