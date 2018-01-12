# SMC for multivariate SDEs.

# 1. assume level-1 Euler discretization is good enough.
# 2. no filtering, only smoothing.

# will need to create an object from the class template.

#--- quick tests ---------------------------------------------------------------

require(testthat)
require(Rcpp)
require(RcppSMC)

sourceCpp(file = "sdeSMC.cpp")
sourceCpp(file = "sdeSMCtest.cpp")

pf_test() # returns nDims

#--- ok now see about particle updates (no loglik calculation yet...) ----------

context("exponential OU model -- partilce filter")

require(msde)
source("smc-functions.R")

# eou model setup
data.names <- c("X", "V")
param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
emod <- sde.make.model("eouModel.h",
                       data.names = data.names, param.names = param.names)

# simulate some data
theta0 <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63) # true parameter values
nObs <- 100 # number of observations
nDims <- emod$ndims # number of dimensions
dT <- 1/252 # time between observations (1 year has about 252 trading days)
Y0 <- c(X = rnorm(1), V = rnorm(1)) # initial SDE values
esim <- sde.sim(emod, x0 = Y0, theta = theta0,
                nobs = nObs, # nObs steps forward
                dt = dT, # observation time interval specified by users
                dt.sim = dT/10) # internal observation time is dt.sim but only one out of every dt/dt.sim simulation steps is kept in the output
Yt <- esim$data # extract the simulated SDE values (X, V), Yt is an nObs x nDims matrix

# initialize the sde model
# number of particles
nPart <- 50
# normal draws, simulated Brownian motions in small time intervals
Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)
# initialize
# assume all volatilities, i.e. column V in Yt, are missing/unobservable
einit <- sde.init(emod, x = Yt, dt = dT, theta = theta0,
                  nvar.obs = 1, # number of observed variables per timepoint/row in data Yt
                  #nvar.obs = sample(nDims, nObs, replace = TRUE),
                  m = 1) # assume no artificial missing points placed between observations

Yup <- matrix(NA, nObs, nDims*nPart)
lwgt <- matrix(NA, nObs, nPart)
for(ipart in 1:nPart) {
  ind <- (ipart-1)*nDims+(1:nDims) # column index for Xt, Vt corresponding to particle ipart
  tmp <- smc.update(Yt, Z = Z[,ind], # update each particle
                    dt = einit$dt.m, nvar.obs = einit$nvar.obs.m,
                    theta = einit$params,
                    dr = function(x, theta) sde.drift(emod, x, theta)[1,],
                    df = function(x, theta) sde.diff(emod, x, theta)[1,])
  Yup[,ind] <- tmp$X
  lwgt[,ipart] <- tmp$lwgt
}

# ok now in C++: without SMCTC
tmp1 <- pf_update(initParams = einit$params, initData = t(Yt),
                 dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
                 NormalDraws = t(Z))
Yup1 <- t(tmp1$X)
lwgt1 <- t(tmp1$lwgt)

range(Yup - Yup1)
range(lwgt - lwgt1)

# ----- test smoothing results -----
test_that("Yup.R == Yup.cpp", {
    expect_equal(Yup, Yup1)
})
# ----- test log-weights -----
test_that("log-weights.R == log-weights.cpp", {
    expect_equal(lwgt, lwgt1)
})

#--- with SMCTC. -----------------------------------------------------------

tmp2 <- pf_eval(initParams = einit$params, initData = t(Yt),
                dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
                NormalDraws = t(Z))

Yup2 <- t(tmp2$X)
lwgt2 <- t(tmp2$lwgt) # normalized log-weights

# normalize R log-weights
# To normalize the weghts in lwgt to make it comparable with tmp2$lgwt
# the outer apply(.., 1, func...) will implictly change the dimension of lwgt
# since each row it extracts will be treated as a column vector in func(x){...}
lwgtn <- apply(apply(lwgt, 2, cumsum), 1, function(x) {
  mx <- max(x)
  x - (log(sum(exp(x - mx))) + mx) # for avoiding enumerical overflow
})
lwgtn <- t(lwgtn)


range(Yup - Yup2)
range(lwgtn - lwgt2)


# ----- test smoothing results -----
test_that("Yup.R == Yup.smctc", {
    expect_equal(Yup, Yup2)
})
# ----- test log-weights -----
test_that("log-weights.R == log-weights.smctc", {
    expect_equal(lwgt2, tmp2$lwgt)
})

# todo:
# 1. wrap pf_eval as sde.pf(...) in R using msde interface (add to package)
# 2. write a testthat unit test formalizing the tests above.


