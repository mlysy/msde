##---------- test file for eou.pf.R --------------------------------------------------------------

#library(msde)
context("exponential OU model -- partilce filter")
source("smc-testfunctions.R")

##---------- eou model settings -------------------------------------------------------------------
# eou model setup
emod <- sde.examples("eou")
data.names <- emod$data.names
param.names <- emod$param.names

# simulate the data
theta0 <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63) # true parameter values
nObs <- 100 # number of observations
nDims <- emod$ndims # number of dimensions
dT <- 1/252 # time between observations (1 year has about 252 trading days)
Y0 <- c(X = rnorm(1), V = rnorm(1)) # initial SDE values
esim <- sde.sim(emod, x0 = Y0, theta = theta0,
                nobs = nObs, # nObs steps forward
                dt = dT, # observation time interval specified by users
                dt.sim = dT/10) # internal observation time
Yt <- esim$data # extract the simulated SDE values (X, V), Yt is an nObs x nDims matrix

# initialize the sde model
nPart <- 50 # number of particles
Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims) # normal draws
# initialization
einit <- sde.init(emod, x = Yt, dt = dT, theta = theta0, nvar.obs = 1, m = 1)

##---------- particle filter in R ------------------------------------------------------------------

pf.R <- pf.fun(einit,
               dr = function(x, theta) sde.drift(emod, x, theta)[1,],
               df = function(x, theta) sde.diff(emod, x, theta)[1,],
               Z = Z)

##---------- particle filter with SMCTC -----------------------------------------------------------

pf <- sde.pf(model = emod, init = einit, npart = nPart, Z = Z)

test_that("pf.R == pf.cpp", {
  expect_equal(pf, pf.R)
})

##---------- testing ------------------------------------------------------------------------------
#----- test smoothing results -----
test_that("Yup.R == Yup.smctc", {
    expect_equal(Yup, Yup2)
})
#----- test log-weights -----
test_that("log-weights.R == log-weights.smctc", {
    expect_equal(lwgt2, tmp2$lwgt)
})
