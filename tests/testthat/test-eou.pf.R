##---------- test file for eou.pf.R --------------------------------------------------------------

#library(msde)
context("exponential OU model -- partilce filter")
source("smc-functions.R")

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

# normalize R log-weights
# To normalize the weghts in lwgt to make it comparable with tmp2$lgwt
# the outer apply(.., 1, func...) will implictly change the dimension of lwgt
# since each row it extracts will be treated as a column vector in func(x){...}
lwgtn <- apply(apply(lwgt, 2, cumsum), 1, function(x) {
  mx <- max(x)
  x - (log(sum(exp(x - mx))) + mx) # for avoiding enumerical overflow
})
lwgtn <- t(lwgtn)

##---------- particle filter with SMCTC -----------------------------------------------------------
tmp2 <- sde.pf(model = emod, init = einit, npart = nPart, Z = Z)
Yup2 <- tmp2$X
lwgt2 <- tmp2$lwgt # normalized log-weights

##---------- testing ------------------------------------------------------------------------------
#----- test smoothing results -----
test_that("Yup.R == Yup.smctc", {
    expect_equal(Yup, Yup2)
})
#----- test log-weights -----
test_that("log-weights.R == log-weights.smctc", {
    expect_equal(lwgt2, tmp2$lwgt)
})
