##---------- test file for eou.pf.R --------------------------------------------------------------

# library(msde)
context("exponential OU model")
# source("smc-testfunctions.R")

# eou model
ModelFile <- "eouModel.h"
param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
data.names <- c("X", "V")
model <- sde.make.model(ModelFile = ModelFile,
                        param.names = param.names,
                        data.names = data.names)

# eou model drift and diffusion
drift.fun <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .5 * exp(x[,2]), -(theta[,2] * x[,2] + theta[,3]))
}

diff.fun <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta) 
  df <- matrix(NA, nrow(x), 4)
  df[,1] <- exp(x[,2]) # exp(V)
  df[,2] <- theta[,5] * theta[,4] * exp(.5 * x[,2]) # rho * sigma * exp(.5*V)
  df[,3] <- df[,2]
  df[,4] <- theta[,4]^2 # sigma^2
  t(apply(df, 1, function(tmp) chol(matrix(tmp,2,2)))) # use sd scale in R
}

# generate eou model data/parameters
randx <- function(nreps) {
  X0 <- c(X = rnorm(1), V = rnorm(1))
  if(nreps > 1) X0 <- apply(t(replicate(nreps, X0)), 2, jitter)
  X0
}
randt <- function(nreps) {
  Theta <- c(alpha = .1, gamma = 1, eta = .3, sigma = .2, rho = -.63)
  if(nreps > 1) Theta <- apply(t(replicate(nreps, Theta)), 2, jitter)
  Theta
}

validx <- function(x, theta) return(TRUE)

source("msde-test_debug.R", local = TRUE)


# ## ------- Old testing steps below ---------
# # simulate the data
# theta0 <- c(alpha = .1, gamma = 4.8, eta = 0.1, sigma = .1, rho = -.63) # true parameter values
# nObs <- 100 # number of observations
# nDims <- model$ndims # number of dimensions
# dT <- 1/252 # time between observations (1 year has about 252 trading days)
# Y0 <- c(X = rnorm(1), V = rnorm(1)) # initial SDE values
# esim <- sde.sim(model, x0 = Y0, theta = theta0,
#                 nobs = nObs, # nObs steps forward
#                 dt = dT, # observation time interval specified by users
#                 dt.sim = dT/10) # internal observation time
# Yt <- esim$data # extract the simulated SDE values (X, V), Yt is an nObs x nDims matrix

# # initialize the sde model
# nPart <- 50 # number of particles
# Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims) # normal draws
# # initialization
# einit <- sde.init(model, x = Yt, dt = dT, theta = theta0, nvar.obs = 1, m = 1)

# ##---------- particle filter in R ------------------------------------------------------------------

# pf.R <- pf.fun(einit,
#                dr = function(x, theta) sde.drift(model, x, theta)[1,],
#                df = function(x, theta) sde.diff(model, x, theta)[1,],
#                Z = Z)

# ##---------- particle filter with SMCTC -----------------------------------------------------------

# pf <- sde.pf(model = model, init = einit, npart = nPart, Z = Z)

# test_that("pf.R == pf.cpp", {
#   expect_equal(pf, pf.R)
# })

# ##---------- testing ------------------------------------------------------------------------------
# #----- test smoothing results -----
# test_that("Yup.R == Yup.smctc", {
#     expect_equal(Yup, Yup2)
# })
# #----- test log-weights -----
# test_that("log-weights.R == log-weights.smctc", {
#     expect_equal(lwgt2, tmp2$lwgt)
# })
