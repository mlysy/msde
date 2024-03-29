##---------- test file for eou.pf.R --------------------------------------------------------------

# library(msde)
context("exponential OU model")
# source("smc-testfunctions.R")

# eou model
ModelFile <- "EOU.h"
param.names <- c("alpha", "gamma", "mu", "sigma", "rho")
data.names <- c("X", "V")
model <- sde.make.model(ModelFile = ModelFile,
                        param.names = param.names,
                        data.names = data.names)

# eou model drift and diffusion
drift.fun <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .5 * exp(2*x[,2]), theta[,2] * (theta[,3] - x[,2]))
}

diff.fun <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  df <- matrix(NA, nrow(x), 4)
  df[,1] <- exp(2*x[,2]) # exp(V)
  df[,2] <- theta[,5] * theta[,4] * exp(x[,2]) # rho * sigma * exp(.5*V)
  df[,3] <- df[,2]
  df[,4] <- theta[,4]^2 # sigma^2
  t(apply(df, 1, function(tmp) chol(matrix(tmp,2,2)))) # use sd scale in R
}

# generate eou model data/parameters
randx <- function(nreps) {
  X0 <- c(X = 3.46, V = -2.793)
  # c(X = rnorm(1), V = rnorm(1))
  if(nreps > 1) X0 <- apply(t(replicate(nreps, X0)), 2, jitter)
  X0
}
randt <- function(nreps) {
  Theta <- c(alpha = -0.47, gamma = 5.20, mu = -1.51, sigma = 1.09, rho = -0.9)
  ## c(alpha = .1, gamma = 1, eta = .3, sigma = .2, rho = -.63)
  if(nreps > 1) Theta <- apply(t(replicate(nreps, Theta)), 2, jitter)
  Theta
}

validx <- function(x, theta) return(TRUE)

source("msde-test_debug.R", local = TRUE)
