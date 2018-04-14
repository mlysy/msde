##---------- test file for EOU NA issue --------------------------------------------------------------
# test the EOU model provided by Ethan
# d xt = (alpha - exp(2*vt)/2) dt + exp(vt) d Bxt
# d vt = gamma*(mu - vt) dt + sigma d Bzt
# d Bxt d Bzt = rho dt

# library(msde)
context("exponential OU model - NA issue")
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
  # df <- matrix(NA, nrow(x), 4)
  # df[,1] <- exp(2*x[,2]) # exp(2V)
  # df[,2] <- theta[,5] * theta[,4] * exp(x[,2]) # rho * sigma * exp(V)
  # df[,3] <- df[,2]
  # df[,4] <- theta[,4]^2 # sigma^2
  # t(apply(df, 1, function(tmp) chol(matrix(tmp,2,2)))) # use sd scale in R
  df <- matrix(0, nrow(x), 4)
  df[,1] <- exp(x[,2])
  df[,3] <- theta[,4]
  df[,4] <- sqrt(1.0 - theta[,5]*theta[,5])*df[,3]
  df[,3] <- df[,3] * theta[,5]
  return(df)
}

# generate eou model data/parameters
randx <- function(nreps) {
  X0 <- c(X = 3, V = -2)
  #if(nreps > 1) X0 <- apply(t(replicate(nreps, X0)), 2, jitter)
  if(nreps > 1) X0 <- t(replicate(nreps, X0)) # used for no jitter
  return(X0)
}
randt <- function(nreps) {
  Theta <- c(alpha = 0.1, gamma = 5, mu = -2.6, sigma = 1.2, rho = -0.6)
  #if(nreps > 1) Theta <- apply(t(replicate(nreps, Theta)), 2, jitter)
  if(nreps > 1) Theta <- t(replicate(nreps, Theta)) # used for no jitter
  return(Theta)
}

validx <- function(x, theta) return(TRUE)

source("msde-test_debug.R", local = TRUE)
