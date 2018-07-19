# particle mcmc for sde model
# currently only support particle marginal metropolis hastings sampler
# TODO: add more pmcmc samplers, complete unit tests

# context("exponential OU model -- partilce mcmc")

require(Rcpp)
require(msde)

source("pmcmc-functions.R")
sourceCpp(file = "sdeSMC.cpp")



# eou model setup
data.names <- c("X", "V")
param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
emod <- sde.make.model("eouModel.h",
                       data.names = data.names, param.names = param.names)

# simulate synthetic data
theta <- c(alpha = .5, gamma = 3.14159, eta = .1, sigma = .32, rho = -.63) # used as true theta
nObs  <- 100 # number of observations
nDims <- emod$ndims # number of dimensions 
dT    <- 1/252 # time between observations (1 year has about 252 trading days)
Y0    <- c(X = rnorm(1), V = rnorm(1)) # initial SDE values
esim  <- sde.sim(emod, x0 = Y0, theta = theta,
                nobs = nObs, # nObs steps forward
                dt = dT, # observation time interval specified by users
                dt.sim = dT/10) # internal observation time is dt.sim but only one out of every dt/dt.sim simulation steps is kept in the output
Yt <- esim$data # extract the simulated SDE values (X, V), Yt is an nObs x nDims matrix

# iteration = 0
# initialize the sde model for PMCMC
# set initial parameters arbitrarily (but should be valid)
theta0 <- c(alpha = 1, gamma = 5, eta = 0.1, sigma = 1.6, rho = 0)
# number of particles
nPart <- 50
# normal draws, simulated Brownian motions in small time intervals
Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)
# assume all volatilities, i.e. column V in Yt, are missing/unobservable
einit <- sde.init(emod, x = Yt, dt = dT, theta = theta0,
                  nvar.obs = 1, # number of observed variables per timepoint/row in data Yt
                  m = 1) # assume no artificial missing points placed between observations
# run the initial SMC in C++: without SMCTC
tmp0 <- pf_update(initParams = einit$params, initData = t(Yt),
                  dT = einit$dt.m, nDimsPerObs = einit$nvar.obs.m,
                  NormalDraws = t(Z))



lwgt <- t(tmp0$lwgt)
logYT <- logMarginT(lwgt, nPart, nObs)

# total number of iterations
M <- 1000
nParams <- length(theta0)
acc <- rep(NA, M) # accept or reject indicator vector
thetaMatrix <- matrix(NA, M+1, nParams)
colnames(thetaMatrix) <- emod$param.names
thetaMatrix[1, ] <- theta0
# start iteration of PMCMC
for(ii in 2:(M+1)) {
    theta_old <- thetaMatrix[ii-1, ]
    # simple random walk proposal for theta
    # assume independence between different parameters
    theta_prop <- rnorm(nParams, mean = theta_old, sd = c(.2, .05, .1, .01, .03))
    
    names(theta_prop) <- emod$param.names
    
    # run an SMC based on theta_prop
    tmp_einit <- sde.init(emod, x = Yt, dt = dT, theta = theta_prop,
                      nvar.obs = 1, # number of observed variables per timepoint/row in data Yt
                      m = 1) # assume no artificial missing points placed between observations
    tmp_smc <- pf_update(initParams = tmp_einit$params, initData = t(Yt),
                         dT = tmp_einit$dt.m, nDimsPerObs = tmp_einit$nvar.obs.m,
                         NormalDraws = t(Z))
    # calculate log p_theta_prop(y_T) 
    lwgt <- t(tmp_smc$lwgt)
    logYT_prop <- logMarginT(lwgt, nPart, nObs)
    
    # calculate the log acceptance ratio
    # remember we use log density 
    # we have assumed prior of theta to be 1 for simplicity
    logRatio <- logYT_prop - logYT + 
        sum(dnorm(theta_old, mean = theta_prop, sd = rep(.1, nParams), log = TRUE)) -
        sum(dnorm(theta_prop, mean = theta_old, sd = rep(.1, nParams), log = TRUE))
    
    if (logRatio > log(runif(1))) {
        # accept the proposal
        thetaMatrix[ii, ] <- theta_prop
        logYT <- logYT_prop
        acc[ii-1] <- 1
        # we don't care the update of X
        # we just want to estimate the parameter
    } else {
        thetaMatrix[ii, ] <- theta_old
        #logYT <- logYT
        acc[ii-1] <- 0
    }
}

# acceptance rate
sum(acc)/M

# posterior mean of theta

colMeans(thetaMatrix)

gamma_post <- thetaMatrix

# Metropolis-within-Gibbs sampler given by msde
nsamples <- 1e3
burn <- 2e2
hyper <- NULL # flat prior
hest.post <- sde.post(model = emod, init = einit, hyper = hyper,
                      nsamples = nsamples, burn = burn)



# to avoid invalid parameter proposal
# there is unsolved issue
if(theta_prop[2] <= 0 || theta_prop[4] <= 0 || theta_prop[5] < -1 || theta_prop > -1) {
    # directly reject the proposal
    thetaMatrix[ii, ] <- theta_old
    #logYT <- logYT
    acc[ii-1] <- 0
    next
}
