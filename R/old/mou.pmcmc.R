# multivariate OU model based on PMCMC

# require(msde)
source("pmcmc-functions.R")

# bivariate OU model
bmod <- sde.examples("biou")

# simulate some data

# true parameter values
Gamma0 <- .1 * crossprod(matrix(rnorm(4),2,2))
Lambda0 <- rnorm(2)
Phi0 <- crossprod(matrix(rnorm(4),2,2))
Psi0 <- chol(Phi0) # precompiled model uses the Cholesky scale
theta0 <- c(Gamma0, Lambda0, Psi0[c(1,3,4)])
names(theta0) <- bmod$param.names
# initial value
Y0 <- rnorm(2)
names(Y0) <- bmod$data.names

# simulation
dT <- runif(1, max = .1) # time step
nObs <- 10
bsim <- sde.sim(bmod, x0 = Y0, theta = theta0,
                dt = dT, dt.sim = dT, nobs = nObs)
YObs <- bsim$data

# inference via MCMC
binit <- sde.init(bmod, x = YObs, dt = dT, theta = theta0,
                  nvar.obs = 1) # second component is unobserved
# only Lambda1 is unknown
fixed.params <- rep(TRUE, bmod$nparams)
names(fixed.params) <- bmod$param.names
fixed.params["Lambda1"] <- FALSE
# prior on (Lambda1, Y_0)
hyper <- list(mu = c(0,0), Sigma = diag(2))
names(hyper$mu) <- bmod$data.names
dimnames(hyper$Sigma) <- rep(list(bmod$data.names), 2)

# posterior sampling
nsamples <- 1e5
burn <- 1e3
bpost <- sde.post(bmod, binit, hyper = hyper,
                  fixed.params = fixed.params,
                  nsamples = nsamples, burn = burn)
L1.mcmc <- bpost$params[,"Lambda1"]

# analytic posterior
L1.seq <- seq(min(L1.mcmc), max(L1.mcmc), len = 500)
L1.loglik <- sapply(L1.seq, function(l1) {
  lambda <- Lambda0
  lambda[1] <- l1
  mou.loglik(X = YObs, dt = dT, nvar.obs = 1,
             Gamma = Gamma0, Lambda = lambda, Phi = Phi0,
             mu0 = hyper$mu, Sigma0 = hyper$Sigma)
})

# normalize density
L1.Kalman <- exp(L1.loglik - max(L1.loglik))
L1.Kalman <- L1.Kalman/sum(L1.Kalman)/(L1.seq[2]-L1.seq[1])

# compare
hist(L1.mcmc, breaks = 100, freq = FALSE,
     main = expression(p(Lambda[1]*" | "*bold(Y)[1])),
     xlab = expression(Lambda[1]))
lines(L1.seq, L1.Kalman, col = "red")
legend("topright", legend = c("Analytic", "MCMC"),
       pch = c(NA, 22), lty = c(1, NA), col = c("red", "black"))


## Inference via particle MCMC
# iteration = 0, initialize the sde model for PMCMC, using the same binit as MCMC
# get initial log-weights & log marginal likelihood via particle filtering
nPart <- 50
pf <- sde.pf(model = bmod, init = binit, npart = nPart,
            resample = "multi", threshold = 0.5,
            history = TRUE)
lwgt <- pf$lwgt
logYt <- logMarginY(lwgt, nPart, nObs)

# posterior sampling
M <- 1000
nParams <- length(theta0)
acc <- rep(NA, M) # accept or reject indicator vector
thetaMatrix <- matrix(NA, M+1, nParams)
colnames(thetaMatrix) <- bmod$param.names
thetaMatrix[1, ] <- theta0
# start iteration of PMCMC
for(ii in 2:(M+1)) {
    theta_old <- thetaMatrix[ii-1, ]
    # simple random walk proposal for theta
    # only change the unknown (non-fixed) Lambda1
    sd <- rep(NA, nParams)
    sd[1:nParams][fixed.params] <- 0
    sd[1:nParams][!fixed.params] <- 0.2
    theta_prop <- rnorm(nParams, mean = theta_old, sd = sd)
    
    names(theta_prop) <- bmod$param.names
    
    # run the particle filter based on theta_prop
    # assume no artificial missing points placed between observations, m = 1 by default
    tmp_init <- sde.init(bmod, x = YObs, dt = dT, theta = theta_prop, nvar.obs = 1)
    tmp_pf <- sde.pf(model = bmod, init = tmp_init, npart = nPart,
                    resample = "multi", threshold = 0.5, history = TRUE)
    # calculate log p_theta_prop(y_T) 
    lwgt <- tmp_pf$lwgt
    logYt_prop <- logMarginY(lwgt, nPart, nObs)
    
    # calculate the log acceptance ratio
    # remember we use log density 
    # we have assumed prior of theta to be 1 for simplicity
    logRatio <- logYt_prop - logYt + 
        sum(dnorm(theta_old[!fixed.params], mean = theta_prop[!fixed.params], sd = sd[!fixed.params], log = TRUE)) -
        sum(dnorm(theta_prop[!fixed.params], mean = theta_old[!fixed.params], sd = sd[!fixed.params], log = TRUE))
    
    if (logRatio > log(runif(1))) {
        # accept the proposal
        thetaMatrix[ii, ] <- theta_prop
        logYt <- logYt_prop
        acc[ii-1] <- 1
    } else {
        thetaMatrix[ii, ] <- theta_old
        #logYT <- logYT
        acc[ii-1] <- 0
    }
}

# acceptance rate
sum(acc)/M

# posterior mean of theta
L1.pmcmc <- thetaMatrix[ ,!fixed.params]
colMeans(thetaMatrix)
