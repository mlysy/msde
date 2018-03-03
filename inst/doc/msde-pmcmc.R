## ----seed----------------------------------------------------------------
set.seed(123)

## ----bmod----------------------------------------------------------------
library(msde)
bmod <- sde.examples("biou")

## ----init-values---------------------------------------------------------
# parameter values
Gamma0 <- .1 * crossprod(matrix(rnorm(4),2,2))
Lambda0 <- rnorm(2)
Phi0 <- crossprod(matrix(rnorm(4),2,2))
Psi0 <- chol(Phi0) # precompiled model uses the Cholesky scale
theta0 <- c(Gamma0, Lambda0, Psi0[c(1,3,4)])
names(theta0) <- bmod$param.names
# initial sde value
Y0 <- rnorm(2)
names(Y0) <- bmod$data.names

## ----sim-----------------------------------------------------------------
# simulation
dT <- runif(1, max = .1) # time step
nObs <- 10
bsim <- sde.sim(bmod, x0 = Y0, theta = theta0,
                dt = dT, dt.sim = dT, nobs = nObs)
YObs <- bsim$data

## ----initialization------------------------------------------------------
# initialization before MCMC
binit <- sde.init(bmod, x = YObs, dt = dT, theta = theta0,
                  nvar.obs = 1) # second component is unobserved

## ----fixed-params--------------------------------------------------------
# only Lambda1 is unknown
fixed.params <- rep(TRUE, bmod$nparams)
names(fixed.params) <- bmod$param.names
fixed.params["Lambda1"] <- FALSE

## ----prior---------------------------------------------------------------
# prior on (Lambda1, Y_0)
hyper <- list(mu = c(0,0), Sigma = diag(2))
names(hyper$mu) <- bmod$data.names
dimnames(hyper$Sigma) <- rep(list(bmod$data.names), 2)

## ----L1-mcmc-------------------------------------------------------------
L1.mcmc <- bpost$params[,"Lambda1"]

## ----mcmc-plot, fig.width = 10, fig.height = 5, out.width = "90%"--------
# compare MCMC with Kalman filter
hist(L1.mcmc, breaks = 100, freq = FALSE,
     main = expression(p(Lambda[1]*" | "*bold(Y)[1])),
     xlab = expression(Lambda[1]))
lines(L1.seq, L1.Kalman, col = "red")
legend("topright", legend = c("Analytic", "MCMC"),
       pch = c(NA, 22), lty = c(1, NA), col = c("red", "black"))

## ----accept--------------------------------------------------------------
# check the acceptance rate
accept <- ppost$accept
print(accept)

## ------------------------------------------------------------------------
# posterior mean of theta
L1.pmcmc <- ppost$params[ ,!fixed.params]

## ----pmcmc-plot, fig.width = 10, fig.height = 5, out.width = "90%"-------
# compare particle MCMC with Kalman filter
hist(L1.pmcmc, breaks = 100, freq = FALSE,
     main = expression(p(Lambda[1]*" | "*bold(Y)[1])),
     xlab = expression(Lambda[1]))
lines(L1.seq, L1.Kalman, col = "red")
legend("topright", legend = c("Analytic", "PMCMC"),
       pch = c(NA, 22), lty = c(1, NA), col = c("red", "black"))

