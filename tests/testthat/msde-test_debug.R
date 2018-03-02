# script to test generic drift and diffusion functions on cholSD scale.
# calculations are checked against msde C++ code

# inherits: param.names, data.names
# drift.fun, diff.fun
# model (sde.model with default prior)

ndims <- model$ndims
nparams <- model$nparams

source("msde-testfunctions.R")
source("smc-testfunctions.R")

#--- test drift and diffusion --------------------------------------------------

nreps <- 10
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE))
ncases <- nrow(cases)

# drift
test_that("drift.R == drift.cpp", {
  mxd <- matrix(NA, ncases, 2)
  for(ii in 1:ncases) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps, sx, st, randx, randt)
    dr <- sde.drift(model = model, x = init$X, theta = init$Theta)
    dr.R <- drift.fun(x = init$X.R, theta = init$Theta.R)
    if(sx && st) dr.R <- dr.R[1,]
    mxd[ii,] <- max.diff(dr, dr.R)
    expect_equal(mxd[ii,], c(0,0))
  }
})

# diffusion
test_that("diff.R == diff.cpp", {
  mxd <- matrix(NA, ncases, 2)
  for(ii in 1:ncases) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps, sx, st, randx, randt)
    df <- sde.diff(model = model, x = init$X, theta = init$Theta)
    df.R <- diff.fun(x = init$X.R, theta = init$Theta.R)
    if(sx && st) df.R <- df.R[1,]
    mxd[ii,] <- max.diff(df, df.R)
    expect_equal(mxd[ii,], c(0, 0))
  }
})

#--- test simulation -----------------------------------------------------------

SEED <- sample(1000, 1)
dT <- runif(1)
nreps <- 10
nobs <- 8
burn <- 3
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE),
                     burn = c(0, burn), nreps = c(1, nreps), rr = c(1, 2))
ncases <- nrow(cases)

test_that("sim.R == sim.cpp", {
  mxd <- matrix(NA, ncases, 2)
  for(ii in 1:ncases) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    burn <- cases$burn[ii]
    nreps <- cases$nreps[ii]
    rr <- cases$rr[ii]
    init <- input.init(nreps, sx, st, randx, randt)
    set.seed(seed = SEED)
    sim <- sde.sim(model = model, x0 = init$X, theta = init$Theta,
                   dt = dT, dt.sim = dT/rr, nobs = nobs,
                   burn = burn, nreps = nreps, verbose = FALSE)$data
    sim.R <- array(NA, dim = c(nobs, ndims, nreps))
    set.seed(seed = SEED)
    for(jj in 1:nreps) {
      sim.R[,,jj] <- sim.fun(nobs = nobs+burn, dt = dT, rr = rr,
                             x0 = init$X.R[jj,],
                             theta = init$Theta.R[jj,],
                             dr = drift.fun, df = diff.fun,
                             validx = validx)[burn+1:nobs,]
    }
    mxd[ii,] <- max.diff(sim, drop(sim.R))
    expect_equal(mxd[ii,], c(0, 0))
  }
})

#--- test log-likelihood -------------------------------------------------------

cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE))
ncases <- nrow(cases)

test_that("ll.R == ll.cpp", {
  mxd <- matrix(NA, ncases, 2)
  for(ii in 1:ncases) {
    dT <- runif(1)
    nobs <- sample(5:20, 1)
    nreps <- sample(10:20, 1)
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps = c(nobs, nreps), sx, st, randx, randt)
    ll <- sde.loglik(model = model, x = init$X, theta = init$Theta, dt = dT)
    ll.R <- rep(NA, nreps)
    for(jj in 1:nreps) {
      ll.R[jj] <- loglik.fun(x = init$X.R[,,jj], theta = init$Theta.R[jj,],
                             dt = dT, dr = drift.fun, df = diff.fun)
    }
    if(sx && st) {
      ll.R <- ll.R[1]
    }
    mxd[ii,] <- max.diff(ll, ll.R)
    expect_equal(mxd[ii,], c(0, 0), tolerance = 1e-6, scale = 1)
  }
})

#--- test default prior --------------------------------------------------------

nreps <- 10
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE),
                     ntheta = 0:nparams, nx = 0:ndims)
ncases <- nrow(cases)

test_that("lpi.R == lpi.cpp", {
  mxd <- matrix(NA, ncases, 2)
  for(ii in 1:ncases) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps = nreps, sx = sx, st = st, randx ,randt)
    ntheta <- cases$ntheta[ii]
    nx <- cases$nx[ii]
    nrv <- sum(ntheta, nx)
    if(nrv > 0) {
      hnames <- NULL
      if(ntheta > 0) hnames <- c(hnames, sample(model$param.names, ntheta))
      if(nx > 0) hnames <- c(hnames, sample(model$data.names, nx))
      hnames <- sample(hnames)
      mu <- rnorm(nrv)
      names(mu) <- hnames
      Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
      dimnames(Sigma) <- list(hnames, hnames)
      lpi <- sde.prior(model = model, theta = init$Theta, x = init$X,
                       hyper = list(mu = mu, Sigma = Sigma))
      lpi.R <- rep(NA, nreps)
      for(jj in 1:nreps) {
        xx <- c(init$Theta.R[jj,], init$X.R[jj,])
        lpi.R[jj] <- lmvn(x = xx[hnames], mean = mu[hnames],
                          cholsd = chol(Sigma)[hnames,hnames])
      }
    } else {
      lpi <- sde.prior(model = model, theta = init$Theta, x = init$X,
                       hyper = NULL)
      lpi.R <- rep(0, nreps)
    }
    if(sx && st) lpi.R <- lpi.R[1]
    mxd[ii,] <- max.diff(lpi, lpi.R)
    expect_equal(mxd[ii,2], 0)
  }
})

#--- test particle filter ------------------------------------------------------
## Note: after we remove the normal draw input Z from sde.pf,
## we cannot compare pf.R == pf.cpp since two results based on two different Z
## I am not sure if there is a way to reference/point to the Z in sdeFilter constructor
## from R version pf.fun
## I have already kept the successfully tested version (after we fix the discrepancy and before we remove Z)
## as a new branch called "pf-test"

## nreps <- 1
## cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE),
##                     single.history = c(TRUE, FALSE), single.rr = c(5,10))
## ncases <- nrow(cases)

ntest <- 10
test_that("pf.R == pf.cpp", {
  mxd <- matrix(NA, ntest, 4)
  for(ii in 1:ntest) {
    ## sx <- cases$single.x[ii]
    ## st <- cases$single.theta[ii]
    ## history <- cases$single.history[ii]
    ## rr <- cases$single.rr[ii]
    # setup
    nObs <- sample(50:100,1) # number of observations
    nPart <- sample(10:50,1) # number of particles
    nDims <- ndims # number of dimensions
    dT <- runif(1, min = 0, max = 0.2) # time between observations (too large dT will cause testing failure in lotvol model)
    mm <- 1 #sample(1:2, 1)
    history <- as.logical(rbinom(1,1,.5))
    init <- input.init(nreps = 1, sx = TRUE, st = TRUE, randx ,randt)
    msim <- sde.sim(model, x0 = init$X, theta = init$Theta,
                    nobs = nObs, dt = dT, dt.sim = dT)
    # initialization
    # m = 1 implies no missing data time points between two observations
    minit <- sde.init(model, x = msim$data, dt = dT,
                      theta = init$Theta,
                      nvar.obs = sample(nDims, nObs, replace = TRUE), m = mm)
    # normal draws
    Z <- matrix(rnorm(nPart*nDims*(nObs-1)), nObs-1, nPart*nDims)
    # pf in R
    pf.R <- pf.fun(minit, dr = drift.fun, df = diff.fun, Z = Z,
                   history = history)
    # pf in C++ (for debugging, disable the resampling)
    # Z input for sde.pf should be a 3-d array of dimensions (ncomp - 1) x ndims x npart
    Z <- array(c(Z), c(nObs-1, nDims, nPart))
    pf <- sde.pf(model = model, init = minit, npart = nPart,
                resample = "multi", threshold = -1,
                Z = Z, history = history)
    # comparison
    mxd[ii,] <- c(max.diff(pf$data, pf.R$data), max.diff(pf$lwgt, pf.R$lwgt))
    expect_equal(mxd[ii,], rep(0, 4), tolerance = 1e-6, scale = 1)
  }
})
