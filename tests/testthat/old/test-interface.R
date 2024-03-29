#--- test msdeHeaders interface ------------------------------------------------

library(msdeHeaders)
source("msde-testfunctions.R")
context("msde interface")

# compile heston model
norho <- FALSE
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
if(norho) param.names <- param.names[1:4]
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = ifelse(norho, "hestModel_diag.h",
                                          "hestModel_var.h"),
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

#--- test drift and diffusion --------------------------------------------------

nreps <- 10
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE))

# drift
test_that("drift.R == drift.cpp", {
  mxd <- sapply(1:nrow(cases), function(ii) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps, sx, st, norho)
    dr <- sde.drift(model = hmod, x = init$X, theta = init$Theta)
    dr.R <- hest.dr(x = init$X.R, theta = init$Theta.R)
    if(sx && st) dr.R <- dr.R[1,]
    max.diff(dr, dr.R)
  })
  sapply(mxd, expect_equal, expected = 0)
})

# diffusion
test_that("diff.R == diff.cpp", {
  mxd <- sapply(1:nrow(cases), function(ii) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps, sx, st, norho)
    df <- sde.diff(model = hmod, x = init$X, theta = init$Theta)
    df.R <- hest.df(x = init$X.R, theta = init$Theta.R, norho = norho)
    if(sx && st) df.R <- df.R[1,]
    max.diff(df, df.R)
  })
  sapply(mxd, expect_equal, expected = 0)
})

#--- test simulation -----------------------------------------------------------

SEED <- sample(1000, 1)
dT <- runif(1)
nreps <- 10
nobs <- 8
burn <- 3
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE),
                     burn = c(0, burn), nreps = c(1, nreps), rr = c(1, 2))

test_that("sim.R == sim.cpp", {
  mxd <- sapply(1:nrow(cases), function(ii) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    burn <- cases$burn[ii]
    nreps <- cases$nreps[ii]
    rr <- cases$rr[ii]
    init <- input.init(nreps, sx, st, norho)
    set.seed(seed = SEED)
    sim <- sde.sim(model = hmod, x0 = init$X, theta = init$Theta,
                   dt = dT, dt.sim = dT/rr, nobs = nobs,
                   burn = burn, nreps = nreps, verbose = FALSE)$data
    sim.R <- array(NA, dim = c(nobs, ndims, nreps))
    set.seed(seed = SEED)
    for(jj in 1:nreps) {
      sim.R[,,jj] <- hest.sim(nobs = nobs+burn, dt = dT, rr = rr,
                              x0 = init$X.R[jj,], norho = norho,
                              theta = init$Theta.R[jj,])[burn+1:nobs,]
    }
    max.diff(sim, sim.R)
  })
  sapply(mxd, expect_equal, expected = 0)
})

#--- test log-likelihood -------------------------------------------------------

dT <- runif(1)
nreps <- 10
nobs <- 8
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE))

test_that("ll.R == ll.cpp", {
  mxd <- sapply(1:nrow(cases), function(ii) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps = c(nobs, nreps), sx = sx, st = st, norho)
    ll <- sde.loglik(model = hmod, x = init$X, theta = init$Theta, dt = dT)
    ll.R <- rep(NA, nreps)
    for(jj in 1:nreps) {
      ll.R[jj] <- hest.ll(x = init$X.R[,,jj], theta = init$Theta.R[jj,],
                          dt = dT, norho = norho)
    }
    if(sx && st) {
      ll.R <- ll.R[1]
    }
    max.diff(ll, ll.R)
  })
  sapply(mxd, expect_equal, expected = 0)
})

#--- test default prior --------------------------------------------------------

nreps <- 10
cases <- expand.grid(single.x = c(TRUE, FALSE), single.theta = c(TRUE, FALSE),
                     ntheta = 0:nparams, nx = 0:ndims)

test_that("lpi.R == lpi.cpp", {
  mxd <- sapply(1:nrow(cases), function(ii) {
    sx <- cases$single.x[ii]
    st <- cases$single.theta[ii]
    init <- input.init(nreps = nreps, sx = sx, st = st, norho)
    ntheta <- cases$ntheta[ii]
    nx <- cases$nx[ii]
    nrv <- sum(ntheta, nx)
    if(nrv > 0) {
      hnames <- NULL
      if(ntheta > 0) hnames <- c(hnames, sample(hmod$param.names, ntheta))
      if(nx > 0) hnames <- c(hnames, sample(hmod$data.names, nx))
      hnames <- sample(hnames)
      mu <- rnorm(nrv)
      names(mu) <- hnames
      Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
      dimnames(Sigma) <- list(hnames, hnames)
      lpi <- sde.prior(model = hmod, theta = init$Theta, x = init$X,
                       hyper = list(mu = mu, Sigma = Sigma))
      lpi.R <- rep(NA, nreps)
      for(jj in 1:nreps) {
        xx <- c(init$Theta.R[jj,], init$X.R[jj,])
        lpi.R[jj] <- lmvn(x = xx[hnames], mean = mu[hnames],
                          cholsd = chol(Sigma)[hnames,hnames])
      }
    } else {
      lpi <- sde.prior(model = hmod, theta = init$Theta, x = init$X,
                       hyper = NULL)
      lpi.R <- rep(0, nreps)
    }
    if(sx && st) lpi.R <- lpi.R[1]
    max.diff(lpi, lpi.R)
  })
  sapply(mxd[,2], expect_equal, expected = 0)
})
