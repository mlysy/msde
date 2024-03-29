#--- speed/accuracy comparisons with old version of package ----------------

# the old msde package was written without OOP and (*) merged all cpp/h
# files into 1 before compile on-the-fly with Rcpp
# the new version is written with OOP and checks the speed effect of (*)

require(msde) # old version
require(msde2) # new version
pkg.names <- c("msde", "msde2")
max.diff <- function(x1, x2) {
  c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/abs(x1)))
}
pkg1 <- function(fun) {
  eval(parse(text = paste0(pkg.names[1], "::", fun)))
}
pkg2 <- function(fun) {
  eval(parse(text = paste0(pkg.names[2], "::", fun)))
}

# the test model is called heston's model, which is an sde with
# two components y = (x,z), and drift (2x1 vector) and diffusion
# (2x2 upper tri chol matrix) coefficients given by the following
# functions.  for more info check vignette("msde-vignette")
hest.dr <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .125 * Z^2, theta[,3]/Z - .5*theta[,2]*Z)
}
hest.df <- function(X, Z, theta) {
  if(!is.matrix(theta)) theta <- t(theta)
  cv <- .5*theta[,5]*theta[,4]*Z
  ans <- cbind(.25 * Z^2, cv, cv, theta[,4]^2)
  # output as Nx4 matrix, shitty coding to save time
  t(apply(ans, 1, function(x) chol(matrix(x,2,2))))
}

# create model objects (i.e., compile c++ code)

# old version (msde)
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- pkg1("sde.make.model")
hmod <- hmod(ModelFile = "hestModel.h", param.names = param.names,
             data.names = data.names)
## hmod <- pkg1("sde.make.model")(list = hestList,
##                                showOutput = TRUE, rebuild = TRUE)

# new version (msdeHeaders -- interface still in progress)
#ndims <- 2
#nparams <- 5
## param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
## data.names <- c("X", "Z")
hmod2 <- pkg2("sde.make.model")
hmod2 <- hmod2(ModelFile = "hestModel.h",
               param.names = param.names,
               data.names = data.names, OpenMP = FALSE,
               showOutput = TRUE)
ndims <- hmod2$ndims
nparams <- hmod2$nparams


#--- check drift/diffusion in the c++ code ---------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# R code
dr.R <- hest.dr(X = X0[,1], Z = X0[,2], theta = Theta)
df.R <- hest.df(X = X0[,1], Z = X0[,2], theta = Theta)
# msde
dr.p1 <- pkg1("sde.drift")(model = hmod, x = X0, theta = Theta)
df.p1 <- pkg1("sde.diff")(model = hmod, x = X0, theta = Theta)
# msdeHeaders
dr.p2 <- pkg2("sde.drift")(model = hmod2, x = X0, theta = Theta)
df.p2 <- pkg2("sde.diff")(model = hmod2, x = X0, theta = Theta)

# check calcs
tmp <- list(dr.p1, dr.p2)
names(tmp) <- pkg.names
sapply(tmp, function(x) max.diff(dr.R, x))
tmp <- list(df.p1, df.p2)
names(tmp) <- pkg.names
sapply(tmp, function(x) max.diff(df.R, x))

# ok.  speed comparisons
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 1e7
Theta <- matrix(theta, nReps, nparams, byrow = TRUE)
colnames(Theta) <- param.names
X0 <- matrix(x0, nReps, ndims, byrow = TRUE)
colnames(X0) <- data.names

time.p1 <- system.time({
  dr.p1 <- pkg1("sde.drift")(model = hmod, x = X0, theta = Theta)
  df.p1 <- pkg1("sde.diff")(model = hmod, x = X0, theta = Theta)
})
time.p2 <- system.time({
  dr.p2 <- pkg2("sde.drift")(model = hmod2, x = X0, theta = Theta)
  df.p2 <- pkg2("sde.diff")(model = hmod2, x = X0, theta = Theta)
})

# nice boost
tmp <- c(time.p1[3], time.p2[3])
names(tmp) <- pkg.names
tmp/tmp[1]

#--- loglikelihood evaluations ---------------------------------------------

nReps <- 1000
nObs <- 10000
dT <- 1/252

# initialize
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# generate data
## hsim <- pkg1("sde.sim")(model = hmod, init.data = X0, params = Theta,
##                         dt = dT, dt.sim = dT, N = nObs, nreps = nReps)
hsim <- pkg1("sde.sim")
hsim <- hsim(model = hmod, x0 = X0, theta = Theta,
             dt = dT, dt.sim = dT, nobs = nObs, nreps = nReps)

# loglikelihood calcs
time.p1 <- system.time({
  ll.p1 <- pkg1("sde.loglik")(model = hmod, x = hsim$data, dt = dT,
                              theta = Theta)
})
ncores <- 1
time.p2 <- system.time({
  ll.p2 <- pkg2("sde.loglik")(model = hmod2, x = hsim$data,
                              dt = dT, theta = Theta, ncores = ncores)
})

max.diff(ll.p1, ll.p2)

tmp <- c(time.p1[3], time.p2[3])
names(tmp) <- pkg.names
tmp/tmp[1]


#--- forward simulation ----------------------------------------------------

# I've checked correctness of code in test-sde.sim.  it's a bit sloppy but
# it's there.
# if you want to do speed comparisons with exactly the same output
# (otherwise output is stochastic) then change this to true
same.rnd <- TRUE
SEED <- 3843

nReps <- 100
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# each of these takes about 3s
nObs <- 1e5
burn <- 1
dT <- 1/252
dt.sim <- dT/2
if(same.rnd) set.seed(SEED)
time.p1 <- system.time({
  ## hsim.p1 <- pkg1("sde.sim")(model = hmod, init.data = X0, params = Theta,
  ##                            dt = dT, dt.sim = dt.sim,
  ##                            N = nObs+1, burn = burn, nreps = nReps)
  hsim.p1 <- pkg1("sde.sim")(model = hmod, x0 = X0, theta = Theta,
                             dt = dT, dt.sim = dt.sim,
                             nobs = nObs, burn = burn, nreps = nReps)
})
if(same.rnd) set.seed(SEED)
time.p2 <- system.time({
  hsim.p2 <- pkg2("sde.sim")(model = hmod2, x0 = X0, theta = Theta,
                             dt = dT, dt.sim = dt.sim,
                             nobs = nObs, burn = burn, nreps = nReps)
})
## sim1 <- aperm(hsim.p1$data, c(2,3,1))[-1,,]
sim1 <- hsim.p1$data
sim2 <- hsim.p2$data

if(same.rnd) {
  ## sapply(names(hsim.p1),
  ##        function(nm) max.diff(hsim.p1[[nm]], hsim.p2[[nm]]))
  max.diff(sim1, sim2)
}

tmp <- c(time.p1[3], time.p2[3])
names(tmp) <- pkg.names
tmp/tmp[1]


#--- MCMC draws ------------------------------------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

same.rnd <- TRUE
SEED <- 2531

# simulate data
nObs <- 1e3
dT <- 1/252
hsim <- pkg1("sde.sim")(model = hmod, x0 = x0, theta = theta,
                        dt = dT, dt.sim = dT/100, nobs = nObs, nreps = 1)

# initialize MCMC rv's
m <- 2 # degree of Euler approximation: dt_euler = dt/m
nvar.obs <- c(2, rep(1, nObs-1)) # all but first vol unobserved
init1 <- pkg1("sde.init")(hmod, x = hsim$data, dt = dT, m = m,
                          nvar.obs = nvar.obs, theta = theta)
init2 <- pkg2("sde.init")(hmod2, x = hsim$data, dt = dT, m = m,
                          nvar.obs = nvar.obs, theta = theta)

# prior
prior <- list(mu = c(.1, .35, 1.0, .5, -.81),
              Sigma = crossprod(matrix(rnorm(25),5)))
prior$Sigma <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$Sigma)
prior$Sigma <- prior$Sigma %*% sqrt(diag(c(.1, 8, .15, .002, .002)))
names(prior$mu) <- param.names
colnames(prior$Sigma) <- param.names
rownames(prior$Sigma) <- param.names
#prior2 <- list(mu = prior$Mu, Sigma = prior$V)
prior2 <- prior

# mcmc specs
rw.jump.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
names(rw.jump.sd) <- param.names
update.params <- TRUE
update.data <- TRUE
nsamples <- ifelse(update.data, 1e4, 4e4)
burn <- 0

if(same.rnd) set.seed(SEED)
time.p1 <- system.time({
  hpost.p1 <- pkg1("sde.post")(model = hmod, init = init1,
                               nsamples = nsamples, burn = burn,
                               hyper = prior,
                               mwg.sd = rw.jump.sd,
                               update.params = update.params,
                               update.data = update.data)
})

if(same.rnd) set.seed(SEED)
time.p2 <- system.time({
  hpost.p2 <- pkg2("sde.post")(model = hmod2, init = init2,
                               nsamples = nsamples, burn = burn,
                               hyper = prior2,
                               mwg.sd = rw.jump.sd,
                               update.params = update.params,
                               update.data = update.data)
})

if(same.rnd) {
  sapply(names(hpost.p1)[1:2],
         function(nm) max.diff(hpost.p1[[nm]], hpost.p2[[nm]]))
}


# gross
tmp <- c(time.p1[3], time.p2[3])
names(tmp) <- pkg.names
tmp/tmp[1]
