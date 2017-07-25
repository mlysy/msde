#--- test omp implementation ---------------------------------------------------

# check that omp = T/F models can be simultaneously loaded
# and give same results


devtools::document()
devtools::install()

require(msde)

# build model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hfile <- sde.examples("hest", file.only = TRUE)

# with omp
hmod <- sde.make.model(ModelFile = hfile,
                       param.names = param.names,
                       data.names = data.names,
                       OpenMP = TRUE)

# without omp
hmod2 <- sde.make.model(ModelFile = hfile,
                        param.names = param.names,
                        data.names = data.names,
                        OpenMP = FALSE)


# NOTE: can't compile both versions in the same R session,
# since Rcpp doesn't create a different .so, so get a segmentation fault.
## # without omp
## hmod2 <- sde.make.model(ModelFile = "hestModel.h",
##                         param.names = param.names,
##                         data.names = data.names,
##                         OpenMP = FALSE, showOutput = TRUE, rebuild = TRUE)
ndims <- hmod$ndims
nparams <- hmod$nparams

#--- drift test ----------------------------------------------------------------

# ncores for sde.drift now DISABLED.
# best to parallelize this from within R

# generate data
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
nReps <- 1e7
Theta <- matrix(theta, nReps, nparams, byrow = TRUE)
colnames(Theta) <- param.names
X0 <- matrix(x0, nReps, ndims, byrow = TRUE)
colnames(X0) <- data.names

## ncores <- 8 # set to < 0 or NA to disable omp
## system.time({
##   dr <- sde.drift(model = hmod, x = X0, theta = Theta,
##                   ncores = ncores, debug = FALSE)
## })

#--- loglik test ---------------------------------------------------------------

# parallel is better when nObs is very large and nReps is relatively small
# serial is better when nObs is small, presumably because memory allocation
# is bottleneck

nReps <- 100
nObs <- 100
dT <- 1/252

# initialize
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# generate data
hsim <- sde.sim(model = hmod, x0 = X0, theta = Theta,
                dt = dT, dt.sim = dT, nobs = nObs, nreps = nReps)


# multicore
ncores <- sample(1:8,1)
system.time({
  ll1 <- sde.loglik(model = hmod, x = hsim$data, dt = dT,
                    theta = Theta, ncores = ncores)
})

# single core
system.time({
  sde.loglik(model = hmod2, x = hsim$data, dt = dT,
             theta = Theta)
})

range(ll1-ll2)

#--- MCMC test -----------------------------------------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# simulate data
nObs <- 1e3
dT <- 1/252
hsim <- sde.sim(model = hmod, x0 = x0, theta = theta,
                dt = dT, dt.sim = dT/100, nobs = nObs, nreps = 1)

# prior
hyper <- list(mu = c(.1, .35, 1.0, .5, -.81),
              Sigma = crossprod(matrix(rnorm(25),5)))
hyper$Sigma <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(hyper$Sigma)
hyper$Sigma <- hyper$Sigma %*% sqrt(diag(c(.1, 8, .15, .002, .002)))
names(hyper$mu) <- param.names
colnames(hyper$Sigma) <- param.names
rownames(hyper$Sigma) <- param.names


# initialize
init <- sde.init(model = hmod, x = hsim$data, dt = dT, m = 2,
                 nvar.obs = 1, theta = theta)

update.params <- TRUE
update.data <- TRUE
nsamples <- 1e3 # ifelse(update.data, 2e4, 4e4)
burn <- 0
SEED <- sample(999,1)

# multicore
ncores <- 8
set.seed(SEED)
system.time({
  hpost <- sde.post(model = hmod, init = init,
                    nsamples = nsamples, burn = burn,
                    hyper = hyper, ncores = ncores,
                    update.params = update.params,
                    update.data = update.data)
})

# single core
set.seed(SEED)
system.time({
  hpost2 <- sde.post(model = hmod2, init = init,
                     nsamples = nsamples, burn = burn,
                     hyper = hyper,
                     update.params = update.params,
                     update.data = update.data)
})

# ok great
identical(hpost, hpost2)

#--- indexing test -------------------------------------------------------------

Rcpp::sourceCpp(file = "seqTest.cpp")

nobs <- 10
ndims <- 2
par.index <- sample(ndims, nobs, replace = TRUE)
names(par.index) <- 1:nobs-1
par.index[par.index < ndims]
seqTest(par.index, ndims)

#--- further parallel testing --------------------------------------------------

Sys.setenv(PKG_CXXFLAGS = "-fopenmp")
Sys.setenv(PKG_LIBS = "-fopenmp")

Rcpp::sourceCpp(file = "ompTest.cpp", showOutput = TRUE, rebuild = TRUE)

ompTest(nCores = 2)
