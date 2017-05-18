#--- test omp implementation ---------------------------------------------------

devtools::document()
devtools::install()

require(msdeHeaders)

# build model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names,
                       omp = TRUE, showOutput = TRUE, rebuild = TRUE)
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

nReps <- 1
nObs <- 1e7
dT <- 1/252

# initialize
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# generate data
hsim <- sde.sim(model = hmod, init.data = X0, params = Theta,
                dt = dT, dt.sim = dT, N = nObs, nreps = nReps)

ncores <- 1
system.time({
  ll <- replicate(n = 5, expr = {
    sde.loglik(model = hmod, x = hsim$data, dt = dT,
               theta = Theta, ncores = ncores)
  })
})

#--- indexing test -------------------------------------------------------------

Rcpp::sourceCpp(file = "seqTest.cpp")

nobs <- 10
ndims <- 2
par.index <- sample(ndims, nobs, replace = TRUE)
names(par.index) <- 1:nobs-1
par.index[par.index < ndims]
seqTest(par.index, ndims)
