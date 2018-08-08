#--- quick check for errors ----------------------------------------------------

require(msde)

hmod <- sde.examples(model = "hest")
ndims <- hmod$ndims
nparams <- hmod$nparams
data.names <- hmod$data.names
param.names <- hmod$param.names

# posterior inference
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# simulate data
nObs <- 13
dT <- 1/252
hsim <- sde.sim(model = hmod, x0 = x0, theta = theta,
                dt = dT, dt.sim = dT/100, nobs = nObs, nreps = 1)
# check input parsing
m <- 2 # degree of Euler approximation: dt_euler = dt/m
nvar.obs <- c(2, rep(1, nObs-1)) # all but first vol unobserved
init <- sde.init(model = hmod, x = hsim$data, dt = dT, m = m,
                 nvar.obs = nvar.obs, theta = theta)

# prior
prior <- list(mu = c(.1, .35, 1.0, .5, -.81),
              Sigma = crossprod(matrix(rnorm(25),5)))
prior$Sigma <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$Sigma)
prior$Sigma <- prior$Sigma %*% sqrt(diag(c(.1, 8, .15, .002, .002)))
names(prior$mu) <- param.names
colnames(prior$Sigma) <- param.names
rownames(prior$Sigma) <- param.names
# mcmc specs
#mwg.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
#names(mwg.sd) <- param.names
mwg.sd <- NULL
update.params <- TRUE
update.data <- TRUE
nsamples <- 10 # ifelse(update.data, 2e4, 4e4)
data.out <- list(isamples = sample(nsamples, 1),
                 idims = sample(ndims,1), icomp = sample(nrow(init$data), 5))
burn <- 100
SEED <- 549

hpost <- sde.post(model = hmod,
                  init = init,
                  nsamples = nsamples, burn = burn,
                  hyper = prior,
                  mwg.sd = mwg.sd, adapt = TRUE,
                  data.out = 1:nsamples,
                  update.params = update.params,
                  update.data = update.data, verbose = TRUE)

#--- check posterior -----------------------------------------------------------

# ok want to check inputs from sde.post vs sde.post2

devtools::document()
devtools::install()

require(msde)

# build model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

# posterior inference
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# simulate data
nObs <- 13
dT <- 1/252
hsim <- sde.sim(model = hmod, x0 = x0, theta = theta,
                dt = dT, dt.sim = dT/100, nobs = nObs, nreps = 1)
# check input parsing
m <- 2 # degree of Euler approximation: dt_euler = dt/m
nvar.obs <- c(2, rep(1, nObs-1)) # all but first vol unobserved
init <- sde.init(model = hmod, x = hsim$data, dt = dT, m = m,
                 nvar.obs = nvar.obs, theta = theta)

# prior
prior <- list(mu = c(.1, .35, 1.0, .5, -.81),
              Sigma = crossprod(matrix(rnorm(25),5)))
prior$Sigma <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$Sigma)
prior$Sigma <- prior$Sigma %*% sqrt(diag(c(.1, 8, .15, .002, .002)))
names(prior$mu) <- param.names
colnames(prior$Sigma) <- param.names
rownames(prior$Sigma) <- param.names
# mcmc specs
#mwg.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
#names(mwg.sd) <- param.names
mwg.sd <- NULL
update.params <- TRUE
update.data <- TRUE
nsamples <- 1e6 # ifelse(update.data, 2e4, 4e4)
data.out <- list(isamples = sample(nsamples, 1),
                 idims = sample(ndims,1), icomp = sample(nrow(init$data), 5))
burn <- 100
SEED <- 549

set.seed(SEED)
system.time({
  hpost1 <- sde.post(model = hmod,
                     init = init,
                     nsamples = nsamples, burn = burn,
                     hyper = prior,
                     mwg.sd = mwg.sd, adapt = TRUE,
                     data.out = 1:nsamples,
                     update.params = update.params,
                     update.data = update.data, verbose = TRUE)
})

set.seed(SEED)
hpost2 <- sde.post(model = hmod,
                   init = init,
                   nsamples = nsamples, burn = burn,
                   hyper = prior, debug = FALSE,
                   mwg.sd = mwg.sd, adapt = TRUE,
                   data.out = data.out,
                   update.params = update.params,
                   update.data = update.data)

dout <- hpost2$data.out
all(hpost1$data[dout$icomp,
                dout$idims,
                dout$isamples, drop=FALSE] == hpost2$data)
