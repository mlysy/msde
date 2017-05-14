#--- check posterior -----------------------------------------------------------

# ok want to check inputs from sde.post vs sde.post2

devtools::document()
devtools::install()

require(msdeHeaders)

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
same.rnd <- TRUE
SEED <- 2531

# simulate data
nObs <- 1e3
dT <- 1/252
hsim <- sde.sim(model = hmod, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)
# check input parsing
m <- 2 # degree of Euler approximation: dt_euler = dt/m
par.index <- c(2, rep(nObs-1, 1)) # all but first vol unobserved
init <- sde.init(data = hsim$data, dt = dT, m = m,
                 par.index = par.index)
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
nsamples <- 1e3 # ifelse(update.data, 2e4, 4e4)
burn <- 1e2

if(same.rnd) set.seed(SEED)
hpost1 <- sde.post(model = hmod,
                   init.data = init, init.params = theta,
                   nsamples = nsamples, burn = burn,
                   hyper.params = prior, debug = FALSE,
                   mwg.sd = mwg.sd, adapt = TRUE,
                   update.params = update.params,
                   update.data = update.data)

if(same.rnd) set.seed(SEED)
hpost2 <- sde.post2(model = hmod,
                    init = init, init.params = theta,
                    nsamples = nsamples, burn = burn,
                    prior = prior, debug = FALSE,
                    rw.jump.sd = mwg.sd,
                    update.data = update.data, update.params = update.params)

sapply(names(hpost1), function(ii) identical(hpost1[[ii]], hpost2[[ii]]))
