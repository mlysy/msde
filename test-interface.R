#--- msde interface ------------------------------------------------------------

require(msdeHeaders)
# TODO: shouldn't need to rebuild = TRUE
# TODO: prior specification
max.diff <- function(x1, x2) {
  c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/abs(x1)))
}

#--- example 1: heston model ---------------------------------------------------

param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hest.dr <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cbind(theta[,1] - .125 * x[,2]^2, theta[,3]/x[,2] - .5*theta[,2]*x[,2])
}
hest.df <- function(x, theta) {
  if(!is.matrix(x)) x <- t(x)
  if(!is.matrix(theta)) theta <- t(theta)
  cv <- .5*theta[,5]*theta[,4]*x[,2]
  ans <- cbind(.25 * x[,2]^2, cv, cv, theta[,4]^2)
  t(apply(ans, 1, function(x) chol(matrix(x,2,2))))
}

# compile model
# TODO: delete rebuild dependency.
# for now must always rebuild
hmod <- sde.make.model(ModelFile = "hestModel.h", rebuild = TRUE,
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

# test drift and diffusion functions
# assuming there are no errors on our end, this is all the user should
# have to do
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)
#X0[1,2] <- -5
#Theta[1,5] <- 5

# R functions
dr.R <- hest.dr(x = X0, theta = Theta)
df.R <- hest.df(x = X0, theta = Theta)
# msde functions (C++)
dr.C <- sde.drift(model = hmod, x = X0, theta = Theta)
df.C <- sde.diff(model = hmod, x = X0, theta = Theta)
# difference
max.diff(x1 = dr.R, x2 = dr.C)
max.diff(x1 = df.R, x2 = df.C)

#--- posterior inference -------------------------------------------------------

# simulate data
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
nObs <- 1e3
dT <- 1/252
hsim <- sde.sim(model = hmod, init.data = x0, params = theta,
                dt = dT, dt.sim = dT/100, N = nObs, nreps = 1)

# initialize data (can no longer provide data directly to sde.post)
m <- 2 # degree of Euler approximation: dt_euler = dt/m
#par.index <- c(2, rep(nObs-1, 1)) # all but first vol unobserved
par.index <- sample(1:2, nObs, replace = TRUE) # randomly observed vol
init.data <- sde.init(data = hsim$data, dt = dT, m = m,
                      par.index = par.index)

# define prior
# default prior is a list with _named_ arguments mu and Sigma
# the names correspond to which elements of (theta, x0) are multivariate normal
# the remaining elements are given flat prior
# NOTE: this is _NOT_ what's described in the document online!
prior <- list(mu = c(.1, .35, 1.0, .5, -.81),
              Sigma = crossprod(matrix(rnorm(25),5)))
prior$Sigma <- sqrt(diag(c(.1, 8, .15, .002, .002))) %*% cov2cor(prior$Sigma)
prior$Sigma <- prior$Sigma %*% sqrt(diag(c(.1, 8, .15, .002, .002)))
names(prior$mu) <- param.names
colnames(prior$Sigma) <- param.names
rownames(prior$Sigma) <- param.names

# mcmc specification
# NEW: adaptive MCMC should tune jump-size automatically!
# acceptance rates should be close to 45%
#mwg.sd <- c(.1, 1, .1, .01, .01) # random walk metropolis for params
#names(mwg.sd) <- param.names
mwg.sd <- NULL
adapt <- TRUE
update.params <- TRUE
update.data <- TRUE
nsamples <- 1e3
burn <- 1e2

# run posterior sampler
hpost <- sde.post(model = hmod,
                  init.data = init.data, init.params = theta,
                  nsamples = nsamples, burn = burn,
                  hyper.params = prior,
                  mwg.sd = mwg.sd, adapt = adapt,
                  update.params = update.params,
                  update.data = update.data)


#--- example 2: Prokaryotic Regulatory Gene Network ----------------------------

param.names <- paste0("gamma", 1:8)
data.names <- c("R", "P", "Q", "D")
pgnet.dr <- function(x, theta, K = 10) {
  x <- matrix(x, ncol = 4)
  theta <- exp(theta)
  if(!is.matrix(theta)) theta <- matrix(theta, nrow = 1)
  y <- matrix(NA, nrow(x), 4)
  y[,4] <- theta[,2] * (K - x[,4]) - theta[,1] * x[,4] * x[,3]
  y[,2] <- theta[,5] * x[,2] * (x[,2]-1)
  y[,3] <- y[,4] + .5 * y[,2]
  y[,1] <- theta[,6] * x[,3]
  y[,3] <- y[,3] - y[,1]
  y[,2] <- 2 * y[,1] - y[,2] + theta[,4] * x[,1] - theta[,8] * x[,2]
  y[,1] <- theta[,3] * x[,4] - theta[,7] * x[,1]
  y
}
pgnet.df <- function(x, theta, K = 10) {
  x <- matrix(x, ncol = 4)
  theta <- exp(theta)
  if(!is.matrix(theta)) theta <- matrix(theta, nrow = 1)
  y <- matrix(0, nrow(x), 16)
  y[,1] <- sqrt(theta[,3] * x[,4] + theta[,7] * x[,1])
  y[,6] <- theta[,8] * x[,2] + 4*theta[,6] * x[,3] + theta[,4] * x[,1] +
    2*theta[,5] * x[,2] * (x[,2]-1)
  y[,15] <- theta[,1] * x[,4] * x[,3] + theta[,2] * (K - x[,4])
  y[,10] <- -2*theta[,6] * x[,3] - theta[,5] * x[,2] * (x[,2]-1)
  y[,16] <- theta[,6] * x[,3] + y[,15] + .5*theta[,5] * x[,2] * (x[,2]-1)
  y[,11] <- y[,16] - y[,10]^2 / y[,6]
  y[,16] <- sqrt(y[,15] - y[,15]^2 / y[,11])
  y[,11] <- sqrt(y[,11])
  y[,15] <- y[,15] / y[,11]
  y[,6] <- sqrt(y[,6])
  y[,10] <- y[,10] / y[,6]
  if(nrow(y) == 1) y <- matrix(y, 4, 4)
  y
}

# test code
pgmod <- sde.make.model(ModelFile = "pgnetModel.h",
                       param.names = param.names,
                       data.names = data.names, rebuild = TRUE)
ndims <- pgmod$ndims
nparams <- pgmod$nparams

theta <- log(c(.1, .7, .35, .2, .1, .9, .3, .1))
names(theta) <- param.names
x0 <- c(R = 5, P = 5, Q = 5, D = 5)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# R code
dr.R <- pgnet.dr(x = X0, theta = Theta)
df.R <- pgnet.df(x = X0, theta = Theta)
# msde
dr.C <- sde.drift(model = pgmod, x = X0, theta = Theta)
df.C <- sde.diff(model = pgmod, x = X0, theta = Theta)

max.diff(x1 = dr.R, x2 = dr.C)
max.diff(x1 = df.R, x2 = df.C)
