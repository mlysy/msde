#--- check the multivariate normal prior ---------------------------------------

# this is a basic check, need to interface later

require(msdeHeaders)
require(mvtnorm)

# Rcpp::sourceCpp(file.path(msdeHeaders:::.msdeCppPath, "sdePrior.cpp"))
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

# generate data
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# hyperparameters
nrv <- 7
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))

nac <- 0
iacc <- (nrv-nac+1):nrv
if(nac > 0) {
  lp1 <- dmvnorm(x = cbind(X0, Theta)[,iacc,drop=FALSE],
                 mean = mu[iacc], sigma = Sigma[iacc,iacc,drop=FALSE],
                 log = TRUE)
  phi <- list(mu = mu[iacc], Sigma = Sigma[iacc,iacc])
} else {
  lp1 <- rep(0, nReps)
  phi <- NULL
}
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)
lp1-lp2


#--- check prior.spec ----------------------------------------------------------

nrv <- 20
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
phi <- list(mu = mu, Sigma = Sigma)
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)

nrv <- 6
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
phi <- list(mu, Sigma)
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)


nrv <- 4
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
phi <- list(mu, raindeer = Sigma)
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)

nrv <- 6
mu <- rnorm(nrv-1)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
phi <- list(mu = mu, Sigma = Sigma)
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)


nrv <- 6
mu <- rnorm(nrv-1)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
phi <- list(mu = NULL, Sigma = NULL)
lp2 <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = phi)


