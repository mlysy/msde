#--- check the multivariate normal prior ---------------------------------------

devtools::document()
devtools::install()

require(msdeHeaders)
require(mvtnorm)

param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

nrv <- 1
rv.names <- sample(c(param.names, data.names), nrv, replace = FALSE)
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
names(mu) <- rv.names
colnames(Sigma) <- rv.names
rownames(Sigma) <- rv.names
prior.args <- list(mu = mu, Sigma = Sigma)

# check parser
mvn.prior.spec(prior.args = prior.args, debug = FALSE,
               param.names = param.names, data.names = data.names)

# check evaluator
nrv <- 0
if(nrv > 0) {
  rv.names <- sample(c(param.names, data.names), nrv, replace = FALSE)
  mu <- rnorm(nrv)
  Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
  names(mu) <- rv.names
  colnames(Sigma) <- rv.names
  rownames(Sigma) <- rv.names
  prior.args <- list(mu = mu, Sigma = Sigma)
} else prior.args <- NULL

# generate data
nReps <- 10
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# R
if(nrv > 0) {
lpR <- dmvnorm(x = cbind(X0, Theta)[,rv.names,drop=FALSE],
               mean = mu[rv.names],
               sigma = Sigma[rv.names,rv.names,drop=FALSE],
               log = TRUE)
} else lpR <- rep(0, nReps)
# C++
lpC <- sde.prior(model = hmod, x = X0, theta = Theta,
                 phi = prior.args)
range(diff(lpR-lpC))

#--- old prior -----------------------------------------------------------------

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


