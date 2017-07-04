#--- check the multivariate normal prior ---------------------------------------

devtools::document()
devtools::install()

require(msde)
require(mvtnorm)

param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

nrv <- 2
rv.names <- sample(c(param.names, data.names), nrv, replace = FALSE)
mu <- rnorm(nrv)
Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
names(mu) <- rv.names
colnames(Sigma) <- rv.names
rownames(Sigma) <- rv.names
hyper <- list(mu = mu, Sigma = Sigma)

# check parser
mvn.hyper.check(hyper = hyper,
                param.names = param.names, data.names = data.names)

# check evaluator
nrv <- 3
if(nrv > 0) {
  rv.names <- sample(c(param.names, data.names), nrv, replace = FALSE)
  mu <- rnorm(nrv)
  Sigma <- crossprod(matrix(rnorm(nrv^2),nrv,nrv))
  names(mu) <- rv.names
  colnames(Sigma) <- rv.names
  rownames(Sigma) <- rv.names
  hyper <- list(mu = mu, Sigma = Sigma)
} else hyper <- NULL

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
                 hyper = hyper)
range(diff(lpR-lpC))

#--- custom prior --------------------------------------------------------------

# random matrix of size nreps x length(x) from vector x
jit.vec <- function(x, nreps) {
  apply(t(replicate(n = nreps, expr = x, simplify = "matrix")), 2, jitter)
}
# maximum absolute and relative error between two arrays
max.diff <- function(x1, x2) {
  c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/max(abs(x1), 1e-8)))
}
# input check
lvcheck <- function(hyper, param.names, data.names) {
  if(is.null(names(hyper)) ||
     !identical(sort(names(hyper)), c("mu", "sigma"))) {
    stop("hyper must be a list with elements mu and sigma.")
  }
  mu <- hyper$mu
  if(length(mu) == 1) mu <- rep(mu, 4)
  if(!is.numeric(mu) || length(mu) != 4) {
    stop("mu must be a numeric scalar or vector of length four.")
  }
  sig <- hyper$sigma
  if(length(sig) == 1) sig <- rep(sig, 4)
  if(!is.numeric(sig) || length(sig) != 4 || !all(sig > 0)) {
    stop("sigma must be a positive scalar or vector of length four.")
  }
  list(mu, sig)
}

data.names <- c("H", "L")
param.names <- c("alpha", "beta", "gamma")
lvmod <- sde.make.model(ModelFile = "lotvolModel.h",
                        PriorFile = "lotvolPrior.h", # prior specification
                        hyper.check = lvcheck, # prior input checking
                        data.names = data.names,
                        param.names = param.names)


# generate some test values
nreta <- 12
x0 <- c(H = 71, L = 79)
theta0 <- c(alpha = .5, beta = .0025, gamma = .3)
X <- jit.vec(x0, nreta)
Theta <- jit.vec(theta0, nreta)
Eta <- cbind(Theta, L = X[,"L"])
nrphi <- 5
Phi <- lapply(1:nrphi, function(ii) list(mu = rnorm(4), sigma = rexp(4)))

# prior check

# R version
lpi.R <- matrix(NA, nreta, nrphi)
for(ii in 1:nrphi) {
  lpi.R[,ii] <- colSums(dlnorm(x =  t(Eta),
                               meanlog = Phi[[ii]]$mu,
                               sdlog = Phi[[ii]]$sigma, log = TRUE))
}

# C++ version
lpi.cpp <- matrix(NA, nreta, nrphi)
for(ii in 1:nrphi) {
  lpi.cpp[,ii] <- sde.prior(model = lvmod, theta = Theta, x = X,
                            hyper = Phi[[ii]])
}

# compare
max.diff(lpi.R, lpi.cpp)


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


