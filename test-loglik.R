#--- speed/accuracy comparisons with old version of package ----------------

# the old msde package was written without OOP and (*) merged all cpp/h
# files into 1 before compile on-the-fly with Rcpp
# the new version is written with OOP and checks the speed effect of (*)

require(msde) # old version
require(msdeHeaders) # new version
pkg.names <- c("msde", "msdeHeaders")
max.diff <- function(x1, x2) {
  c(abs = max(abs(x1-x2)), rel = max(abs(x1-x2)/abs(x1)))
}
pkg1 <- function(fun) {
  eval(parse(text = paste0("msde::", fun)))
}
pkg2 <- function(fun) {
  eval(parse(text = paste0("msdeHeaders::", fun)))
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
hmod <- pkg1("sde.make.model")(list = hestList,
                               showOutput = TRUE, rebuild = TRUE)

# new version (msdeHeaders -- interface still in progress)
#ndims <- 2
#nparams <- 5
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
#hmod2 <- file.path(pkg2(":.msdeCppPath"), "hestModel.cpp")
hmod2 <- pkg2("sde.make.model")(ModelFile = "hestModel.h",
                                param.names = param.names,
                                data.names = data.names,
                                showOutput = TRUE, rebuild = TRUE)
ndims <- hmod2$ndims
nparams <- hmod2$nparams

#--- loglikelihood speed comp ----------------------------------------------

nReps <- 1e3
nObs <- 1e3
dT <- 1/252

# initialize
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# generate data
hsim <- pkg1("sde.sim")(model = hmod, init.data = X0, params = Theta,
  dt = dT, dt.sim = dT, N = nObs, nreps = nReps)
X1 <- hsim$data
X2 <- aperm(X1, c(2, 3, 1))

# loglikelihood calcs
time.p1 <- system.time({
  ll.p1 <- pkg1("sde.loglik")(model = hmod, x = X1, dt = dT,
                              theta = Theta)
})
time.p2 <- system.time({
  ll.p2 <- pkg2("sde.loglik")(model = hmod2, debug = FALSE,
                              x = X2, dt = dT,
                              theta = Theta)
})
max.diff(ll.p1, ll.p2)
tmp <- c(time.p1[3], time.p2[3])
names(tmp) <- pkg.names
tmp/tmp[1]
