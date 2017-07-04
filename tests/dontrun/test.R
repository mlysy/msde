#--- basic tests of the msde header-only library ---------------------------

require(msde)

# does the package work?

hest.model <- sde.make.model(data.names = c("X", "Z"),
                             param.names = c("alpha", "gamma", "beta",
                               "sigma", "rho"),
                             verbose = TRUE, rebuild = TRUE, debug = FALSE)

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

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

nReps <- 10
Theta <- apply(t(replicate(nReps, theta)), 2, jitter)
X0 <- apply(t(replicate(nReps, x0)), 2, jitter)

# R code
dr.R <- hest.dr(X = X0[,1], Z = X0[,2], theta = Theta)
df.R <- hest.df(X = X0[,1], Z = X0[,2], theta = Theta)

# c++ code
dr.cpp <- sde.drift(model = hest.model, x = X0, theta = Theta)
df.cpp <- sde.diff(model = hest.model, x = X0, theta = Theta)

range(dr.R-dr.cpp)
range(df.R-df.cpp)

#--- simulation ------------------------------------------------------------

theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)

# each of these takes about 15s on my lappy
nObs <- 1e3
nReps <- 1e3
dT <- 1/252
if(same.rnd) set.seed(SEED)
time.m <- system.time({
  hsim.m <- pkg1("sde.sim")(model = hmod, init.data = x0, params = theta,
                            dt = dT, dt.sim = dT/100,
                            N = nObs, nreps = nReps)
})
if(same.rnd) set.seed(SEED)
time.mt <- system.time({
  hsim.mt <- pkg2("sde.sim")(model = hmod2, init.data = x0, params = theta,
                             dt = dT, dt.sim = dT/100,
                             N = nObs, nreps = nReps)
})
if(same.rnd) set.seed(SEED)
time.mto <- system.time({
  hsim.mto <- hest.sim(model = hmod2, init.data = x0, params = theta,
                       dt = dT, dt.sim = dT/100, N = nObs, nreps = nReps)
})
