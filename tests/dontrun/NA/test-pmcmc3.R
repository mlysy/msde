# a script to reproduce the testthat tests

require(msde)
source("smc-testfunctions.R")

# create model
dn <- c("X",  "V")
pn <- c("alpha", "gamma", "mu", "sigma", "rho")
eou <- sde.make.model("EOU.h", data.names = dn, param.names = pn)

# data and parameters
## theta <- c(alpha = -0.47, gamma = 5.20, mu = -1.51, sigma = 1.09, rho = -0.99)
## x <- c(X = 3.46, V = -2.793)
theta0 <- c(-0.473861, 5.197600, -1.514236, 1.094254, -0.9)#93718)
x0 <- c(3.455181, -2.793101)
v0 <- x0[2]
load("eou1_1000.rdata")
dT <- sim$dt
Xt <- sim$data[,,456]
Xt[1,2] <- v0
sde.valid.params(eou, theta0)
sde.valid.data(eou, x0, theta0)

# particle filter
ncomp <- nrow(Xt)
npart <- 100
ndims <- 2
init <- sde.init(model, x = Xt[1:ncomp,], dt = dT,
                 theta = theta0,
                 nvar.obs = 1)
ncomp <- nrow(init$data)
# normal draws
Z <- matrix(rnorm(npart*ndims*(ncomp-1)), ncomp-1, npart*ndims)

# pf in R
pf.R <- pf.fun(init, dr = drift.fun, df = diff.fun, Z = Z,
               history = TRUE)

# pf in C++
pf.cpp <- sde.pf(model = model, init = init, npart = npart,
                 resample = "multi", threshold = .5,
                 Z = array(c(Z), c(ncomp-1, ndims, npart)),
                 history = FALSE)
anyNA(pf.cpp$lwgt)

range(pf.R$lwgt - pf.cpp$lwgt)
