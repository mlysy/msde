# test pmcmc based on Ethan's EOU model
require(msde)
source("../../../inst/proj/sde.pmcmc.R")
# ----- extract simulation data -----
load("eou1_1000.rdata")
# ----- initialization -----
dn <- c("X",  "V")
pn <- c("alpha", "gamma", "mu", "sigma", "rho")
nvar <- 1 # number of observed variables
eou <- sde.make.model("EOU.h", data.names = dn, param.names = pn)
theta0 <- sim$params
dT <- sim$dt
path <- sample(1:1000, 1) # select a path
data_obs <- sim$data[,,path]
init <- sde.init(eou, x = data_obs, dt = dT, theta = theta0, nvar.obs = nvar)
# ----- prior on (X, V), V is unknown -----
# fixed.params <- rep(TRUE, eou$nparams)
# names(fixed.params) <- eou$param.names
# hyper <- list(mu = c(0, 0), Sigma = diag(1, 2))
# names(hyper$mu) <- c("X", "V")
# dimnames(hyper$Sigma) <- rep(list(c("X", "V")), 2)
fixed.params <- rep(FALSE, eou$nparams)
hyper <- NULL # flat prior
# ----- MCMC posterior sampling -----
# We do not need to run the following code
# Instead we directly use the data Ethan provides
nsamples <- 1e5
burn <- 1e3
# post <- sde.post(eou, init, hyper = hyper,
#                 fixed.params = rep(FALSE,5),
#                 nsamples = nsamples, burn = burn)
# L.mcmc <- post$params
# ----- particle MCMC posterior sampling -----
npart <- 100
rw.sd <- 1
ppost <- sde.pmcmc(eou, init, hyper,
                  nsamples, burn, rw.sd, fixed.params,
                  last.miss.out = TRUE, adapt = TRUE,
                  npart, resample = "multi", threshold = 0.5)
# ------------------- NA issue ---------------------------------------------------------
# encounter NA lwgt when ii = -980
path <- 456
data_obs <- sim$data[,,path]
init <- sde.init(eou, x = data_obs, dt = dT, theta = theta0, nvar.obs = nvar)
theta.prop = c(-0.473861, 5.197600, -1.514236, 1.094254, -0.993718)
x.prop = c(3.455181, -2.793101)
sde.valid.params(eou, theta.prop)
sde.valid.data(eou, x.prop, theta.prop)
# for testing
test_init <- init
test_init$params <- theta.prop
test_init$data[1,] <- x.prop
pf <- sde.pf(eou, test_init, npart, resample = "multi", threshold = 0.5, history = FALSE)
print(pf$lwgt)
# ---------------------------------------------------------------------------------------

# output
accept <- ppost$acc
print(accept)
L.pmcmc <- ppost$params
# ----- comparison of histograms -----
L.seq <- seq(min(L.mcmc, L.pmcmc), max(L.mcmc, L.pmcmc), len = 500)
hist(L.mcmc, breaks = 100, freq = FALSE,
     main = expression(p(V*" | "*bold(Y)[1])),
     xlab = expression(V))
lines(L.seq, L.pmcmc, col = "red")
legend("topright", legend = c("MCMC", "PMCMC"),
       pch = c(NA, 22), lty = c(1, NA), col = c("black", "red"))
