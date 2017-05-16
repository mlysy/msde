#--- test omp implementation ---------------------------------------------------

devtools::document()
devtools::install()

require(msdeHeaders)

# build model
param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names,
                       omp = TRUE, showOutput = TRUE, rebuild = TRUE)
ndims <- hmod$ndims
nparams <- hmod$nparams

# generate data
theta <- c(alpha = 0.1, gamma = 1, beta = 0.8, sigma = 0.6, rho = -0.8)
x0 <- c(X = log(1000), Z = 0.1)
nReps <- 1e7
Theta <- matrix(theta, nReps, nparams, byrow = TRUE)
colnames(Theta) <- param.names
X0 <- matrix(x0, nReps, ndims, byrow = TRUE)
colnames(X0) <- data.names


# ncores for sde.drift now DISABLED.
# best to parallelize this from within R
## ncores <- 8 # set to < 0 or NA to disable omp
## system.time({
##   dr <- sde.drift(model = hmod, x = X0, theta = Theta,
##                   ncores = ncores, debug = FALSE)
## })
