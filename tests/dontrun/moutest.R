library(msde)
source("mOU-Kalman-revised.R")
require("mvtnorm")
source("test-post_function.R")
# Create Lambda
ndims <- 1
Lambda <- 1.5205

Gamma <- matrix(-0.064364,byrow=TRUE,nrow=ndims,ncol=ndims)

Psi <- 1.4469
theta <- c(theta = runif(1,min=0.1,max=0.9))

ndims <- length(Lambda) # length of Lambda == number of dimension

	# Initialize param.names
param.names <- names(theta)

	# Initialize data.names
data.names <- c("X1")

	# Initialize model in C++
message("Initializing Model in C++")
moumod <- sde.make.model(ModelFile = "mOUModel.h",
                         data.names = data.names,
                         param.names = param.names,
                         rebuild=TRUE)

message("Model initialization complete.")

	# check init
	# Make a list of sde.init object
	# generate data
X0 <- runif(ndims, min = -1, max = 1)
dT <- 1/250

		# generate nobs, burn
burn <- 0
nobs.sim <- 2000
mou.sim <- sde.sim(model = moumod, x0 = X0, theta = theta,
                   dt = dT, dt.sim = dT, nobs = nobs.sim, burn=burn)
