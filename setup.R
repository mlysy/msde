#--- compile --------------------------------------------------------------

require(Rcpp)
#require(devtools)

devtools::document()
devtools::install()

devtools::build()

#--- compile test ----------------------------------------------------------

require(msde)

param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(param.names = param.names,
                       data.names = data.names,
                       showOutput = TRUE, rebuild = TRUE)


ndims <- hmod$ndims
nparams <- hmod$nparams

#--- header-only version of msde -----------------------

# TODO: still somewhat slower than no-op version
# interface for param.names and data.names
# check a few more models (ghest, gl05, lotka-volterra)

# Hierarchy:

# LinAlgUtils.h
#
# * chol_decomp
# * v_mult
# * U_mult

# mvnUtils.h
#
# * xmvn
# * zmvn
# * lmvn
# * propMV, a class for mean and variance vectors and associated temporary storage for the *mvn calculations.

# sdeModel.h
#
# this is a class which defines the sde model.
# public members:
# nParams, nDims, sdeDr, sdeDf, isValidData, isValidParams
#
# all but sdeDr and sdeDf are static members.  but these two "mean" and
# "variance" members should be constructed for each data point.
# this is because sometimes, sdeDr and sdeDf will require large-ish
# matrices for temporary storage of intermediate calculations.  so each
# data point should have access to its own storage, and i think it should
# be dynamically allocated.  see sdeData

# sdeLogLik.h
#
# this is a class which defines the complete data and whatever it needs
# to evaluate the log-density.
# the log-density of an sde is essentially the sum of conditional normals,
# each having its own mean and variance for every point.
# i.e.,
# logDens(x | theta) = sum dnorm(x_n+1 | mu(x_n, theta), sd(x_n, theta)).
#
# so the constructor should create an array of sdeModels, each with enough
# storage to compute log-densities.
#
# public members:
# * nComp, tSeq, dT, sqrtDT
# * Drift(x, t, theta, i), Diff(x, t, theta, i)
# * EulerMV(x, theta, i)
# * loglik(x, theta)

# sdePrior.h


# sdeMCMC.h
#
# each MCMC object contains a "transition operator" member and whatever
# memory it needs to do this.
# so for instance:
# * MissGibbs creates storage propMV for each data point.  Its update
#   function is update(x, theta, rvJump, accept)
# * ParamVanilla creates propMV just for the parameters.  it has the same
#   update function.
# * another proposal might compute regression-based updates...
#
# so the constructor should contain a pointer to an sdeLogLik object, the
# fixedParams logical array, and some information on the prior...
