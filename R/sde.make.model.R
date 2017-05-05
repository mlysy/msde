#' Create an SDE Model Object
#'
#' @param ModelFile Path to the header file where the sde model is defined.
#' @param PriorFile Path to the header file where the prior is defined.
#' @param data.names Optional vector of names for the components of the sde.  Defaults to \code{X1,...,Xd}.
#' @param param.names Optional vector of names for the parameters of the sde.  Defaults to \code{theta1,...,thetap}.
#' @param ... additional parameters that are passed to \code{sourceCpp} when compiling the c++ code
#'@return An \code{sde.model} object, which is a list containing the following elements:
#' \itemize{
#' \item \code{ndims}, \code{nparams}: the number of sde components and parameters.
#' \item \code{data.names}, \code{param.names}: the names of the sde components and parameters.
#' \item \code{drift}, \code{diff}: functions to evaluate the drift and diffusion.
#' \item \code{loglik}: the loglikelihood.
#' \item \code{logprior}: the custom log prior.
#' \item \code{sim}: function for simulating data.
#' \item \code{post}: MCMC sampler for posterior distribution.
#' }
#' @details The functions \code{sim}, \code{post}, \code{drift}, \code{diff}, \code{logpior}, and \code{loglik} should never be called directly. Instead use \code{sde.sim}, \code{sde.post} \code{sde.diff}, \code{sde.drift} and \code{sde.loglik}.
#'
#' The code is compiled by copying the \code{ModelFile} to the \code{tmpdir} directory, along with a wrapper \code{.cpp} file to be compiled by \code{Rcpp::sourceCpp}.
#' @note TODO: fix \code{rebuild = TRUE}.  This could be done with \code{utils::changedFiles}.
#'@export
sde.make.model <- function(ModelFile, PriorFile = "default",
                           data.names, param.names, prior.spec,
                           ..., debug = FALSE) {
  if(debug) browser()
  #.msdeCppPath <- "C:/Users/Jerome/Documents/R/library/msdeHeaders/cppTemplates"
  sde.model <- list()
  # prior specification
  if(PriorFile == "default") {
    PriorFile <- file.path(.msdeCppPath, "mvnPrior.h")
    if(!missing(prior.spec)) {
      warning("Custom prior.spec ignored for default prior.")
    }
    prior.spec <- mvn.prior.spec
  } else {
    if(missing(prior.spec)) {
      stop("Must provide prior.spec for custom prior.")
    }
  }
  file.copy(from = PriorFile,
            to = file.path(tempdir(), "sdePrior.h"),
            overwrite = TRUE, copy.date = TRUE)
  # compile C++ code
  file.copy(from = file.path(.msdeCppPath, "sdeUtils.cpp"),
            to = file.path(tempdir(), "sdeUtils.cpp"),
            overwrite = TRUE, copy.date = TRUE)
  ## if(missing(ModelFile)) {
  ##   ModelFile <- file.path(.msdeCppPath, "sdeModel.h")
  ## }
  file.copy(from = ModelFile,
            to = file.path(tempdir(), "sdeModel.h"),
            overwrite = TRUE, copy.date = TRUE)
  sourceCpp(file = file.path(tempdir(), "sdeUtils.cpp"),
            env = environment(), ...)
  environment(sde.model$sim) <- globalenv()
  #environment(sde.model$post) <- globalenv()
  environment(sde.model$drift) <- globalenv()
  environment(sde.model$diff) <- globalenv()
  environment(sde.model$loglik) <- globalenv()
  environment(sde.model$logprior) <- globalenv()
  # extract ndims and nparams
  ndims <- ndims()
  nparams <- nparams()
  # parameter and data names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) stop("Incorrect data.names.")
  if(length(param.names) != nparams) stop("Incorrect param.names.")
  sde.model <- c(sde.model,
                 list(ndims = ndims, nparams = nparams,
                      data.names = data.names, param.names = param.names,
                      prior.spec = prior.spec))
  # output
  class(sde.model) <- "sde.model"
  sde.model
}
