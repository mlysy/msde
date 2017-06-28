#' Create an SDE Model Object
#'
#' @param ModelFile Path to the header file where the sde model is defined.
#' @param PriorFile Path to the header file where the prior is defined.  See \code{\link{sde.prior}} for details.
#' @param data.names Optional vector of names for the components of the sde.  Defaults to \code{X1,...,Xd}.
#' @param param.names Optional vector of names for the parameters of the sde.  Defaults to \code{theta1,...,thetap}.
#' @param hyper.check A function with arguments \code{prior.args}, \code{param.names}, and \code{data.names} used for passing the model hyper parameters to the C++ code.  See \code{\link{mvn.hyper.check}} for details.
#' @param OpenMP This sets whether the model is compiled for usage in \code{OpenMP}.
#' @param ... additional parameters that are passed to \code{Rcpp::sourceCpp} when compiling the C++ code.
#'@return An \code{sde.model} object, which is a list containing the following elements:
#' \itemize{
#' \item \code{ndims}, \code{nparams}: the number of sde components and parameters.
#' \item \code{data.names}, \code{param.names}: the names of the sde components and parameters.
#' \item \code{drift}, \code{diff}: functions to evaluate the drift and diffusion.
#' \item \code{loglik}: the loglikelihood.
#' \item \code{logprior}: the custom log prior.
#' \item \code{sim}: function for simulating data.
#' \item \code{post}: MCMC sampler for posterior distribution.
#' \item \code{omp}: a logical flag for whether or not the model was compiled for multicore functionality with \code{OpenMP}.
#' }
#' @details The functions \code{sim}, \code{post}, \code{drift}, \code{diff}, \code{logpior}, and \code{loglik} should never be called directly. Instead use \code{sde.sim}, \code{sde.post} \code{sde.diff}, \code{sde.drift} and \code{sde.loglik}.
#'
#' The code is compiled by copying the \code{ModelFile} to the \code{tmpdir} directory, along with a wrapper \code{.cpp} file to be compiled by \code{Rcpp::sourceCpp}.
#' @examples
#' \donttest{
#' hex <- example.models("hest")
#' sde.make.model(ModelFile = hex$ModelFile,
#'                param.names = hex$param.names,
#'                data.names = hex$data.names)
#'
#' }
#'@export
sde.make.model <- function(ModelFile, PriorFile = "default",
                           data.names, param.names, hyper.check,
                           OpenMP = FALSE, ...) {
  sde.model <- list()
  # prior specification
  if(PriorFile == "default") {
    PriorFile <- file.path(.msdeCppPath, "mvnPrior.h")
    if(!missing(hyper.check)) {
      warning("Custom hyper.check ignored for default prior.")
    }
    hyper.check <- mvn.hyper.check
  } else {
    if(missing(hyper.check)) {
      stop("Must provide hyper.check for custom prior.")
    }
    if(!identical(formalArgs(hyper.check),
                  c("hyper", "param.names", "data.names"))) {
      stop("hyper.check must have formal arguments: hyper, param.names, data.names.")
    }
  }
  # compile C++ code
  rebuild <- .copy.cpp.files(ModelFile, PriorFile)
  cpp.args <- list(...)
  if(is.null(cpp.args$rebuild) || !cpp.args$rebuild) {
    cpp.args$rebuild <- rebuild
  }
  cpp.args <- c(list(file = file.path(tempdir(), "msdeExports.cpp"),
                     env = environment()), cpp.args)
  # OpenMP support
  if(OpenMP) old.env <- .omp.set()
  # if(debug) browser()
  do.call(sourceCpp, cpp.args)
  if(OpenMP) .omp.unset(env = old.env)
  environment(sde.model$sim) <- globalenv()
  environment(sde.model$post) <- globalenv()
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
  if(length(data.names) != ndims) {
    stop("data.names has wrong length.")
  }
  if(length(param.names) != nparams) {
    stop("param.names has wrong length.")
  }
  sde.model <- c(sde.model,
                 list(ndims = ndims, nparams = nparams,
                      data.names = data.names, param.names = param.names,
                      hyper.check = hyper.check, omp = OpenMP))
  # output
  class(sde.model) <- "sde.model"
  sde.model
}

# flags appropriate errors and returns T/F of whether
# both files have any changes
.copy.cpp.files <- function(ModelFile, PriorFile) {
  rebuild <- FALSE
  # prior
  fname <- file.path(tempdir(), "sdePrior.h")
  if(file.exists(fname)) {
    fold <- readLines(con = fname)
    fnew <- readLines(con = PriorFile)
    rebuild <- !identical(fold, fnew)
  }
  flag <- file.copy(from = PriorFile,
                    to = fname,
                    overwrite = TRUE, copy.date = TRUE)
  if(!flag) {
    stop("PriorFile \"", PriorFile, "\" not found.")
  }
  # model
  fname <- file.path(tempdir(), "sdeModel.h")
  if(file.exists(fname)) {
    fold <- readLines(con = fname)
    fnew <- readLines(con = ModelFile)
    rebuild <- rebuild || !identical(fold, fnew)
  }
  flag <- file.copy(from = ModelFile,
                    to = fname,
                    overwrite = TRUE, copy.date = TRUE)
  if(!flag) {
    stop("ModelFile \"", ModelFile, "\" not found.")
  }
  # export file
  if(!file.exists(file.path(.msdeCppPath, "msdeExports.cpp"))) {
    msdeExports <- c(readLines(file.path(.msdeCppPath, "sdeUtils.cpp")),
                     readLines(file.path(.msdeCppPath, "sdeSim.cpp")),
                     readLines(file.path(.msdeCppPath, "sdePost.cpp")))
    cat(msdeExports, sep = "\n",
        file = file.path(.msdeCppPath, "msdeExports.cpp"))
  }
  file.copy(from = file.path(.msdeCppPath, "msdeExports.cpp"),
            to = file.path(tempdir(), "msdeExports.cpp"),
            overwrite = TRUE, copy.date = TRUE)
  rebuild
}

#--- omp set and unset ---------------------------------------------------------

# adds -fopenmp flags to PKG_CXXFLAGS and PKG_LIBS

.omp.set <- function() {
  cxx <- Sys.getenv(x = "PKG_CXXFLAGS", unset = NA)
  libs <- Sys.getenv(x = "PKG_LIBS", unset = NA)
  env <- list(cxx = cxx, libs = libs)
  Sys.setenv(PKG_CXXFLAGS = ifelse(is.na(cxx),
                                "-fopenmp", paste("-fopenmp", cxx)))
  Sys.setenv(PKG_LIBS = ifelse(is.na(libs),
                               "-fopenmp", paste("-fopenmp", libs)))
  env
}

.omp.unset <- function(env) {
  if(is.na(env$cxx)) {
    Sys.unsetenv(x = "PKG_CXXFLAGS")
  } else {
    Sys.setenv(PKG_CXXFLAGS = env$cxx)
  }
  if(is.na(env$libs)) {
    Sys.unsetenv(x = "PKG_LIBS")
  } else {
    Sys.setenv(PKG_LIBS = env$libs)
  }
}
