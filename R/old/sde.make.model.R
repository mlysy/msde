#'@name sde.make.model
#'@title Model Creation
#'@description Calling \code{sde.make.model} will compile the c++ code provided in its parameters along with
#'             c++ code included in this package in order to produce a shared library.
#'             The functions of this shared library should be accessed by passing the result of \code{sde.make.model}
#'             into any of the functions \code{sde.sim}, \code{sde.post}, \code{sde.diff}, \code{sde.drift},
#'             or \code{sde.loglik}
#'@param list \code{NULL} by default. Rather than passing model parameters to \code{sde.make.model} individually,
#'             the parameters may be passed as a list. Any combination of model parameters may be provided in this
#'             list (excluding debug and ...).
#'             Values from this list will only be used if no coresponding value is provided as a parameter to \code{sde.make.model}
#'@param ndims an integer giving the number of dimensions for the model
#'@param nparams an integer giving the number of parameters for the model
#'@param definitionFilePath The full file path of a file containing c++ definitions for \code{sdeDr}, \code{sdeDf},
#'                          \code{isValidData} and \code{isValidParams}\cr
#'                          If \code{definitionFilePath} is provided, any values for \code{sdeDr}, \code{sdeDf},
#'                          \code{isValidData} or \code{isValidParams} will be ignored.
#'                          The file must contain definitions for all 4 functions.
#'@param sdeDr A string of c++ code giving a function definition for the model's drift function. \cr
#'             The function's signature must be: void sdeDr(double dr[], double x[], double t, double params[])
#'@param sdeDf A string of c++ code giving a function definition for the model's diffusion function. \cr
#'             The function's signature must be: void sdeDf(double df[], double x[], double sqrtT, double params[])
#'@param isValidData A string of c++ code giving a function definition for the model's data validation function.\cr
#'             The function's signature must be: int isValidData(double x[]) \cr
#'             The return value must be 1 for valid data, 0 for invalid data
#'@param isValidParams A string of c++ code giving a function definition for the model's parameter validation function. \cr
#'             The function's signature must be: int isValidParams(double params[]) \cr
#'             The return value must be 1 for valid parameters, 0 for invalid parameters
#'@param logPrior A string of c++ code defining a custom prior for the model.\cr
#'             This parameter is optional.\cr
#'             If it is not provided, the values of CustomPriorClass, CustomPriorConstructor and .CustomPriorFunction
#'             will be used. \cr
#'             If it is provided, it must consist of a complete class definition meeting the criteria specified
#'             in the package vignette. \cr
#'             Custom priors are only relevant for those users who do not wish to use any of the provided priors
#'             (flat, normal and Gaussian Copula)
#'@param data.names a vector of names for the columns of data that will be used by the model
#'@param param.names a vector of names for the parameters that will be used by the model
#'@param ... additional parameters that are passed to sourceCpp when compiling the c++ code
#'@param debug a boolean (\code{FALSE} by default) if set to \code{TRUE}, will cause the function to open a browser mid-call
#'@return an sde.model containing:\cr \cr
#'                          ndims, \cr
#'                          nparams, \cr
#'                          logPrior, \cr
#'                          data.names, \cr
#'                          param.names, \cr
#'                          sim (simulation function for the model), \cr
#'                          post (posterior inference function for the model), \cr
#'                          drift (drift function for the model), \cr
#'                          diff (diffusion function for the model), \cr
#'                          loglik (loglikelihood function for the model) \cr \cr
#'
#'                          The functions \code{sim}, \code{post}, \code{drift}, \code{diff} and \code{loglik}
#'                          should never be called directly. Instead use \code{sde.sim}, \code{sde.post},
#'                          \code{sde.diff}, \code{sde.drift} and \code{sde.loglik}
#'@examples
#'hest.model <- sde.make.model(list = hestList, cpp.out = TRUE)
#'names(hest.model)
#'# preview the C++ code
#'message(paste(hest.model$cpp.code[1:10], collapse = "\n"))
#'@export
sde.make.model <- function(list = NULL, ndims, nparams, data.names, param.names, custom.names,
                           sdeDr, sdeDf, isValidData, isValidParams, logCustomPrior,
                           cpp.out = TRUE, ..., debug = FALSE) {
  # list supplies values not otherwise specified
  nm <- c("ndims", "nparams", "data.names", "param.names", "custom.names",
          "sdeDr", "sdeDf", "isValidData", "isValidParams", "logCustomPrior")
  for(ii in nm) {
    eval(substitute(if(missing(ii) && !is.null(list$ii))
                    ii <- list$ii, list(ii = as.symbol(ii))))
  }
  # default data and parameter names
  if(missing(data.names)) data.names <- paste0("X", 1:ndims)
  if(missing(param.names)) param.names <- paste0("theta", 1:nparams)
  if(length(data.names) != ndims) stop("Incorrect data.names.")
  if(length(param.names) != nparams) stop("Incorrect param.names.")
  # check for custom prior
  if(missing(logCustomPrior)) logCustomPrior <- NULL
  if(missing(custom.names)) custom.names <- NULL
  sde.model <- list(ndims = ndims, nparams = nparams,
                    data.names = data.names, param.names = param.names,
                    sdeDr = sdeDr, sdeDf = sdeDf, isValidData = isValidData, isValidParams = isValidParams)
  if(!is.null(logCustomPrior)) {
    sde.model <- c(sde.model, list(logCustomPrior = logCustomPrior))
    if(!is.null(custom.names)) sde.model <- c(sde.model, list(custom.names = custom.names))
  } else {
    logCustomPrior <- "double CustomPrior::logPrior(double params[], double x[]) {
  return(0.0);
}"
  }
  # create c++ code
  h.nm <- c("Priors.h", "sdeCore.h", "ConvenienceFunctions.h", "sdeAPI-Rcpp.h")
  cpp.nm <- c("ConvenienceFunctions.cpp", "sdeCore.cpp", "Priors.cpp", "sdeAPI-Rcpp.cpp")
  hLines <- sapply(h.nm, function(nm) {
    con <- file(file.path(.msdeCppPath, nm), "r")
    zz <- readLines(con)
    close(con)
    zz
  })
  cppLines <- sapply(cpp.nm, function(nm) {
    con <- file(file.path(.msdeCppPath, nm), "r")
    zz <- readLines(con)
    close(con)
    zz
  })
  userLines <- c(sdeDr, sdeDf, isValidData, isValidParams, logCustomPrior)
  cpp.code <- c(paste0("const int nParams = ", nparams, ";"),
                paste0("const int nDims = ", ndims, ";"),
                unlist(hLines), userLines, unlist(cppLines),"\n")
  if(cpp.out) sde.model <- c(sde.model, list(cpp.code = cpp.code))

  # compile c++ code
  if(debug) browser()
  sourceCpp(code = paste(cpp.code, collapse = "\n"), env = environment(), ...)
  environment(sde.model$sim) <- globalenv()
  environment(sde.model$post) <- globalenv()
  environment(sde.model$drift) <- globalenv()
  environment(sde.model$diff) <- globalenv()
  environment(sde.model$loglik) <- globalenv()

  class(sde.model) <- "sde.model"
  sde.model
}
