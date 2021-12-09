#' R6 class wrapping `sdeRobj<{{sdeModel}}, {{sdePrior}}>`.
#'
#' @details This class is intended to be hidden, i.e., wouldn't be exported from a package.  Rather, it is the instantiated object which would potentially be exported, since it would be called with the `sde.{method}` generics provided by **msde**.
#'
#' That being said, the package would probably want to export the class constructor, so as to instantiate new models...
#' @noRd
{{R6ClassName}} <- R6::R6Class(

  classname = "sde.model",

  private = list(
    .ptr = NULL,
    .ndims = NULL,
    .nparams = NULL,
    .data.names = NULL,
    .param.names = NULL,
    .omp = NULL
  ),

  active = list(

    #' @field ptr Pointer to C++ `sdeRobj<{{sdeModel}}, {{sdePrior}}>`.
    ptr = function(value) {
      if(missing(value)) {
        private$.ptr
      } else {
        stop("$ptr is read-only.", call. = FALSE)
      }
    },

    #' @field ndim Number of SDE components.
    ndim = function(value) {
      if(missing(value)) {
        private$.ndim
      } else {
        stop("$ndim is read-only.", call. = FALSE)
      }
    },

    #' @field nparams Number of model parameters.
    nparams = function(value) {
      if(missing(value)) {
        private$.nparams
      } else {
        stop("$nparams is read-only.", call. = FALSE)
      }
    },

    #' @field data.names Names of the SDE components.
    data.names = function(value) {
      if(missing(value)) {
        private$.data.names
      } else {
        if(length(value) != self$ndims) {
          stop("$data.names has wrong length.")
        } else {
          private$.data.names <- value
        }
      }
    },

    #' @field param.names Names of the model parameters.
    param.names = function(value) {
      if(missing(value)) {
        private$.param.names
      } else {
        if(length(value) != self$nparams) {
          stop("$param.names has wrong length.")
        } else {
          private$.param.names <- value
        }
      }
    },

    #' @field omp A logical flag for whether or not the model was compiled for multicore functionality with \code{OpenMP}.
    omp = function(value) {
      if(missing(value)) {
        private$.omp
      } else {
        stop("$omp is read-only.", call. = FALSE)
      }
    }
  ),

  public = list(
    isData = {{RClassName}}_isData,
    isParams = {{RClassName}}_isParams,
    Drift = {{RClassName}}_Drift,
    Diff = {{RClassName}}_Diff,
    Loglik = {{RClassName}}_Loglik,
    Prior = {{RClassName}}_Prior,
    Post = {{RClassName}}_Post,
    hyper.check = NULL,

    initialize = function(data.names, param.names, hyper.check,
                          OpenMP = FALSE) {
      # initialize sdeRobj<{{sdeModel}}, {{sdePrior}}> object in C++
      private$.ptr <- {{RClassName}}_ctor()
      # set ndims and nparams
      private$.ndims <- {{RClassName}}_nDims(self$ptr)
      private$.nparams <- {{RClassName}}_nParams(self$ptr)
      # parameter and data names
      if(missing(data.names)) data.names <- paste0("X", 1:self$ndims)
      if(missing(param.names)) param.names <- paste0("theta", 1:self$nparams)
      self$data.names <- data.names
      self$param.names <- param.names
      # remaining members
      self$hyper.check <- hyper.check
      private$.omp <- OpenMP
    }
  )
)
