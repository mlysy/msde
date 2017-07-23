#' R Methods
#'
#' @name Methods
#' @export
Robj.eval <- function(obj, y) {
  if(class(obj) != "Robj") stop("obj must be an Robj object.")
  .Cobj_eval(obj$ptr, y)
}

#' @name Methods
#' @export
Robj.set.x <- function(obj, x) {
  if(class(obj) != "Robj") stop("obj must be an Robj object.")
  .Cobj_setx(obj$ptr, x)
}

#' @name Methods
#' @export
Robj.get.x <- function(obj) {
  if(class(obj) != "Robj") stop("obj must be an Robj object.")
  .Cobj_getx(obj$ptr)
}

#' @name Methods
#' @export
sde.size <- function(model) {
  if(class(model) != "sde.model") stop("model must be an sde.model object.")
  nparams <- .sde_nParams(model$ptr)
  ndims <- .sde_nDims(model$ptr)
  c(nparams = nparams, ndims = ndims)
}
