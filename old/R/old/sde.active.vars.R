#' Extract the active variables from an SDE model
#'
#' @param vars Character vector.
#' @param fixed.params Logical vector, which parameters are fixed.
#' @param nmiss Number of unobserved SDE components at given time point.
#' @return Logical vector of length \code{nparams+ndims} indicating which of the parameters are active in the order \code{(theta, x)}.  If both \code{vars} and \code{fixed.params}/\code{nmiss} are given, an error is thrown if their active sets differ.
#' @export
sde.active.vars <- function(model, vars, fixed.params, nmiss, debug = FALSE) {
  param.names <- model$param.names
  nparams <- model$nparams
  data.names <- model$data.names
  ndims <- model$ndims
  var.names <- c(param.names, data.names)
  var.mode <- !missing(vars)
  fix.mode <- !missing(fixed.params) && !missing(nmiss)
  if(var.mode) {
    if(!all(vars %in% var.names)) {
      stop("vars not in param.names or data.names.")
    }
    if(anyDuplicated(vars)) {
      stop("vars must be unique.")
    }
  }
  if(debug) browser()
  if(fix.mode) {
    vars2 <- c(param.names[!fixed.params], data.names[ndims:1 <= nmiss])
    if(!var.mode) vars <- vars2
  }
  var.id <- var.names %in% vars
  if(var.mode && fix.mode) {
    if(!all(var.id == (var.names %in% vars2))) {
      stop("vars and (fixed.params,nmiss) specify different active sets. ")
    }
  }
  names(var.id) <- var.names
  var.id
}
