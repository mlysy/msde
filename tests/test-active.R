#--- test active model variables -----------------------------------------------

require(msdeHeaders)
#require(mvtnorm)

param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
data.names <- c("X", "Z")
hmod <- sde.make.model(ModelFile = "hestModel.h",
                       param.names = param.names,
                       data.names = data.names)
ndims <- hmod$ndims
nparams <- hmod$nparams

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

vars <- c("alpha", "sigma", "beta")
fixed.params <- c(FALSE, TRUE, FALSE, FALSE, TRUE)
nmiss <- 0
sde.active.vars(model = hmod, vars = vars,
                fixed.params = fixed.params, nmiss = nmiss)

