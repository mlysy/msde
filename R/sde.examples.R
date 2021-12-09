#' Example SDE models.
#'
#' Provides sample \code{C++} code for several SDE models.
#' @param model Character string giving the name of a sample model.  Possible values are: \code{hest}, \code{pgnet}, \code{lotvol}, \code{biou}, \code{eou}.  See Details.
#' @param file.only If \code{TRUE} returns only the path to the header file containing the \code{sdeModel} object implementation.
#' @return An \code{sde.model} object, or the path to the C++ model header file.
#' @details All pre-compiled models are with the default prior and with \code{OpenMP} disabled.  A full description of the example models can be found in the package vignette; to view it run \code{vignette("msde-exmodels")}.
#' @seealso \code{\link{sde.make.model}} for \code{sde.model} objects, \code{\link{mvn.hyper.check}} for specification of the default prior.
#' @examples
#' # Heston's model
#' hmod <- sde.examples("hest") # load pre-compiled model
#'
#' # inspect model's C++ code
#' hfile <- sde.examples("hest", file.only = TRUE)
#' cat(readLines(hfile), sep = "\n")
#'
#' \dontrun{
#' # compile it from scratch
#' param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
#' data.names <- c("X", "Z")
#' hmod <- sde.make.model(ModelFile = hfile,
#'                        param.names = param.names,
#'                        data.names = data.names)
#' }
#' @export
sde.examples <- function(model = c("hest", "pgnet", "lotvol", "biou", "eou"),
                         file.only = FALSE) {
  model <- match.arg(model)
  if(model == "hest") {
    ModelFile <- file.path(.msde_include_path, "hestModel.h")
    param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
    data.names <- c("X", "Z")
    cobj <- new(msde_hestModel)
  } else if(model == "pgnet") {
    ModelFile <- file.path(.msde_include_path, "pgnetModel.h")
    param.names <- paste0("theta", 1:8)
    data.names <- c("R", "P", "Q", "D")
    cobj <- new(msde_pgnetModel)
  } else if(model == "biou") {
    ModelFile <- file.path(.msde_include_path, "biouModel.h")
    param.names <- c("Gamma11", "Gamma21", "Gamma12", "Gamma22",
                     "Lambda1", "Lambda2",
                     "Psi11", "Psi12", "Psi22")
    data.names <- c("Y1","Y2")
    cobj <- new(msde_biouModel)
  } else if(model == "lotvol") {
    ModelFile <- file.path(.msde_include_path, "lotvolModel.h")
    param.names <- c("alpha", "beta", "gamma")
    data.names <- c("H", "L")
    cobj <- new(msde_lotvolModel)
  } else if(model == "eou") {
    ModelFile <- file.path(.msde_include_path, "eouModel.h")
    param.names <- c("alpha", "gamma", "eta", "sigma", "rho")
    data.names <- c("X", "V")
    cobj <- new(msde_eouModel)
  }
  if(!file.only) {
    sde.model <- list(cobj = cobj,
                      ndims = length(data.names), nparams = length(param.names),
                      data.names = data.names, param.names = param.names,
                      hyper.check = mvn.hyper.check, omp = FALSE)
    class(sde.model) <- "sde.model"
  } else sde.model <- ModelFile
  sde.model
}
