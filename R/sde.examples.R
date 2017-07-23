#' Example SDE models.
#'
#' Provides sample \code{C++} code for several SDE models.
#' @param model Character string giving the name of a sample model.  Possible values are: \code{hest}, \code{pgnet}, \code{lotvol}, \code{biou}.  See Details.
#' @param ModelFile.only If \code{TRUE} returns only the path to the header file containing the \code{sdeModel} object implementation.
#' @return An \code{sde.model} object (see \code{\link{sde.make.model}}).
#' @details A full description of the sample models can be found in the package vignette; to view it run \code{vignette("msde-exmodels")}.
#' @examples
#' \donttest{
#' # Heston's model
#' hex <- example.models("hest")
#' cat(readLines(hex$ModelFile), sep = "\n") # view model's C++ code
#'
#' # compile the model
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#' }
#' @export
sde.examples <- function(model = c("hest", "pgnet", "lotvol", "biou"),
                         ModelFile.only = FALSE) {
  model <- match.arg(model)
  if(model == "hest") {
    ModelFile <- file.path(.msde_examples_path, "hestModel.h")
    param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
    data.names <- c("X", "Z")
    sptr <- .hest_MakeModel()
  } else if(model == "pgnet") {
    ModelFile <- file.path(.msde_examples_path, "pgnetModel.h")
    param.names <- paste0("theta", 1:8)
    data.names <- c("R", "P", "Q", "D")
    sptr <- .pgnet_MakeModel()
  } else if(model == "biou") {
    ModelFile <- file.path(.msde_examples_path, "biouModel.h")
    param.names <- c("Gamma11", "Gamma21", "Gamma12", "Gamma22",
                     "Lambda1", "Lambda2",
                     "Psi11", "Psi12", "Psi22")
    data.names <- c("X1","X2")
    sptr <- .biou_MakeModel()
  } else if(model == "lotvol") {
    ModelFile <- file.path(.msde_examples_path, "lotvolModel.h")
    param.names <- c("alpha", "beta", "gamma")
    data.names <- c("H", "L")
    sptr <- .lotvol_MakeModel()
  }
  if(!ModelFile.only) {
    sde.model <- list(ptr = sptr,
                      ndims = length(data.names), nparams = length(param.names),
                      data.names = data.names, param.names = param.names,
                      hyper.check = mvn.hyper.check, omp = FALSE)
    class(sde.model) <- "sde.model"
  } else sde.model <- ModelFile
  sde.model
}
