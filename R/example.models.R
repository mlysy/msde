#' Example SDE models
#'
#' @param model Character string giving the name of a sample model.  Possible values are: \code{hest}, \code{pgnet}.  See Details.
#' @return A list with elements
#' \describe{
#' \item{\code{ModelFile}}{The full path to the example's header file}
#' \item{\code{param.names}}{The names of the model parameters}
#' \item{\code{data.names}}{The names of the SDE model components}
#' }
#' @examples
#' \donttest{
#' # compile Heston's model
#' hex <- example.models("hest")
#' hmod <- sde.make.model(ModelFile = hex$ModelFile,
#'                        param.names = hex$param.names,
#'                        data.names = hex$data.names)
#' }
#' @export
example.models <- function(model = c("hest", "pgnet")) {
  model <- match.arg(model)
  if(model == "hest") {
    ModelFile <- file.path(.msdeCppPath, "hestModel.h")
    param.names <- c("alpha", "gamma", "beta", "sigma", "rho")
    data.names <- c("X", "Z")
  } else if(model == "pgnet") {
    ModelFile <- file.path(.msdeCppPath, "pgnetModel.h")
    param.names <- paste0("gamma", 1:8)
    data.names <- c("R", "P", "Q", "D")
  }
  list(ModelFile = ModelFile, param.names = param.names,
       data.names = data.names)
}
