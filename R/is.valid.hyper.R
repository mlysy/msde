is.valid.hyper <- function(phi) {
  if(!is.list(phi)) stop("phi must be a list.")
  is.valid.phi <- sapply(phi, function(x) {
    is.null(x) || (is.double(x) & is.vector(x))
  })
  is.valid.phi <- all(is.valid.phi)
  if(!is.valid.phi) {
    stop("Each element of phi must be NULL or a double vector.")
  }
  is.valid.phi
}
