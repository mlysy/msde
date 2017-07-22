#' R constructor
#'
#' @importFrom Rcpp sourceCpp
#' @export
Robj.make <- function(type = c("A", "B")) {
  type <- match.arg(type)
  fname <- paste0("Derived", type, ".h")
  file.copy(from = fname, to = file.path(tempdir(), "Derived.h"),
            overwrite = TRUE, copy.date = TRUE)
  fname <- paste0("Derived", type, ".cpp")
  file.copy(from = file.path(.msde_src_path, "Derived.cpp"),
            to = file.path(tempdir(), fname),
            overwrite = TRUE, copy.date = TRUE)
  sourceCpp(file.path(tempdir(), fname), env = environment())
  obj <- list(ptr = construct())
  class(obj) <- "Robj"
  obj
}
