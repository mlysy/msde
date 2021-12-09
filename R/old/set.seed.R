# reset fast norm_rand parity.
#' @keywords internal
#' @export
set.seed <- function(seed, kind = NULL, normal.kind = NULL) {
  if(is.null(kind)) kind <- RNGkind()[1]
  if(is.null(normal.kind)) normal.kind <- RNGkind()[2]
  if(normal.kind == "user-supplied") {
    base::set.seed(1, kind = "user-supplied", normal.kind = "user-supplied")
  }
  base::set.seed(seed, kind, normal.kind)
}
