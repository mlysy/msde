.onLoad <- function(libname, pkgname) {
  assign(".msde_src_path",
         file.path(libname, pkgname, "tools"),
         envir = parent.env(environment()))
  assign(".msde_examples_path",
         file.path(libname, pkgname, "examples"),
         envir = parent.env(environment()))
}
