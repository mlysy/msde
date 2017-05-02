.onLoad <-
function(libname, pkgname) {
  assign(".msdeCppPath", file.path(libname, pkgname, "cppTemplates"), envir = parent.env(environment()))
}
