library(testthat)
library(msdeHeaders)

Sys.unsetenv("R_TESTS")
test_check("msdeHeaders")
