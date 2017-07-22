library(testthat)
library(msde2)

Sys.unsetenv("R_TESTS")
test_check("msde")
