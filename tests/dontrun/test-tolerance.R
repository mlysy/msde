# check that tests have sufficient tolerance.

require(msde)

system.time({
  ans <- sapply(1:20, function(ii) {
    message("test ", ii)
    testthat::test_package("msde", reporter = "progress")
  })
})

# now check results
res <- apply(ans, 2, function(aa) {
  sapply(unlist(sapply(aa, function(ls) ls$results), recursive = FALSE),
         function(ls) ls$message)
})
all(res %in% c("success", "PF test too long"))
