# The setting of R_TESTS exists to work around an R bug. See
# https://github.com/hadley/testthat/issues/144
# We should remove it when the issue is resolved.
Sys.unsetenv("R_TESTS")
library(testthat)
library(lemad2)

test_check("lemad")
