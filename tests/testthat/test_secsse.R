context("test_secsse")

test_that("secsse gives the same result as hisse", {
  secsse_test_hisse()
})

test_that("secsse gives the same result as GeoSSE", {
  secsse_test_geosse()  
})

test_that("trying a short ML search: secsse_ml & Parallel procedure", {
 secsse_test_ml()
})  

test_that("trying a short ML search: secsse_ml_func_def_pars", {
  secsse_test_ml2()
})

test_that("trying a short ML search: cla_secsse", {
  secsse_test_ml3()
})
