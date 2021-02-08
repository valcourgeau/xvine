test_that("preprocessing/uniform_ecdf/", {
  n <- 1000
  n_col <- 10
  x <- matrix(runif(n = n, min = 0, max = 1), ncol = n_col)
  x_ecdf <- uniform_ecdf(x)
  testthat::expect_equal(dim(x_ecdf), c(n/n_col, n_col))
})
