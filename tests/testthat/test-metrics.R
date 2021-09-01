test_that("metrics/binning/", {
  n <- 1000
  x <- runif(n)
  x_b <- binning(x)
  testthat::expect_equal(sum(x_b), n)
})

test_that("metrics/binning/dims", {
  n <- 1000
  x <- runif(n)
  y <- runif(2*n)
  x_b <- binning(x)
  y_b <- binning(y)

  testthat::expect_equal(sum(x_b), n)
  testthat::expect_equal(sum(y_b), 2*n)
})

test_that("metrics/coercive_binning/dims", {
  n <- 1000
  x <- runif(n)
  y <- runif(2*n)
  b <- coercive_binning(x, y)

  testthat::expect_equal(dim(b), c(2, n/10))
})

test_that("metrics/distr_distance/kl", {
  n <- 10000
  x <- runif(n)
  y <- runif(n)
  testthat::expect_equal(
    distr_distance(x, y, method='kullback-leibler', per_bin=100),
    0.01,
    tolerance = 0.02
  )
})

test_that("metrics/distr_distance/euclidean", {
  n <- 10000
  x <- runif(n)
  y <- runif(n)
  testthat::expect_equal(
    distr_distance(x, y, method='euclidean', per_bin=100),
    0.01,
    tolerance = 0.02
  )
})

test_that("metrics/apply_distr_distance/dims", {
  n <- 50000
  d <- 4
  x <- matrix(runif(n), nrow = d)
  y <- matrix(runif(n), nrow = d)
  distr_dist <- apply_distr_distance(x, y, method='kullback-leibler', per_bin=100)

  testthat::expect_equal(length(distr_dist), d)
  testthat::expect_equal(
    distr_dist, rep(0.01, d), tolerance = 0.02
  )
})

test_that("metrics/apply_distr_distance/dims_3d", {
  n <- 10000
  d <- 4
  k <- 3
  set.seed(52)
  x <- array(runif(n), c(k, d, n))
  y <- array(runif(n), c(k, d, n))
  distr_dist <- apply_distr_distance(x, y, method='kullback-leibler', per_bin=100)
  testthat::expect_equal(dim(distr_dist), c(k, d))
})


test_that("metrics/apply_distr_distance/transpose", {
  n <- 50000
  d <- 4
  x <- matrix(runif(n), ncol = d)
  y <- matrix(runif(n+100), ncol = d)
  d_dist <- apply_distr_distance(
    x, y, method='kullback-leibler', per_bin=100
  )

  testthat::expect_equal(length(d_dist), d)
  testthat::expect_equal(
    d_dist, rep(0.01, d), tolerance = 0.02
  )
})
