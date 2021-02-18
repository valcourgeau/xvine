test_that("preprocessing/uniform_ecdf/", {
  n <- 1000
  n_col <- 10
  x <- matrix(runif(n = n, min = 0, max = 1), ncol = n_col)
  x_ecdf <- uniform_ecdf(x)
  testthat::expect_equal(dim(x_ecdf), c(n/n_col, n_col))
})

test_that("preprocessing/partial_ecdf/", {
  n <- 1000
  quant <- .75
  x <- runif(n = n, min = 0, max = 1)

  x_partial_ecdf <- partial_ecdf(x, .75)[['data']]

  testthat::expect_equal(length(x_partial_ecdf), n * quant)
  testthat::expect_lte(max(x_partial_ecdf), 1.0)
  testthat::expect_gte(max(x_partial_ecdf), 0.0)
})

test_that("preprocessing/partial_gpd/", {
  n <- 1000
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- evir::rgpd(n = n, xi = xi, beta = beta)

  x_partial_ecdf <- partial_gpd(x, .75)[['data']]

  testthat::expect_equal(length(x_partial_ecdf), n * (1-quant))
  testthat::expect_lte(max(x_partial_ecdf), 1.0)
  testthat::expect_gte(max(x_partial_ecdf), 0.0)
})

test_that("preprocessing/integral_transform/", {
  n <- 10000
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- evir::rgpd(n = n, xi = xi, beta = beta)

  x_it <- integral_transform(x, u0=quant)[['data']]
  testthat::expect_equal(length(x_it), n)
  testthat::expect_lte(max(x_it), 1.0)
  testthat::expect_gte(max(x_it), 0.0)
  testthat::expect_gte(ks.test(x_it, "punif", 0, 1)$p.value, .05)
})

test_that("preprocessing/apply_integral_transform/", {
  n <- 1000
  ncol <- 10
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- matrix(evir::rgpd(n = n*ncol, xi = xi, beta = beta), ncol = ncol, nrow = n)

  u0s <- rep(quant, ncol)
  x_it <- apply_integral_transform(x, u0s=u0s)
  testthat::expect_equal(dim(x_it$data), c(n, ncol))
  testthat::expect_lte(max(x_it$data), 1.0)
  testthat::expect_gte(max(x_it$data), 0.0)
})

test_that("preprocessing/apply_integral_transform/three_dims", {
  n <- 1000
  n_col <- 10
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- array(evir::rgpd(n = n*n_col*n_col, xi = xi, beta = beta), c(n_col, n_col, n))

  u0s <- rep(quant, n_col)
  x_it <- apply_integral_transform(x, u0s=u0s)
  testthat::expect_equal(dim(x_it$data), c(n_col, n_col, n))
  testthat::expect_lte(max(x_it$data), 1.0)
  testthat::expect_gte(max(x_it$data), 0.0)
})

test_that("preprocessing/reverse_integral_transform/", {
  n <- 10000
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- evir::rgpd(n = n, xi = xi, beta = beta)

  x_it <- integral_transform(x, u0=quant)[['data']]
  x_rit <- reverse_integral_transform(x=x_it, x_source=x, u0=quant, shape = xi, scale = beta)[['data']]
  testthat::expect_equal(length(x_rit), n)
  testthat::expect_equal(max(x_rit)/max(x), 1, tolerance=0.30)
  testthat::expect_equal(min(x_rit), min(x), tolerance=1e-3)
})

test_that("preprocessing/apply_reverse_integral_transform/", {
  n <- 1000
  n_col <- 10
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)

  u0s <- rep(quant, n_col)
  x_it <- apply_integral_transform(x, u0s=u0s)$data
  x_rit <- apply_reverse_integral_transform(
    x_it, data_source=x, u0s = rep(quant, n_col), shapes = rep(xi, n_col), scales = rep(beta, n_col)
  )
  x_rit <- x_rit$data
  testthat::expect_equal(dim(x_rit), c(n, n_col))
  testthat::expect_equal(max(x_rit)/max(x), 1, tolerance=0.30)
  testthat::expect_equal(min(x_rit), min(x), tolerance=5e-3)
})

test_that("preprocessing/apply_reverse_integral_transform/three_dims", {
  n <- 10000
  n_col <- 10
  quant <- .75
  xi <- 0.1
  beta <- 1.
  x <- array(evir::rgpd(n = n*n_col*n_col, xi = xi, beta = beta), c(n_col, n_col, n))

  u0s <- rep(quant, n_col)
  x_it <- apply_integral_transform(x, u0s=u0s)$data
  x_rit <- apply_reverse_integral_transform(
    x_it, data_source=x, u0s = rep(quant, n_col), shapes = rep(xi, n_col), scales = rep(beta, n_col)
  )
  x_rit <- x_rit$data
  testthat::expect_equal(dim(x_rit), c(n_col, n_col, n))
  testthat::expect_equal(max(x_rit)/max(x), 1, tolerance=0.30)
  testthat::expect_equal(min(x_rit), min(x), tolerance=5e-3)
})
