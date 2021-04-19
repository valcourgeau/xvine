test_that("threshold/forward_stop/", {
  set.seed(42)
  n <- 50000
  m <- 9
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- evir::rgpd(n = n, xi = xi, beta = beta)
  fs <- forward_stop(
    p_values = seq_len(m)/10
  )
  testthat::expect_equal(length(fs$idx), 1)
  testthat::expect_equal(length(fs$stopping_values), m)
})

test_that("threshold/gpd_goodness/", {
  set.seed(42)
  n <- 50000
  m <- 9
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- evir::rgpd(n = n, xi = xi, beta = beta)
  gof_results <- gpd_goodness(
    data = data,
    u0s = seq_len(m)/10,
    gof_test = eva::gpdAd
  )
  testthat::expect_equal(length(gof_results$p.values), m)
  testthat::expect_equal(length(gof_results$statistics), m)
  testthat::expect_equal(length(gof_results$optimal_idx), 1)
})

test_that("threshold/apply_gpd_goodness/", {
  set.seed(42)
  n <- 10000
  n_col <- 5
  m <- 9
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col)
  gof_results <- apply_gpd_goodness(
    data = data,
    u0s = seq_len(m)/10,
    gof_test = eva::gpdAd
  )
  testthat::expect_equal(dim(gof_results$p.values), c(m, n_col))
  testthat::expect_equal(gof_results$optimal_idx, rep(1, n_col))
  testthat::expect_lte(max(gof_results$p.values, na.rm = T), 1.)
  testthat::expect_gte(min(gof_results$p.values, na.rm = T), 0.)
})

