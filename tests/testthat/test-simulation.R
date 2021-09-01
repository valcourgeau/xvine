test_that("postprocessing/model_simulation/fit_markov_svines/", {
  n <- 5000
  n_col <- 3
  k.m <- 2
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 10

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_svines(
    data = data_it$data, k.markov=k.m, family_set="archimedean", selcrit="mbicv"
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
  testthat::expect_equal(as.vector(dim(v_sims)), c(k.m, n_col, n_sims))
  testthat::expect_lte(max(v_sims), 1.0)
  testthat::expect_gte(max(v_sims), 0.0)
})

test_that("postprocessing/model_simulation/fit_markov_rvines/", {
  n <- 1000
  n_col <- 3
  k.m <- 2
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 20

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  rvine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=rvine_fit)
  testthat::expect_equal(as.vector(dim(v_sims)), c(k.m, n_col, n_sims))
  testthat::expect_lte(max(v_sims), 1.0)
  testthat::expect_gte(max(v_sims), 0.0)
})

test_that("postprocessing/model_simulation/shape_svines_rvines/", {
  n <- 5000
  n_col <- 3
  k.m <- 2
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 20

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_svines(
    data_it$data, k.markov=2,
    family_set="archimedean", selcrit="mbicv",
  )
  rvine_fit <- fit_markov_rvines(
    data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )

  r_sims <- model_simulation(n=n_sims, model=rvine_fit)
  s_sims <- model_simulation(n=n_sims, model=svine_fit)
  testthat::expect_equal(
    as.vector(dim(r_sims)), as.vector(dim(s_sims))
  )
})

test_that("postprocessing/model_simulation/fit_markov_conditional_rvines/", {
  n <- 1000
  n_col <- 3
  k.m <- 2
  quant <- .5
  xi <- 0.1
  beta <- 1.

  n_sims <- 20

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  cond_rvine_fit <- fit_markov_conditional_rvines(
    data = data_it$data, k.markov=k.m, col_source = 1, u0 = quant,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=cond_rvine_fit)
  testthat::expect_equal(as.vector(dim(v_sims)), c(k.m, n_col, n_sims))
  testthat::expect_lte(max(v_sims), 1.0)
  testthat::expect_gte(max(v_sims), 0.0)
})

test_that("postprocessing/model_simulation/fit_markov_conditional_rvines/above_vine_is_above_u0", {
  n <- 1000
  n_col <- 3
  k.m <- 2
  quant <- .5
  xi <- 0.1
  beta <- 1.

  col_fixed <- 1

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  cond_rvine_fit <- fit_markov_conditional_rvines(
    data = data_it$data, k.markov=k.m, col_source = col_fixed, u0 = quant,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=cond_rvine_fit)
  testthat::expect_equal(sum(v_sims[1,col_fixed,] > quant), n_sims/2, tolerance = 5)
  testthat::expect_equal(sum(v_sims[1,col_fixed,] <= quant), n_sims/2, tolerance = 5)
})

test_that("postprocessing/model_simulation/fit_markov_timewise_rvines/", {
  n <- 1000
  n_col <- 3
  k.m <- 2
  quant <- .5
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  timewise_rvine_fit <- fit_markov_timewise_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=timewise_rvine_fit)
  testthat::expect_equal(dim(v_sims), c(k.m, n_col, n_sims))
})

test_that("postprocessing/model_simulation/fit_markov_conditional_timewise_rvines/", {
  n <- 1000
  n_col <- 3
  k.m <- 2
  quant <- .5
  col_fixed <- 2
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  cond_timewise_rvine_fit <- fit_markov_conditional_timewise_rvines(
    data = data_it$data, k.markov=k.m, col_source = col_fixed, u0 = quant,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=cond_timewise_rvine_fit)
  testthat::expect_equal(dim(v_sims), c(k.m, n_col, n_sims))
  n_above <- sum(apply(v_sims, MARGIN = 3, FUN = function(x){x[1,col_fixed] > quant}))
  testthat::expect_equal(n_above, n_sims/2, tolerance = 3)
})
