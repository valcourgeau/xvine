test_that("causality/wrapper_pn_ct/numeric_results", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
  wrap_ct <- wrapper_pn_ct(v_sims, col_source=1, times=3, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_ct(rep(1, n_col))))

  wrap_ct <- wrapper_pn_ct(v_sims, col_source=1, times=3:5, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_ct(rep(1, n_col))))

  # print(microbenchmark::microbenchmark(wrap_ct(rep(1, n_col))))
})

test_that("causality/wrapper_pn_tt/numeric_results", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)

  w_t <- rep(1,  k.m-1)
  wrap_tt <- wrapper_pn_tt(v_sims, col_target=3, col_source=1, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_tt(w_t)))
  testthat::expect_equal(length(wrap_tt(w_t)), 1)

  wrap_tt <- wrapper_pn_tt(v_sims, col_target=2:3, col_source=1, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_tt(w_t)))
  testthat::expect_equal(length(wrap_tt(w_t)), 1)

  # print(microbenchmark::microbenchmark(wrap_tt(rep(1, n_col))))
})

test_that("causality/wrapper_pn_all/numeric_results", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)

  w_t <- rep(1,  n_col*(k.m-1))
  wrap_all <- wrapper_pn_all(v_sims, col_source=1, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_all(w_t)))
  testthat::expect_equal(length(wrap_all(w_t)), 1)

  wrap_tt <- wrapper_pn_all(v_sims, col_source=1, u0_source=.1, u0_target=.1)
  testthat::expect_true(is.numeric(wrap_all(w_t)))
  testthat::expect_equal(length(wrap_all(w_t)), 1)

  # print(microbenchmark::microbenchmark(wrap_tt(rep(1, n_col))))
})

test_that("causality/wrapper_pn_all/optim", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)

  wrap_all <- wrapper_pn_all(v_sims, col_source=1, u0_source=.1, u0_target=.1)
  sc_all_optim <- optim(
    par = rep(1, n_col * (k.m - 1)),
    fn = function(w){-wrap_all(w)}
  )
  testthat::expect_equal(sc_all_optim$convergence, 0)
})

test_that("causality/wrapper_pn_ct/optim", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)

  wrap_ct <- wrapper_pn_ct(v_sims, times = 1, col_source=1, u0_source=.1, u0_target=.1)
  sc_all_optim <- optim(
    par = rep(1, n_col),
    fn = function(w){-wrap_ct(w)}
  )
  testthat::expect_equal(sc_all_optim$convergence, 0)
})

# test_that("causality/cross/optim", {
#   set.seed(42)
#
#   n <- 5000
#   n_col <- 3
#   k.m <- 5
#   quant <- .95
#   xi <- 0.1
#   beta <- 1.
#
#   n_sims <- 100
#
#   data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
#   data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
#   svine_fit <- fit_markov_rvines(
#     data = data_it$data, k.markov=k.m,
#     copula_controls = list(family_set="archimedean", selcrit="mbicv")
#   )
#   v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
#
#   wrap_ct <- wrapper_pn_ct(v_sims, times = 1, col_source=1, u0_source=.1, u0_target=.1)
#   sc_all_optim <- optim(
#     par = rep(1, n_col),
#     fn = function(w){-wrap_ct(w)}
#   )
#   testthat::expect_equal(sc_all_optim$convergence, 0)
# })


