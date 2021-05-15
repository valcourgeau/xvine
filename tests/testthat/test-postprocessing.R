test_that("postprocessing/proba_single_target/dim_rvines/", {
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
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=T)
  p_sc_st <- proba_single_target(
    v_sims, col_source=1, col_target=1, u0_source=.1, u0_target=.1
  )
  testthat::expect_equal(length(p_sc_st[['factual']]), k.m-1)
  testthat::expect_equal(length(p_sc_st[['counterfactual']]), k.m-1)
})

test_that("postprocessing/proba_single_target/reverse/", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .1
  xi <- 0.1
  beta <- 1.

  n_sims <- 1000

  u0s_fixed <- rep(quant, n_col)
  u0_rev <- evir::qgpd(quant, xi = xi, beta = beta)

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=u0s_fixed)
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
  v_sims_reverse <- apply_reverse_integral_transform(
    data_unif=v_sims, data_source=data, u0s=u0s_fixed, shapes=data_it$par.ests[,1], scales=data_it$par.ests[,2]
  )
  testthat::expect_equal(dim(v_sims), dim(v_sims_reverse$data))

  p_sc_st <- proba_single_target(
    v_sims, col_source=1, col_target=1, u0_source=quant, u0_target=quant
  )

  p_sc_st_rev <- proba_single_target(
    v_sims_reverse$data, col_source=1, col_target=1, u0_source=u0_rev, u0_target=u0_rev
  )
  testthat::expect_equal(length(p_sc_st[['factual']]), k.m-1)
  testthat::expect_equal(length(p_sc_st[['counterfactual']]), k.m-1)
  testthat::expect_equal(p_sc_st[['factual']], p_sc_st_rev[['factual']], tolerance = 0.1)
  testthat::expect_equal(p_sc_st[['counterfactual']], p_sc_st_rev[['counterfactual']], tolerance = 1e-2)
})

test_that("postprocessing/proba_cross_target/dim_rvines/", {
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
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=T)

  w_t <- rep(1, n_col) / n_col
  p_sc_ct <- proba_cross_target(
    v_sims, col_source=1, w_target=w_t, u0_source=.1, u0_target=.1
  )
  testthat::expect_equal(length(p_sc_ct[['factual']]), k.m-1)
  testthat::expect_equal(length(p_sc_ct[['counterfactual']]), k.m-1)
})

test_that("postprocessing/proba_single_target/regression_with_cross_target/", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- 0.1
  xi <- 0.1
  beta <- 1.

  n_sims <- 1000

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
  for(i in 1:n_col){
    p_sc_st <- proba_single_target(
      v_sims, col_source=i, col_target=i, u0_source=.1, u0_target=.1
    )

    w_t <- rep(0, n_col)
    w_t[i] <- 1
    p_sc_ct <- proba_cross_target(
      v_sims, col_source=i, w_target=w_t, u0_source=.1, u0_target=.1
    )
    testthat::expect_equal(p_sc_st[['factual']], p_sc_ct[['factual']], tolerance = 0.1)
    testthat::expect_equal(p_sc_st[['counterfactual']], p_sc_ct[['counterfactual']], tolerance = 0.1)
  }

})

test_that("postprocessing/proba_time_target/dim_rvines/", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 1000

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=T)

  w_t <- rep(1, k.m - 1) / (k.m-1)
  p_sc_tt <- proba_time_target(
    v_sims, col_source=1, w_target=w_t, u0_source=.1, u0_target=.9
  )
  testthat::expect_equal(length(p_sc_tt[['factual']]), n_col)
  testthat::expect_equal(length(p_sc_tt[['counterfactual']]), n_col)
})

test_that("postprocessing/proba_all_target/dim_rvines/", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 1000

  data <- matrix(evir::rgpd(n = n*n_col, xi = xi, beta = beta), ncol = n_col, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, n_col))
  svine_fit <- fit_markov_rvines(
    data = data_it$data, k.markov=k.m,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=T)

  w_t <- matrix(1, nrow = k.m - 1, ncol=n_col) / (n_col * (k.m-1))
  p_sc_at <- proba_all_target(
    v_sims, col_source=1, w_target=w_t, u0_source=.1, u0_target=.9
  )
  testthat::expect_equal(length(p_sc_at[['factual']]), 1)
  testthat::expect_equal(length(p_sc_at[['counterfactual']]), 1)
})
