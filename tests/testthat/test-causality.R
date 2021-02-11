test_that("causality/fit_markov_svines/", {
  set.seed(42)

  n <- 5000
  n_col <- 3
  k.m <- 5
  quant <- .95
  xi <- 0.1
  beta <- 1.

  n_sims <- 100

  data <- matrix(evir::rgpd(n = n*ncol, xi = xi, beta = beta), ncol = ncol, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, ncol))
  svine_fit <- fit_markov_svines(data_it$data, k.markov=k.m, family_set="archimedean", selcrit="mbicv")
  v_sims <- model_simulation(n=n_sims, model=svine_fit, qrng=F)
  for(i in 1:n_col){
    p_sc_st <- proba_single_conditional_single_target(
      v_sims, col_source=i, col_target=i, u0_source=.1, u0_target=.1
    )

    w_t <- rep(0, n_col)
    w_t[i] <- 1
    p_sc_ct <- proba_single_conditional_cross_target(
      v_sims, col_source=i, w_target=w_t, u0_source=.1, u0_target=.1
    )
    testthat::expect_equal(p_sc_st[['factual']], p_sc_ct[['factual']])
    testthat::expect_equal(p_sc_st[['counterfactual']], p_sc_ct[['counterfactual']])
  }
})

