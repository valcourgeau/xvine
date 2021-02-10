test_that("models/fit_markov_svines/", {
  n <- 1000
  ncol <- 3
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- matrix(evir::rgpd(n = n*ncol, xi = xi, beta = beta), ncol = ncol, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, ncol))
  svine_fit <- fit_markov_svines(data_it$data, k.markov=1, family_set="archimedean", selcrit="mbicv")
})


test_that("models/fit_markov_rvines/", {
  n <- 1000
  ncol <- 3
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- matrix(evir::rgpd(n = n*ncol, xi = xi, beta = beta), ncol = ncol, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, ncol))
  rvine_fit <- fit_markov_rvines(
    data_it$data, k.markov=5,
    copula_controls = list(family_set="archimedean", selcrit="mbicv")
  )
})
