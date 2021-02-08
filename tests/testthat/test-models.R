test_that("models/fit_markov_svines/", {
  n <- 1000
  ncol <- 10
  quant <- .75
  xi <- 0.1
  beta <- 1.

  data <- matrix(evir::rgpd(n = n*ncol, xi = xi, beta = beta), ncol = ncol, nrow = n)
  data_it <- apply_integral_transform(data, u0s=rep(quant, ncol))
  svine_fit <- fit_markov_svines(data_it$data, k.markov=1, family_set="archimedean", cores=8)
  print(summary(svine_fit))
})
