
fit_markov_svines <- function(data, k.markov, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)

  n_cores <- parallel::detectCores()
  return(svines::svine(data = data, p = k.markov, cores = n_cores, ...))
}

fit_markov_rvines <- function(data, k.markov, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)

  n_cores <- parallel::detectCores()
  d <- ncol(data)
  dt_stack <- build_stack(data, k.markov)
  vine <- rvinecopulib::vine(
    dt_stack, cores = n_cores,
    margins_controls=list(
      xmin=rep(0, d),
      xmax=rep(1, d)
    ),
    ...
  )
  vine$d <- ncol(data)
  return(vine)
}

fit_markov_conditional_rvines <- function(data, k.markov, col, u0, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)

  n_cores <- parallel::detectCores()
  dt_stack <- build_conditional_stack(data, k.markov, col, u0)
  vine <- rvinecopulib::vine(
    dt_stack, cores = n_cores,
    margins_controls=list(
      xmin=rep(0, d),
      xmax=rep(1, d)
    ),
    ...
  )
  vine$d <- ncol(data)
  return(vine)
}
