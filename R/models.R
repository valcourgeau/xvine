
fit_markov_svines <- function(data, k.markov, ...){
  n_cores <- parallel::detectCores()
  return(svines::svine(data,p = k.markov, cores = n_cores, ...))
}

fit_markov_rvines <- function(data, k.markov, ...){
  n_cores <- parallel::detectCores()
  dt_stack <- build_stack(data, k.markov)
  return(rvinecopulib::vine(dt_stack, cores = n_cores, ...))
}

fit_markov_conditional_rvines <- function(data, k.markov, col, u0, ...){
  n_cores <- parallel::detectCores()
  dt_stack <- build_conditional_stack(data, k.markov, col, u0)
  return(rvinecopulib::vine(dt_stack, cores = n_cores, ...))
}
