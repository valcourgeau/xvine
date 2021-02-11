model_simulation <- function(n, model, ...){
  n_cores <- parallel::detectCores()
  if(inherits(model, "svine_dist")){
    k.markov <- model$copula$p
    vine_sim_vals <- NULL
    while(is.null(vine_sim_vals)){
      tryCatch(
        {vine_sim_vals <- svines::svine_sim(n = k.markov, rep = n, model = model, cores = n_cores, ...)}
      )
    }

    return(vine_sim_vals)
  }else{
    if(inherits(model, "vine")){
      vine_sim_vals <- rvinecopulib::rvine(n = n, vine = model, cores = n_cores, ...)

      # process into k.markov * d * n
      k.markov <- dim(vine_sim_vals)[2] / model$d
      vine_sim_vals <- t(vine_sim_vals) # dims: (d*k.markov) * n
      vine_sim_vals <- array(data = vine_sim_vals, c(model$d, k.markov, n))
      vine_sim_vals <- aperm(vine_sim_vals, c(2, 1, 3))
      return(vine_sim_vals)
    }
  }
}
