model_simulation <- function(n, model, ...){
  n_cores <- parallel::detectCores()
  if(inherits(model, "svine_dist")){
    k.markov <- model$copula$p
    vine_sim_vals <- NULL
    while(is.null(vine_sim_vals)){
      vine_sim_vals <-tryCatch(
        { svines::svine_sim(n = k.markov, rep = n, model = model, cores = n_cores, ...)},
        error=function(cond){return(NULL)}
      )
    }

    return(vine_sim_vals)
  }

  if(inherits(model, "vine")){
    vine_sim_vals <- NULL
    while(is.null(vine_sim_vals)){
      vine_sim_vals <- tryCatch(
        { rvinecopulib::rvine(n = n, vine = model, cores = n_cores, ...)},
        error=function(cond){return(NULL)}
      )
    }

    # process into k.markov * d * n
    k.markov <- dim(vine_sim_vals)[2] / model$d
    vine_sim_vals <- t(vine_sim_vals) # dims: (d*k.markov) * n
    vine_sim_vals <- array(data = vine_sim_vals, c(model$d, k.markov, n))
    vine_sim_vals <- aperm(vine_sim_vals, c(2, 1, 3))
    return(vine_sim_vals)
  }

  if(inherits(model, "cond_vine")){
    n_half <- round(n/2)
    vine_sim_vals_above <- model_simulation(n_half, model = model$above, ...) # (k, d, n_half)
    vine_sim_vals_below <- model_simulation(n-n_half, model = model$below, ...) # (k, d, n-n_half)
    vine_sim_vals <- abind::abind(vine_sim_vals_above, vine_sim_vals_below, along=3) # (k, d, n)
    return(vine_sim_vals)
  }

  if(inherits(model, "timewise_vine")){
    # create base vector
    k.m <- length(model$timewise)

    vine_sim_first <- rvinecopulib::rvine(
      n = n, vine = model$timewise[[1]], cores = n_cores
    ) # this is the simulation of the first vine
    vine_sim_base <- vine_sim_first[, seq_len(model$d)] # keep the start values x_t
    vine_sim_first <- vine_sim_first[, seq(model$d+1, 2*model$d)] # keep x_{t+1} values
    vine_sim_extra <- NULL
    if(k.m > 1){
      vine_sim_extra <- lapply(
        seq(2, k.m),
        function(i){
          # generate x_{t+i}
          runif_extra_dim <- matrix(runif(n*model$d), ncol=model$d)
          vine_sim_extra <- rvinecopulib::inverse_rosenblatt(
            u = cbind(vine_sim_base, runif_extra_dim),
            model = model$timewise[[i]],
            cores = n_cores
          )
          return(vine_sim_extra[,seq(model$d+1, 2*model$d)]) # throw away first model$d components
        }
      )
    }

    vine_sim_vals <- abind::abind(
      vine_sim_base, vine_sim_first, along=3
    ) # concatenates to (n, d , 2)
    vine_sim_vals <- aperm(vine_sim_vals, c(3, 2, 1)) # switch to (2, d, n)

    if(!is.null(vine_sim_extra)){
      vine_sim_extra <- abind::abind(vine_sim_extra, along=3)
      vine_sim_extra <- aperm(vine_sim_extra, c(3, 2, 1))
      assertthat::are_equal(dim(vine_sim_extra), c(k.m-2, model$d, n))
      vine_sim_vals <- abind::abind(vine_sim_vals, vine_sim_extra, along=1) # (k, d, n)
    }

    assertthat::are_equal(dim(vine_sim_vals), c(k.m, model$d, n))

    return(vine_sim_vals)
  }

  if(inherits(model, "cond_timewise_vine")){
    n_half <- round(n/2)
    vine_sim_vals_above <- model_simulation(n = n_half, model$above_timewise) # (k, d, n_half)
    vine_sim_vals_below <- model_simulation(n = n-n_half, model$below_timewise) # (k, d, n-n_half)
    vine_sim_vals <- abind::abind(vine_sim_vals_above, vine_sim_vals_below, along=3)
    return(vine_sim_vals)
  }

  stop('Model type not implemented!')
}
