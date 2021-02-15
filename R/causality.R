proba_necessary_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-p_counterfactual/(p_factual+.Machine$double.eps), 0.)
  )
}

maximise_pn <- function(data, target, col_source, u0_target, u0_source, type = 'time'){
  # target can be times or col_target
  # type can be 'time', 'cross' & 'all'
  stopifnot(type %in% c('time', 'cross', 'all'))

  k.m <- dim(data)[1]
  d <- dim(data)[2]

  if(type == 'time'){
    wrap_pn <- wrapper_pn_tt(data, target, col_source, u0_target, u0_source)
    w_t <- rep(1, k.m - 1)
  }else{
    if(type == 'cross'){
      wrap_pn <- wrapper_pn_ct(data, target, col_source, u0_target, u0_source)
      w_t <- rep(1, d)
    }else{
      if(type == 'all'){
        wrap_pn <- wrapper_pn_all(data, col_source, u0_target, u0_source)
        w_t <- rep(1, d * (k.m - 1))
      }
    }
  }

  optim_routine <- optim(
    par = w_t, fn = function(w){-wrap_pn(w)}
  ) # maximising PN
  if(optim_routine$convergence != 0) warning(paste('maximise_pn: optim routine has not converged for target', target))
  return(list('weights'=optim_routine$par))
}

assemble_pn <- function(pn){
  return(sum(log(pn + .Machine$double.eps)))
}

wrapper_pn_ct <- function(data, times, col_source, u0_target, u0_source){
  # wraps the single conditional cross target
  # weights are standardised using softmax
  # if times is a vector, returns the PN product

  f <- wrapper_ct(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- softmax(w)
      causal_p <- f(w)
      pn <- proba_necessary_causation(
        causal_p$factual,
        causal_p$counterfactual
      )
      return(assemble_pn(pn[times]))
    }
  )
}

wrapper_pn_tt <- function(data, col_target, col_source, u0_target, u0_source){
  # wraps the single conditional time target
  # weights are standardised using
  # if col_target is a vector, returns the PN product

  f <- wrapper_tt(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- softmax(w)
      causal_p <- f(w)
      pn <- proba_necessary_causation(
        causal_p$factual,
        causal_p$counterfactual
      )
      return(assemble_pn(pn[col_target]))
    }
  )
}

wrapper_pn_all <- function(data, col_source, u0_target, u0_source){
  # wraps the single conditional time target
  # weights are standardised using a softmax

  f <- wrapper_all(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- as.vector(w)
      w <- softmax(w)
      causal_p <- f(w)
      pn <- proba_necessary_causation(
        causal_p$factual,
        causal_p$counterfactual
      )
      return(assemble_pn(pn))
    }
  )
}



