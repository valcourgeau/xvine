proba_necessary_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-p_counterfactual/(p_factual+.Machine$double.eps), 0.)
  )
}

proba_sufficient_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-(1.-p_counterfactual)/(1-p_factual+.Machine$double.eps), 0.)
  )
}

proba_necessary_sufficient_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(p_factual-p_counterfactual, 0.)
  )
}

maximise_poc <- function(data, target, col_source, u0_target, u0_source, poc=NULL, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  # target can be times or col_target
  # type can be 'time', 'cross' & 'all'
  stopifnot(type %in% c('time', 'cross', 'all'))

  k.m <- dim(data)[1]
  d <- dim(data)[2]
  if(is.null(poc)){
    poc <- proba_necessary_causation
  }

  if(type == 'time'){
    wrap_pn <- wrapper_pn_tt(data, target, col_source, u0_target, u0_source, poc=poc)
    w_t <- rep(1, k.m - 1)
  }else{
    if(type == 'cross'){
      wrap_pn <- wrapper_pn_ct(data, target, col_source, u0_target, u0_source, poc=poc)
      w_t <- rep(1, d)
    }else{
      if(type == 'all'){
        wrap_pn <- wrapper_pn_all(data, col_source, u0_target, u0_source, poc=poc, lambda=lambda, p=p, cutoff.value=cutoff.value)
        w_t <- rep(1, d * (k.m - 1))
      }
    }
  }

  if(routine == 'optim'){
    optim_routine <- optim(
      par = w_t, fn = function(w){-wrap_pn(w)}
    ) # maximising PN
    if(optim_routine$convergence != 0) warning(paste('maximise_pn: optim routine has not converged for target', target))
    return(list('weights'=softmax(optim_routine$par)))
  }
  if(routine == 'deoptim'){
    optim_routine <- DEoptim::DEoptim(
      lower = rep(-4, length(w_t)),
      upper = rep(0, length(w_t)),
      fn = function(w){-wrap_pn(w)},
      control = DEoptim::DEoptim.control(trace = FALSE)
    ) # maximising PN
    return(list('weights'=softmax(optim_routine$optim$bestmem)))
  }

  stop('Not Implemented')
}

maximise_pn <- function(data, target, col_source, u0_target, u0_source, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  return(
    maximise_poc(
      data, target, col_source, u0_target, u0_source,
      poc=proba_necessary_causation, type=type, routine=routine,
      lambda=lambda, p=p, cutoff.value=cutoff.value
    )
  )
}

maximise_ps <- function(data, target, col_source, u0_target, u0_source, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  return(
    maximise_poc(
      data, target, col_source, u0_target, u0_source,
      poc=proba_sufficient_causation, type=type, routine=routine,
      lambda=lambda, p=p, cutoff.value=cutoff.value)
  )
}

maximise_pns <- function(data, col_source, u0_target, u0_source, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  return(
    maximise_poc(
      data=data, col_source=col_source, u0_target=u0_target, u0_source=u0_source,
      poc=proba_necessary_sufficient_causation, type=type, routine=routine,
      lambda=lambda, p=p, cutoff.value=cutoff.value)
  )
}

assemble_pn <- function(pn){
  return(sum(log(pn + .Machine$double.eps)))
}

wrapper_pn_ct <- function(data, times, col_source, u0_target, u0_source, poc=NULL, lambda=NULL, p=1, cutoff.value=NULL){
  # wraps the single conditional cross target
  # weights are standardised using softmax
  # if times is a vector, returns the PN product

  if(is.null(poc)){
    poc <- proba_necessary_causation
  }

  f <- wrapper_ct(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- softmax(w)
      w <- cutoff(w, cutoff.value=cutoff.value)
      causal_p <- f(w)
      pn <- poc(
        causal_p$factual,
        causal_p$counterfactual
      )

      if(is.null(lambda)){
        return(assemble_pn(pn[times]))
      }else{
        return(assemble_pn(pn[times]) - lambda * vnorm(w, p))
      }
    }
  )
}

wrapper_pn_tt <- function(data, col_target, col_source, u0_target, u0_source, poc=NULL, lambda=NULL, p=1, cutoff.value=NULL){
  # wraps the single conditional time target
  # weights are standardised using
  # if col_target is a vector, returns the PN product
  if(is.null(poc)){poc <- proba_necessary_causation}

  f <- wrapper_tt(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- softmax(w)
      w <- cutoff(w, cutoff.value=cutoff.value)
      causal_p <- f(w)
      pn <- poc(
        causal_p$factual,
        causal_p$counterfactual
      )
      if(is.null(lambda)){
        return(assemble_pn(pn[col_target]))
      }else{
        return(assemble_pn(pn[col_target]) - lambda * vnorm(w, p))
      }
    }
  )
}

wrapper_pn_all <- function(data, col_source, u0_target, u0_source, poc=NULL, lambda=NULL, p=1, cutoff.value=NULL){
  # wraps the single conditional time target
  # weights are standardised using a softmax
  if(is.null(poc)){poc <- proba_necessary_causation}

  f <- wrapper_all(data, col_source, u0_target, u0_source) # wrapper
  return(
    function(w){
      w <- as.vector(w)
      w <- softmax(w)
      w <- cutoff(w, cutoff.value=cutoff.value)
      causal_p <- f(w)
      pn <- poc(
        causal_p$factual,
        causal_p$counterfactual
      )
      if(is.null(lambda)){
        return(assemble_pn(pn))
      }else{
        return(assemble_pn(pn) - lambda * vnorm(w, p))
      }
    }
  )
}
