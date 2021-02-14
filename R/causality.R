proba_necessary_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-p_factual/(p_counterfactual+.Machine$double.eps), 0.)
  )
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
      return(prod(pn[times]))
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
      return(prod(pn[col_target]))
    }
  )
}



