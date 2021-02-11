proba_necessary_causation <- function(p_causal, p_counterfactual){
  return(
    pmax(1.-p_factual/p_counterfactual, 0.)
  )
}

function(w)
