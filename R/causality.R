proba_necessary_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-p_counterfactual/(p_factual+.Machine$double.eps), 0.)
  )
}

proba_sufficient_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(1.-(1.-p_factual)/(1-p_counterfactual+.Machine$double.eps), 0.)
  )
}

proba_necessary_sufficient_causation <- function(p_factual, p_counterfactual){
  return(
    pmax(p_factual-p_counterfactual, 0.)
  )
}

maximise_poc <- function(data, col_source, u0_target, u0_source, target=NULL, poc=NULL, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  # target can be times or col_target
  # type can be 'time', 'cross' & 'all'
  stopifnot(type %in% c('time', 'cross', 'all', 'all-but-source'))
  stopifnot(routine %in% c('optim', 'deoptim', 'optimParallel'))

  k.m <- dim(data)[1]
  d <- dim(data)[2]
  if(is.null(poc)){
    poc <- proba_necessary_causation
  }

  if(type == 'time'){
    wrap_pn <- wrapper_pn_tt(data, target, col_source, u0_target, u0_source, poc=poc)
    w_t <- rep(.5, k.m - 1)
  }else{
    if(type == 'cross'){
      wrap_pn <- wrapper_pn_ct(data, target, col_source, u0_target, u0_source, poc=poc)
      w_t <- rep(.5, d)
    }else{
      if(type == 'all'){
        wrap_pn <- wrapper_pn_all(data, col_source, u0_target, u0_source, poc=poc, lambda=lambda, p=p, cutoff.value=cutoff.value)
        w_t <- rep(.5, d * (k.m - 1))
      }else{
        if(type == 'all-but-source'){
          wrap_pn <- wrapper_pn_all_but_source(data, col_source, u0_target, u0_source, poc=poc, lambda=lambda, p=p, cutoff.value=cutoff.value)
          w_t <- rep(.5, (d-1) * (k.m - 1))
        }
      }
    }
  }

  if(routine == 'optim'){
    optim_routine <- optim(
      par = w_t,
      fn = function(w){-wrap_pn(w)},
      method = 'L-BFGS-B',
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      control = list(trace=1)
    ) # maximising PC
    return(list('weights'=optim_routine$par))
  }

  if(routine == 'optimParallel'){
    poc_env <- new.env()
    poc_env$wrap_fn_e <- wrap_pn
    fn_optim <- function(w){-wrap_fn_e(w)}
    clu <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(clu, c('wrap_fn_e'), poc_env)
    tryCatch({
      optim_routine <- optimParallel::optimParallel(
        par = w_t,
        fn = fn_optim,
        method = 'L-BFGS-B',
        lower = rep(0, length(w_t)),
        upper = rep(1, length(w_t)),
        control = list(trace=1),
        parallel = list(cl=clu, forward=F, loginfo=F)
      ) # maximising PC
    },
    finally={
      parallel::stopCluster(clu)
    })

    return(list('weights'=optim_routine$par))
  }
  if(routine == 'deoptim'){
    optim_routine <- DEoptim::DEoptim(
      lower = rep(0, length(w_t)),
      upper = rep(1, length(w_t)),
      fn = function(w){-wrap_pn(w)},
      control = DEoptim::DEoptim.control(trace = FALSE, parallelType = 1)
    ) # maximising PC
    return(list('weights'=optim_routine$optim$bestmem)) # softmax(optim_routine$optim$bestmem)))
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
      data=data, target=target,col_source=col_source, u0_target, u0_source,
      poc=proba_sufficient_causation, type=type, routine=routine,
      lambda=lambda, p=p, cutoff.value=cutoff.value)
  )
}

maximise_pns <- function(data, target, col_source, u0_target, u0_source, type = 'time', routine='optim', lambda=NULL, p=1, cutoff.value=NULL){
  return(
    maximise_poc(
      data=data, target=target, col_source=col_source, u0_target=u0_target, u0_source=u0_source,
      poc=proba_necessary_sufficient_causation, type=type, routine=routine,
      lambda=lambda, p=p, cutoff.value=cutoff.value)
  )
}

assemble_pn <- function(pn, log.p = TRUE){
  poc_v <- log(pn + .Machine$double.eps)
  if(log.p){
    return(sum(poc_v))
  } else {
    return(sum(exp(poc_v)))
  }
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
      if(sum(abs(w)) < .Machine$double.eps){
        w <- w + .Machine$double.eps
      }
      causal_p <- f(w / sum(w)) # TODO used to be f(w)
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
      w <- as.vector(w)
      if(sum(abs(w)) < .Machine$double.eps){
        w <- w + .Machine$double.eps
      }
      causal_p <- f(w / sum(w)) # TODO used to be f(w)
      poc_value <- poc(
        causal_p$factual,
        causal_p$counterfactual
      )
      if(is.null(lambda)){
        return(assemble_pn(poc_value[col_target]))
      }else{
        return(assemble_pn(poc_value[col_target]) - lambda * vnorm(w, p))
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
      if(sum(abs(w)) < .Machine$double.eps){
        w <- w + .Machine$double.eps
        return(-1000)
      }
      causal_p <- f(w / sum(w)) # TODO used to be f(w)
      poc_value <- poc(causal_p$factual, causal_p$counterfactual)
      print(paste('wrapper_pn_all', 'poc_value', poc_value))
      if(is.null(lambda)){
        return(assemble_pn(poc_value))
      }else{
        return(assemble_pn(poc_value) - lambda * vnorm(w, p))
      }
    }
  )
}

wrapper_pn_all_but_source <- function(data, col_source, u0_target, u0_source, poc=NULL, lambda=NULL, p=1, cutoff.value=NULL){
  # wraps the single conditional time target
  # weights are standardised using a softmax
  if(is.null(poc)){poc <- proba_necessary_causation}

  k <- dim(data)[1]
  d <- dim(data)[2]

  # takes wrapper_pn_all as backend
  wrap_f <- wrapper_pn_all(data, col_source, u0_target, u0_source, poc=poc, lambda=lambda, p=p, cutoff.value=cutoff.value)
  return(
    function(w){
      w <- as.vector(w)

      w2 <- rep(0, d*(k-1))
      idx_to_add <- (1+(col_source-1)*(k-1)):(col_source*(k-1))
      w2[idx_to_add] <- rep(0, k-1)
      w2[-idx_to_add] <- w

      return(wrap_f(w2))
    }
  )
}
