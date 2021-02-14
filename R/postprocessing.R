proba_single_target <- function(data, col_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))
  stopifnot(any(!above_idx))# check some are below

  above_target_data <- data[2:k, col_target, above_idx]

  if(is.vector(above_target_data)){
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- apply(above_target_data, 1, function(x){mean(x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k, col_target, -above_idx]
  if(is.vector(below_target_data)){
    p_counterfactual <- mean(below_target_data > u0_target)
  }else{
    p_counterfactual <- apply(below_target_data, 1, function(x){mean(x > u0_target)}) # filter target
  }

  return(list('factual'=p_factual, 'counterfactual'=p_counterfactual))
}

proba_cross_target <- function(data, w_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  # w_target are not checked to sum to 1
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)
  stopifnot(length(w_target) == d)
  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))
  stopifnot(any(!above_idx))# check some are below

  above_target_data <- data[2:k,, above_idx] # remove time t and keeps all marginals
  if(is.vector(above_target_data)){
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- apply(above_target_data, 1, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k,, -above_idx]
  if(is.vector(below_target_data)){
    p_counterfactual <- mean(below_target_data > u0_target)
  }else{
    p_counterfactual <- apply(below_target_data, 1, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  return(list('factual'=p_factual, 'counterfactual'=p_counterfactual))
}

wrapper_ct <- function(data, col_source, u0_target, u0_source){
  #returns a function only on parameters
  return(
    function(w){
      proba_cross_target(
        data, w, col_source, u0_target, u0_source
      )
    }
  )
}

proba_time_target <- function(data, w_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  # w_target are not checked to sum to 1
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)
  stopifnot(length(w_target) == k - 1)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx)) # check some are above
  stopifnot(any(!above_idx))# check some are below

  above_target_data <- data[2:k,, above_idx] # remove time t and keeps all marginals
  if(is.vector(above_target_data)){
    stop('Not Implemented')
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- apply(above_target_data, 2, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k,, -above_idx]
  if(is.vector(below_target_data)){
    stop('Not Implemented')
    p_counterfactual <- mean(below_target_data > u0_target)
  }else{
    p_counterfactual <- apply(below_target_data, 2, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  return(list('factual'=p_factual, 'counterfactual'=p_counterfactual))
}

wrapper_tt <- function(data, col_source, u0_target, u0_source){
  #returns a function only on parameters
  return(
    function(w){
      proba_time_target(
        data, w, col_source, u0_target, u0_source
      )
    }
  )
}
