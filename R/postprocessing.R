
split_data <- function(data, col_source, u0_source){
  # splits data of col_source below and above u0_source
  # taking the minimum of the two
  k <- dim(data)[1]
  n <- dim(data)[3]
  stopifnot(k >= 2)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))
  stopifnot(any(!above_idx)) # check some are below

  if(sum(above_idx) >= n/2){
    # takes length of complement
    below_target_data <- data[2:k, col_target, !above_idx]
    update_above_idx <- which(above_idx)
    above_idx <- sample(update_above_idx, size = sum(!above_idx), replace = F)
    above_target_data <- data[2:k, col_target, above_idx]
  }else{
    # takes length of above_idx
    above_target_data <- data[2:k, col_target, above_idx]
    below_idx <- which(!above_idx)
    below_idx <- sample(below_idx, size = sum(above_idx), replace = F)
    below_target_data <- data[2:k, col_target, below_idx]
  }

  return(
    list(
      'above'=above_target_data,
      'below'=below_target_data
    )
  )
}

proba_single_target <- function(data, col_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  k <- dim(data)[1]
  d <- dim(data)[2]
  n <- dim(data)[3]
  stopifnot(k >= 2)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))
  stopifnot(any(!above_idx)) # check some are below

  if(sum(above_idx) >= n/2){
    below_target_data <- data[2:k, col_target, !above_idx]
    update_above_idx <- which(above_idx)
    above_idx <- sample(update_above_idx, size = sum(!above_idx), replace = F)
    above_target_data <- data[2:k, col_target, above_idx]
  }else{
    above_target_data <- data[2:k, col_target, above_idx]
    below_idx <- which(!above_idx)
    below_idx <- sample(below_idx, size = sum(above_idx), replace = F)
    below_target_data <- data[2:k, col_target, below_idx]
  }

  # above_target_data <- data[2:k, col_target, above_idx]
  if(is.vector(above_target_data)){
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- apply(above_target_data, 1, function(x){mean(x > u0_target)}) # filter target
  }

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
    stop('Not Implemented')
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- apply(above_target_data, 1, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k,, !above_idx]
  if(is.vector(below_target_data)){
    stop('Not Implemented')
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

  below_target_data <- data[2:k,, !above_idx]
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

proba_all_target <- function(data, w_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  # w_target are not checked to sum to 1
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)
  if(is.matrix(w_target)){
    stopifnot(dim(w_target)[1] == k - 1)
    stopifnot(dim(w_target)[2] == d)
  }else{
    if(is.vector(w_target)){
      # weights are assumed to be by column
      stopifnot(length(w_target) == d * (k - 1))
      w_target <- matrix(w_target, nrow = k - 1, ncol = d, byrow = F)
    }else{
      stop('w_target must be a vector or matrix of size (k-1, d).')
    }
  }

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx)) # check some are above
  stopifnot(any(!above_idx))# check some are below

  above_target_data <- data[2:k,, above_idx] # remove time t and keeps all marginals
  if(is.vector(above_target_data)){
    stop('Not Implemented')
    p_factual <- mean(above_target_data > u0_target)
  }else{
    p_factual <- mean(apply(above_target_data, 3, function(x){sum(w_target * x) > u0_target})) # filter target
  }

  below_target_data <- data[2:k,, !above_idx]
  if(is.vector(below_target_data)){
    stop('Not Implemented')
    p_counterfactual <- mean(below_target_data > u0_target)
  }else{
    p_counterfactual <- mean(apply(below_target_data, 3, function(x){sum(w_target * x) > u0_target})) # filter target
  }

  return(list('factual'=p_factual, 'counterfactual'=p_counterfactual))
}

wrapper_all <- function(data, col_source, u0_target, u0_source){
  #returns a function only on parameters
  return(
    function(w){
      proba_all_target(
        data=data, w_target=w, col_source=col_source, u0_target=u0_target, u0_source=u0_source
      )
    }
  )
}
