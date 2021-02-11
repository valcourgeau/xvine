proba_single_conditional_single_target <- function(data, col_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))

  above_target_data <- data[2:k, col_target, above_idx]

  print(dim(above_target_data))
  if(is.vector(above_target_data)){
    p_above <- mean(above_target_data > u0_target)
  }else{
    p_above <- apply(above_target_data, 1, function(x){mean(x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k, col_target, -above_idx]
  if(is.vector(below_target_data)){
    p_below <- mean(below_target_data > u0_target)
  }else{
    p_below <- apply(below_target_data, 1, function(x){mean(x > u0_target)}) # filter target
  }


  return(cbind(p_above, p_below))
}

proba_single_conditional_cross_target <- function(data, w_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  # w_target are not checked to sum to 1
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)
  stopifnot(length(w_target) == d)
  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))

  above_target_data <- data[2:k,, above_idx] # remove time t and keeps all marginals
  if(is.vector(above_target_data)){
    p_above <- mean(above_target_data > u0_target)
  }else{
    p_above <- apply(above_target_data, 1, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k,, -above_idx]
  if(is.vector(below_target_data)){
    p_below <- mean(below_target_data > u0_target)
  }else{
    p_below <- apply(below_target_data, 1, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  return(cbind(p_above, p_below))
}

proba_single_conditional_time_target <- function(data, w_target, col_source, u0_target, u0_source){
  # data: matrix k * d * n
  # w_target are not checked to sum to 1
  k <- dim(data)[1]
  d <- dim(data)[2]
  stopifnot(k >= 2)
  stopifnot(length(w_target) == k - 1)

  above_idx <- data[1,col_source,] > u0_source
  stopifnot(any(above_idx))

  above_target_data <- data[2:k,, above_idx] # remove time t and keeps all marginals
  if(is.vector(above_target_data)){
    p_above <- mean(above_target_data > u0_target)
  }else{
    p_above <- apply(above_target_data, 2, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  below_target_data <- data[2:k,, -above_idx]
  if(is.vector(below_target_data)){
    p_below <- mean(below_target_data > u0_target)
  }else{
    p_below <- apply(below_target_data, 2, function(x){rowMeans(w_target %*% x > u0_target)}) # filter target
  }

  return(cbind(p_above, p_below))
}
