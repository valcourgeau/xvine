# Implement metrics between distribution from samples


binning <- function(x, ...){
  # x is a vector
  stopifnot(is.vector(x))
  return(hist(x, plot=FALSE, ...)$counts)
}

coercive_binning <- function(x, y, per_bin=10){
  # returns two rows of binned data
  b_from  <- min(min(x), min(y))
  b_to <- max(max(x), max(y))
  n <- min(length(x), length(y))
  binned_data <- rbind(
    binning(x, breaks=seq(from=b_from, to=b_to, length.out=round(n/per_bin)+1)),
    binning(y, breaks=seq(from=b_from, to=b_to, length.out=round(n/per_bin)+1))
  )
  return(binned_data)
}

distr_distance <- function(x, y, method, ...){
  as.numeric(
    philentropy::distance(
      coercive_binning(x, y, ...), method=method, est.prob='empirical'
    )
  )
}

apply_distr_distance <- function(x, y, method, ...){
  # x y are arrays, if 3d, they are (k, d, n)
  # distr_distance are computed on rows
  stopifnot(is.array(x), is.array(y))
  dim_x <- dim(x)
  k <- dim_x[1]
  d <- dim_x[2]

  if(length(dim_x) == 3){
    # we apply apply_distr_distance to second dim
    arr_dist <- vapply(
      seq_len(d), FUN.VALUE = rep(0, k),
      FUN = function(j){
        apply_distr_distance(
          x = x[,j,],
          y = y[,j,],
          method = method,
          ...
        )
      }
    )
    return(arr_dist)
  }

  if(nrow(x) == nrow(y)){
    res <- vapply(
      seq_len(nrow(x)),
      function(i){distr_distance(x[i,], y[i,], method = method, ...)},
      1.0
    )
    return(res)
  }

  if(ncol(x) == ncol(y)){
    return(apply_distr_distance(t(x), t(y), method, ...))
  }

  stop('Not Implemented')
}



