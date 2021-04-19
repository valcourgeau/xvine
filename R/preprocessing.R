# fits and applies ecdf on each column
uniform_ecdf <- function(x){
  ecdf_handle <-  function(x_col){ecdf(x_col)(x_col)}
  if(is.matrix(x)){
    return(apply(
        x, MARGIN = 2,
        FUN = ecdf_handle
      )
    )
  }else{
    if(is.list(x)){
      return(lapply(
          x, FUN = ecdf_handle
        )
      )
    }
  }
}

#apply ecdf on data below quantile u0
partial_ecdf <- function(x, u0){
  # u0 is a quantile value between 0 and 1
  # x is a vector
  y <- x[x <= quantile(x, u0)]
  f <- ecdf(y)
  return(list('data'=f(y), 'ecdf'=f))
}

# GPD MLE for data above quantile u0
partial_gpd <- function(x, u0){
  # u0 is a quantile value between 0 and 1
  # x is a vector
  thres <- quantile(x, u0)
  gpd_fit <- evir::gpd(data = x, threshold = thres)
  return(
    list(
      'data'=evir::pgpd(gpd_fit$data-thres, xi = gpd_fit$par.ests[1], beta = gpd_fit$par.ests[2]),
      'par.ests'=gpd_fit$par.ests,
      'par.ses'=gpd_fit$par.ses
    )
  )
}

# Integral transform
integral_transform <- function(x, u0){
  # u0 is a quantile value between 0 and 1
  # x is a vector
  x_it <- rep(0, length(x))
  thres <- quantile(x, u0)
  p_ecdf <- partial_ecdf(x, u0)
  p_gpd <- partial_gpd(x, u0)

  x_it[x <= thres] <- p_ecdf[['data']] * u0
  x_it[x > thres] <- p_gpd[['data']] * (1. - u0) + u0
  return(
    list(
      'data'=x_it,
      'par.ests'=p_gpd$par.ests,
      'par.ses'=p_gpd$par.ses,
      'u0'=u0
    )
  )
}

# Applying integral transform on each marginal
# TODO add ecdf which can compute quantile
apply_integral_transform <- function(data, u0s){
  stopifnot(ncol(data) == length(u0s))
  has_3_dims <- length(dim(data)) == 3
  if(has_3_dims){
    k <- dim(data)[1]
    d <- dim(data)[2]
    n <- dim(data)[3]
    data_2 <- apply(data, 2, cbind) # (k*n, d)
  }else{
    data_2 <- data
  }

  # apply it on each marginal
  data_it <- lapply(
    1:length(u0s),
    function(i){integral_transform(data_2[,i], u0 = u0s[i])}
  )

  data_it_values <- do.call(cbind, lapply(data_it, function(x){x$data}))
  data_it_ests <- do.call(rbind, lapply(data_it, function(x){x$par.ests}))
  data_it_ses <- do.call(rbind, lapply(data_it, function(x){x$par.ses}))
  data_it_u0 <- do.call(cbind, lapply(data_it, function(x){x$u0}))

  if(has_3_dims){
    data_it_values <- array(data_it_values, c(k, n, d))
    data_it_values <- aperm(data_it_values, c(1, 3, 2))
  }

  return(
    list(
      'data'=data_it_values,
      'par.ests'=data_it_ests,
      'par.ses'=data_it_ses,
      'u0s'=data_it_u0
    )
  )
}

reverse_integral_transform <- function(x, x_source, u0, shape, scale){
  stopifnot(0 <= u0 & u0 <= 1)
  reversed_data <- rep(0, length(x))
  thres <- quantile(x_source, u0)

  reversed_data[x <= u0] <- quantile(
    x_source[x_source <= thres],
    probs=x[x <= u0]/u0
  )

  reversed_data[x > u0] <- evir::qgpd(
    p = (x[x > u0]-u0)/(1-u0), mu = thres, xi = shape, beta = scale
  )
  return(
    list(
      'data'=reversed_data
    )
  )
}

apply_reverse_integral_transform <- function(data_unif, data_source, u0s, shapes, scales){
  stopifnot(ncol(data_unif) == ncol(data_source))
  stopifnot(ncol(data_unif) == length(u0s))
  stopifnot(length(shapes) == ncol(data_unif))
  stopifnot(length(scales) == ncol(data_unif))
  stopifnot(0 <= min(data_unif) & max(data_unif) <= 1)

  has_3_dims <- length(dim(data_unif)) == 3

  if(has_3_dims){
    k <- dim(data_unif)[1]
    d <- dim(data_unif)[2]
    n <- dim(data_unif)[3]
    data_unif_2 <- apply(data_unif, 2, cbind) # concatenate into (k*n, d)
  }else{
    data_unif_2 <- data_unif
  }

  data_source <- apply(data_source, 2, cbind)
  print(dim(data_unif_2))
  # apply it on each marginal
  data_rit <- lapply(
    1:length(u0s),
    function(i){
      reverse_integral_transform(
        x = data_unif_2[,i],
        x_source = data_source[,i],
        u0 = u0s[i],
        shape = shapes[i],
        scale = scales[i]
      )
    }
  )

  data_rit_values <- do.call(cbind, lapply(data_rit, function(x) x$data)) # (k*n, d)

  if(has_3_dims){
    data_rit_values <- array(data_rit_values, c(k, n, d))
    data_rit_values <- aperm(data_rit_values, c(1, 3, 2))
  }

  return(
    list(
      'data'=data_rit_values,
      'data_source'=data_source,
      'u0s'=u0s
    )
  )
}


apply_dual_reverse_integral_transform <- function(data_unif, data_source, u0s, conditional_on, shapes_f, scales_f, shapes_cf, scales_cf){
  stopifnot(ncol(data_unif) == ncol(data_source))
  stopifnot(ncol(data_unif) == length(u0s))
  stopifnot(length(shapes_f) == ncol(data_unif))
  stopifnot(length(shapes_cf) == ncol(data_unif))
  stopifnot(length(scales_f) == ncol(data_unif))
  stopifnot(length(scales_cf) == ncol(data_unif))
  stopifnot(0 <= min(data_unif) & max(data_unif) <= 1)

  has_3_dims <- length(dim(data_unif)) == 3

  if(has_3_dims){
    k <- dim(data_unif)[1]
    d <- dim(data_unif)[2]
    n <- dim(data_unif)[3]
    idx_f <- which(data_unif[1, conditional_on,] > u0s[conditional_on])
    idx_cf <- which(data_unif[1, conditional_on,] <= u0s[conditional_on])

    data_rit_f <- apply_reverse_integral_transform(
      data_unif = data_unif[,, idx_f],
      data_source = data_source,
      u0s = u0s,
      shapes = shapes_f,
      scales = scales_f
    )$data
    print(paste('data_rit_f', dim(data_rit_f)))

    data_rit_cf <- apply_reverse_integral_transform(
      data_unif = data_unif[,,idx_cf],
      data_source = data_source,
      u0s = u0s,
      shapes = shapes_cf,
      scales = scales_cf
    )$data
    print(paste('data_rit_cf', dim(data_rit_cf)))

    tmp <- array(0, dim(data_unif))
    tmp[,,idx_f] <- data_rit_f
    tmp[,,idx_cf] <- data_rit_cf
  }else{
    stop('Not Implemented')
  }


  return(
    list(
      'data'=tmp,
      'data_source'=data_source,
      'u0s'=u0s,
      'scales_f'=scales_f,
      'scales_cf'=scales_cf,
      'shapes_f'=shapes_f,
      'shapes_cf'=shapes_cf
    )
  )
}

# TODO test
reverse_exponential <- function(x){
  # takes uniform(0, 1) and transforms them into Exp(1) RV.
  stopifnot(min(x) >= 0, max(x) <= 1)
  return(list('data'=stats::qexp(x, rate = 1)))
}

# TODO test
apply_reverse_exponential <- function(data_unif){
  stopifnot(0 <= min(data_unif), max(data_unif) <= 1)

  has_3_dims <- length(dim(data_unif)) == 3
  if(has_3_dims){
    k <- dim(data_unif)[1]
    d <- dim(data_unif)[2]
    n <- dim(data_unif)[3]
    data_unif_2 <- apply(data_unif, 2, cbind) # concatenate into (k*n, d)
  }else{
    data_unif_2 <- data_unif
  }

  data_source <- apply(data_unif, 2, cbind)

  # apply it on each marginal
  data_rit <- lapply(
    seq_len(ncol(data_unif)),
    function(i){reverse_exponential(x = data_unif_2[,i])}
  )

  data_rit_values <- do.call(cbind, lapply(data_rit, function(x) x$data)) # (k*n, d)

  if(has_3_dims){
    data_rit_values <- array(data_rit_values, c(k, n, d))
    data_rit_values <- aperm(data_rit_values, c(1, 3, 2))
  }

  return(list('data'=data_rit_values))
}

build_stack <- function(data, k){
  # data is a matrix
  # k is the length of the window

  dt_stack <- zoo::rollapply(
    data, width=k,
    FUN=function(x){c(t(x))},
    by.column=F
  )
  return(dt_stack)
}

build_above_conditional_stack <- function(data, k, col, u0){
  # k-window of data where the colomn col is above quantile(u0)
  # data is a matrix (n, d)
  # k is the length of the window
  # col idx
  # u0 quantile

  thres <- quantile(data[,col], u0)
  idx <- which(data[seq_len(nrow(data)-k), col] > thres)
  dt_stack <- build_stack(data, k)

  return(dt_stack[idx,])
}

build_below_conditional_stack <- function(data, k, col, u0){
  # k-window of data where the colomn col is above quantile(u0)
  # data is a matrix (n, d)
  # k is the length of the window
  # col idx
  # u0 quantile

  thres <- quantile(data[,col], u0)
  idx <- which(data[seq_len(nrow(data)-k), col] <= thres)
  dt_stack <- build_stack(data, k)

  return(dt_stack[idx,])
}

build_timewise_stack <- function(data, k){
  # data is a matrix
  # k is the length of the window
  # returns a list with (x_{t}, x_{t+s}) for 1 <= s <= k
  dt_stack <- lapply(
    seq_len(k-1),
    function(i){
      cbind(data[1:(nrow(data)-k+1),], data[(i+1):(nrow(data)-k+1+i),])
    }
  )

  return(dt_stack)
}

build_above_conditional_timewise_stack <- function(data, k, col, u0){
  # data is a matrix
  # k is the length of the window
  # returns a list with (x_{t}, x_{t+s}) for 1 <= s <= k where x_t above u0 quantile

  thres <- quantile(data[,col], u0)
  idx <- which(data[seq_len(nrow(data)-k), col] > thres)
  dt_stack <- build_timewise_stack(data, k)
  dt_stack <- lapply(dt_stack, function(dt_s) dt_s[idx,])

  return(dt_stack)
}

build_below_conditional_timewise_stack <- function(data, k, col, u0){
  # data is a matrix
  # k is the length of the window
  # returns a list with (x_{t}, x_{t+s}) for 1 <= s <= k where x_t above u0 quantile

  thres <- quantile(data[,col], u0)
  idx <- which(data[seq_len(nrow(data)-k), col] <= thres)
  dt_stack <- build_timewise_stack(data, k)
  dt_stack <- lapply(dt_stack, function(dt_s) dt_s[idx,])

  return(dt_stack)
}
