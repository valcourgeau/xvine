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
  thres <- quantile(x, u0)
  p_ecdf <- partial_ecdf(x, u0)
  p_gpd <- partial_gpd(x, u0)
  x[x <= thres] <- p_ecdf[['data']] * u0
  x[x > thres] <- p_gpd[['data']] * (1. - u0) + u0
  return(
    list(
      'data'=x,
      'par.ests'=p_gpd$par.ests,
      'par.ses'=p_gpd$par.ses,
      'ecdf'=p_ecdf$ecdf,
      'u0'=u0
    )
  )
}

# Applying integral transform on each marginal
apply_integral_transform <- function(data, u0s){
  stopifnot(ncol(data) == length(u0s))
  # apply it on each marginal
  data_it <- lapply(
    1:length(u0s),
    function(i){integral_transform(data[,i], u0 = u0s[i])}
  )

  data_it_values <- do.call(cbind, lapply(data_it, function(x){x$data}))
  data_it_ests <- do.call(rbind, lapply(data_it, function(x){x$par.ests}))
  data_it_ses <- do.call(rbind, lapply(data_it, function(x){x$par.ses}))
  data_it_ecdf <- do.call(cbind, lapply(data_it, function(x){x$ecdf}))
  data_it_u0 <- do.call(cbind, lapply(data_it, function(x){x$u0}))

  return(
    list(
      'data'=data_it_values,
      'par.ests'=data_it_ests,
      'par.ses'=data_it_ses,
      'ecdf'=data_it_ecdf,
      'u0s'=data_it_u0
    )
  )
}

#TODO write tests
reverse_integral_transform <- function(x, x_source, u0, shape, scale){
  stopifnot(0 <= u0 & u0 <= 1)
  reversed_data <- rep(0, length(x))
  thres <- quantile(x_source, u0)

  reversed_data[x <= u0] <- quantile(
    x_source[x_source <= u0],
    probs=x[x <= u0]
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

#TODO write tests
apply_reverse_integral_transform <- function(data_unif, data_source, u0s, shapes, scales){
  stopifnot(ncol(data_unif) == ncol(data_source))
  stopifnot(ncol(data_unif) == length(u0s))
  stopifnot(length(shapes) == ncol(data_unif))
  stopifnot(length(scales) == ncol(data_unif))
  stopifnot(0 <= min(data_unif) & max(data_unif) <= 1)

  # apply it on each marginal
  data_rit <- lapply(
    1:length(u0s),
    function(i){
      reverse_integral_transform(
        x = data_unif[,i],
        x_source = data_source[,i],
        u0 = u0s[i],
        shape = shapes[i],
        scale = scales[i]
      )
    }
  )

  data_rit_values <- do.call(cbind, lapply(data_rit, function(x) x$data))
  return(
    list(
      'data'=data_rit_values,
      'data_source'=data_source,
      'u0s'=u0s
    )
  )
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

build_conditional_stack <- function(data, k, col, u0){
  # k-window of data where the colomn col is above quantile(u0)
  # data is a matrix
  # k is the length of the window
  # col idx
  # u0 quantile

  thres <- quantile(data[,col], u0)
  idx <- which(data[,col] >= thres)
  dt_stack <- build_stack(data, k)

  return(dt_stack[idx,])
}
