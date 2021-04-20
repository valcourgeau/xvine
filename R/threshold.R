
forward_stop <- function(p_values, level=0.05){
  # implements the ForwardStop rule from
  # G’SELL, M. G., WAGER, S., CHOULDECHOVA, A. and TIBSHIRANI, R. (2016).
  # Sequential selection procedures and false discovery rate control.
  # J. R. Stat. Soc. Ser. B. Stat. Methodol. 78
  # 423–444. MR3454203
  k <- length(p_values)
  stopping_vals <- - cumsum(log(pmax(.Machine$double.eps, 1-p_values))) / seq_len(k)
  idx <- which(stopping_vals[stopping_vals < level])
  if(length(idx) > 0){
    idx <- max(idx)
  }else{
    idx <- NA
  }
  return(list(idx=idx, stopping_values=stopping_vals))
}


gpd_goodness <- function(data, u0s, gof_test, level=0.05, bootstrap = FALSE, bootnum = NULL, allowParallel = FALSE, numCores = 1){
  # wraps eva's Formal (Automated) Goodness-of-Fit Testing
  q_data <- quantile(data, u0s)
  gof_res <- lapply(q_data, function(q_val){gof_test(
    data=data[data > q_val],
    bootstrap=bootstrap,
    bootnum=bootnum,
    allowParallel=allowParallel,
    numCores=numCores)
    }
  )

  results <- list(
    quantiles=q_data,
    u0s=u0s,
    statistics=vapply(gof_res, function(x){x$statistic}, 1.0),
    p.values=vapply(gof_res, function(x){x$p.value}, 1.0),
    thetas=lapply(gof_res, function(x){as.vector(x$theta)}),
    effective_bootnums=as.vector(lapply(gof_res, function(x){x$effective_bootnum}))
  )
  fs <- forward_stop(results$p.values, level)
  results$optimal_idx <- fs$idx
  results$optimal_stopping_values <- fs$stopping_values
  results$optimal_u0 <-  if(length(results$optimal_idx) == 1){u0s[results$optimal_idx]}else{1-level}
  results$optimal_threshold <- if(length(results$optimal_idx) == 1){q_data[results$optimal_idx]}else{quantile(data, 1-level)}
  return(results)
}

apply_gpd_goodness <- function(data, u0s, gof_test, level=0.05, bootstrap = FALSE, bootnum = NULL, allowParallel = FALSE, numCores = 1){
  gof_results <- lapply(
    seq_len(ncol(data)),
    function(i){gpd_goodness(data[,i], u0s, gof_test, level, bootstrap, bootnum, allowParallel, numCores)}
  )
  return(list(
    optimal_u0=vapply(gof_results, function(dt){dt$optimal_u0}, 1.0),
    optimal_threshold=vapply(gof_results, function(dt){dt$optimal_threshold}, 1.0),
    optimal_idx=vapply(gof_results, function(dt){dt$optimal_idx}, 1.0),
    p.values=vapply(gof_results, function(dt){dt$p.values}, rep(0, length(u0s))),
    statistics=vapply(gof_results, function(dt){dt$statistics}, rep(0, length(u0s)))
  ))
}
