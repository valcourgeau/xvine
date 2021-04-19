
forward_stop <- function(p_values, level=0.05){
  # implements the ForwardStop rule from
  # G’SELL, M. G., WAGER, S., CHOULDECHOVA, A. and TIBSHIRANI, R. (2016).
  # Sequential selection procedures and false discovery rate control.
  # J. R. Stat. Soc. Ser. B. Stat. Methodol. 78
  # 423–444. MR3454203
  k <- length(p_values)
  stopping_vals <- - cumsum(log(pmax(.Machine$double.eps, 1-p_values))) / k

  return(which.max(stopping_vals < level))
}


