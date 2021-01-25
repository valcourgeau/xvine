uniform_ecdf <- function(x){
  if(is.matrix(x)){
    return(apply(
        x, MARGIN = 2,
        FUN = function(x_col){ecdf(x_col)(x_col)}
      )
    )
  }
}
