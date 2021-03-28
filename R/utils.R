logsumexp <- function(x) {
  y <- max(x)
  return(y + log(sum(exp(x - y))))
}

softmax <- function(x) {
  return(exp(x - logsumexp(x)))
}

logits <- function(w) {
  w <- softmax(w)
  return(log(w))
}

cutoff <- function(w, cutoff.value=NULL) {
  if(!is.null(cutoff.value)){
    w[w < cutoff.value] <- 0
    w <- w / sum(w)
  }

  return(w)
}

vnorm <- function(x, p){
  return(sum(abs(x)^p)^{1./p})
}
