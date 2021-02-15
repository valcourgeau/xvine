logsumexp <- function(x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function(x) {
  exp(x - logsumexp(x))
}

logits <- function(w) {
  w <- softmax(w)
  return(log(w))
}
