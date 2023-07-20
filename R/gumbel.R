rgumbel <- function(n, a, b) {
  stopifnot(b > 0)
  rcpp_qgumbel(runif(n), a, b)
}
