#' @importFrom stats rexp
#' @noRd
rhexp <- function(n, probs, rates) {
  x <- vapply(rates, function(lambda) {
    rexp(n, lambda)
  }, numeric(n))
  i <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
  x[cbind(seq_len(n), i)]
}

