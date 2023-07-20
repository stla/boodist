#' @importFrom stats rnorm runif
#' @noRd
rig <- function(n, mu, lambda) {
  stopifnot(mu > 0)
  stopifnot(lambda > 0)
  nu <- rnorm(n)
  y  <- nu * nu
  x <- 1 + (mu*y - sqrt(4*mu*lambda*y + mu*mu*y*y)) / (2*lamba)
  z <- runif(n)
  ifelse(z < 1/(1+x), x, 1/x) * mu
}

#' @importFrom stats rnorm
#' @noRd
rnig <- function(n, alpha, beta, delta, mu) {
  stopifnot(alpha > beta)
  stopifnot(delta > 0)
  z <- rig(n, delta, sqrt(alpha*alpha - beta*beta))
  rnorm(n, mu + beta*z, sqrt(z))
}
