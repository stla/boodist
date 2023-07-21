#' @importFrom stats rnorm runif
#' @noRd
rig <- function(n, mu, lambda) {
  stopifnot(mu > 0)
  stopifnot(lambda > 0)
  nu <- rnorm(n)
  y  <- nu * nu
  x <- 1 + (mu*y - sqrt(4*mu*lambda*y + mu*mu*y*y)) / (2*lambda)
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

#' @title Inverse Gaussian distribution
#' @description A R6 class to represent an inverse Gaussian distribution.
#' @export
#' @importFrom R6 R6Class
inverseGaussian <- R6Class(
  "inverseGaussian",

  private = list(
    ".mu"     = NA_real_,
    ".lambda" = NA_real_
  ),

  active = list(
    "mu" = function(value) {
      if(missing(value)) {
        return(private[[".mu"]])
      } else {
        private[[".mu"]] <- value
      }
    },
    "lambda" = function(value) {
      if(missing(value)) {
        return(private[[".lambda"]])
      } else {
        private[[".lambda"]] <- value
      }
    }
  ),

  public = list(

    "initialize" = function(mu, lambda) {
      stopifnot(lambda > 0)
      private[[".mu"]]     <- mu
      private[[".lambda"]] <- lambda
    },

    "d" = function(x, log = FALSE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_dig(x, mu, lambda, log)
    },

    "p" = function(q, lower = TRUE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_pig(q, mu, lambda, lower)
    },

    "q" = function(p, lower = TRUE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_qig(p, mu, lambda, lower)
    },

    "r" = function(n) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rig(n, mu, lambda)
    },

    "mean" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_mean(mu, lambda)
    },

    "median" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_median(mu, lambda)
    },

    "mode" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_mode(mu, lambda)
    },

    "variance" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_variance(mu, lambda)
    },

    "skewness" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_skewness(mu, lambda)
    },

    "kurtosis" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_kurtosis(mu, lambda)
    },

    "kurtosisExcess" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_kurtosis_excess(mu, lambda)
    }
  )
)
