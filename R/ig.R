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

#' @title Inverse Gaussian distribution
#' @description A R6 class to represent an inverse Gaussian distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution}{Wikipedia}.
#' @export
#' @importFrom R6 R6Class
InverseGaussian <- R6Class(
  "InverseGaussian",

  private = list(
    ".mu"     = NA_real_,
    ".lambda" = NA_real_
  ),

  active = list(

    #' @field mu Get or set the value of \code{mu}.
    "mu" = function(value) {
      if(missing(value)) {
        return(private[[".mu"]])
      } else {
        stopifnot(value > 0)
        private[[".mu"]] <- value
      }
    },

    #' @field lambda Get or set the value of \code{lambda}.
    "lambda" = function(value) {
      if(missing(value)) {
        return(private[[".lambda"]])
      } else {
        stopifnot(value > 0)
        private[[".lambda"]] <- value
      }
    }
  ),

  public = list(

    #' @description New inverse Gaussian distribution.
    #' @param mu parameter, the mean, \code{>0}
    #' @param lambda shape parameter, \code{>0}
    #' @return An \code{inverseGaussian} object.
    "initialize" = function(mu, lambda) {
      stopifnot(lambda > 0)
      private[[".mu"]]     <- mu
      private[[".lambda"]] <- lambda
    },

    #' @description Density function of the inverse Gaussian distribution.
    #' @param x vector of positive numbers
    #' @param log Boolean, whether to return the logarithm of the density
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x, log = FALSE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_dig(x, mu, lambda, log)
    },

    #' @description Cumulative distribution function of the inverse Gaussian
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_pig(q, mu, lambda, lower)
    },

    #' @description Quantile function of the inverse Gaussian distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rcpp_qig(p, mu, lambda, lower)
    },

    #' @description Sampling from the inverse Gaussian distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      rig(n, mu, lambda)
    },

    #' @description Mean of the inverse Gaussian distribution.
    #' @return The mean of the inverse Gaussian distribution.
    "mean" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_mean(mu, lambda)
    },

    #' @description Median of the inverse Gaussian distribution.
    #' @return The median of the inverse Gaussian distribution.
    "median" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_median(mu, lambda)
    },

    #' @description Mode of the inverse Gaussian distribution.
    #' @return The mode of the inverse Gaussian distribution.
    "mode" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_mode(mu, lambda)
    },

    #' @description Standard deviation of the inverse Gaussian distribution.
    #' @return The standard deviation of the inverse Gaussian distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the inverse Gaussian distribution.
    #' @return The variance of the inverse Gaussian distribution.
    "variance" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_variance(mu, lambda)
    },

    #' @description Skewness of the inverse Gaussian distribution.
    #' @return The skewness of the inverse Gaussian distribution.
    "skewness" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_skewness(mu, lambda)
    },

    #' @description Kurtosis of the inverse Gaussian distribution.
    #' @return The kurtosis of the inverse Gaussian distribution.
    "kurtosis" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_kurtosis(mu, lambda)
    },

    #' @description Kurtosis excess of the inverse Gaussian distribution.
    #' @return The kurtosis excess of the inverse Gaussian distribution.
    "kurtosisExcess" = function() {
      mu     <- private[[".mu"]]
      lambda <- private[[".lambda"]]
      ig_kurtosis_excess(mu, lambda)
    }
  )
)
