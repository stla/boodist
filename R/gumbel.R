#' @importFrom stats runif
#' @noRd
rgumbel <- function(n, a, b) {
  stopifnot(b > 0)
  rcpp_qgumbel(runif(n), a, b)
}


#' @title Gumbel distribution
#' @description A R6 class to represent a Gumbel distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Gumbel_distribution}{Wikipedia}.
#' @export
#' @importFrom R6 R6Class
Gumbel <- R6Class(
  "Gumbel",

  private = list(
    ".a" = NA_real_,
    ".b" = NA_real_
  ),

  active = list(

    #' @field a Get or set the value of \code{a}.
    "a" = function(value) {
      if(missing(value)) {
        return(private[[".a"]])
      } else {
        private[[".a"]] <- value
      }
    },

    #' @field b Get or set the value of \code{b}.
    "b" = function(value) {
      if(missing(value)) {
        return(private[[".b"]])
      } else {
        stopifnot(value > 0)
        private[[".b"]] <- value
      }
    }
  ),

  public = list(

    #' @description New Gumbel distribution.
    #' @param a location parameter
    #' @param b scale parameter, \code{>0}
    #' @return A \code{Gumbel} object.
    "initialize" = function(a, b) {
      stopifnot(b > 0)
      private[[".a"]] <- a
      private[[".b"]] <- b
    },

    #' @description Density function of the Gumbel distribution.
    #' @param x numeric vector
    #' @param log Boolean, whether to return the logarithm of the density
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x, log = FALSE) {
      a <- private[[".a"]]
      b <- private[[".b"]]
      rcpp_dgumbel(x, a, b, log)
    },

    #' @description Cumulative distribution function of the Gumbel
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      a <- private[[".a"]]
      b <- private[[".b"]]
      rcpp_pgumbel(q, a, b, lower)
    },

    #' @description Quantile function of the Gumbel distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      a <- private[[".a"]]
      b <- private[[".b"]]
      rcpp_qgumbel(p, a, b, lower)
    },

    #' @description Sampling from the Gumbel distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      a <- private[[".a"]]
      b <- private[[".b"]]
      rgumbel(n, a, b)
    },

    #' @description Mean of the Gumbel distribution.
    #' @return The mean of the Gumbel distribution.
    "mean" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_mean(a, b)
    },

    #' @description Median of the Gumbel distribution.
    #' @return The median of the Gumbel distribution.
    "median" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_median(a, b)
    },

    #' @description Mode of the Gumbel distribution.
    #' @return The mode of the Gumbel distribution.
    "mode" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_mode(a, b)
    },

    #' @description Standard deviation of the Gumbel distribution.
    #' @return The standard deviation of the Gumbel distribution.
    "sd" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_sd(a, b)
    },

    #' @description Variance of the Gumbel distribution.
    #' @return The variance of the Gumbel distribution.
    "variance" = function() {
      self$sd() * self$sd()
    },

    #' @description Skewness of the Gumbel distribution.
    #' @return The skewness of the Gumbel distribution.
    "skewness" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_skewness(a, b)
    },

    #' @description Kurtosis of the Gumbel distribution.
    #' @return The kurtosis of the Gumbel distribution.
    "kurtosis" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_kurtosis(a, b)
    },

    #' @description Kurtosis excess of the Gumbel distribution.
    #' @return The kurtosis excess of the Gumbel distribution.
    "kurtosisExcess" = function() {
      a <- private[[".a"]]
      b <- private[[".b"]]
      gumbel_kurtosis_excess(a, b)
    }
  )
)
