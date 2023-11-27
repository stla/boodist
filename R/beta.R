#' @title Non-central beta distribution
#' @description A R6 class to represent a non-central beta distribution.
#' @export
#' @importFrom R6 R6Class
#' @importFrom stats rbeta
Beta <- R6Class(
  "Beta",

  private = list(
    ".a"    = NA_real_,
    ".b"    = NA_real_,
    ".delta" = NA_real_
  ),

  active = list(

    #' @field a Get or set the value of \code{a}.
    "a" = function(value) {
      if(missing(value)) {
        return(private[[".a"]])
      } else {
        stopifnot(value > 0)
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
    },

    #' @field delta Get or set the value of \code{delta}.
    "delta" = function(value) {
      if(missing(value)) {
        return(private[[".delta"]])
      } else {
        private[[".delta"]] <- value
      }
    }
  ),

  public = list(

    #' @description New beta distribution.
    #' @param a,b shape parameters, \code{> 0}
    #' @param delta non-centrality parameter, \code{>= 0}
    #' @return A \code{Beta} object.
    "initialize" = function(a, b, delta) {
      stopifnot(a > 0, b > 0, delta >= 0)
      private[[".a"]]     <- a
      private[[".b"]]     <- b
      private[[".delta"]] <- delta
    },

    #' @description Density function of the beta distribution.
    #' @param x numeric vector
    #' @return The density evaluated at \code{x}.
    "d" = function(x) {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      rcpp_dbeta(x, a, b, delta)
    },

    #' @description Cumulative distribution function of the beta
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      rcpp_beta(q, a, b, delta, lower)
    },

    #' @description Quantile function of the beta distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      rcpp_qbeta(p, a, b, delta, lower)
    },

    #' @description Sampling from the beta distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      rbeta(n, a, b, ncp = delta)
    },

    #' @description Mean of the beta distribution.
    #' @return The mean of the beta distribution.
    "mean" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      beta_mean(a, b, delta)
    },

    #' @description Median of the beta distribution.
    #' @return The median of the beta distribution.
    "median" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      beta_median(a, b, delta)
    },

    #' @description Mode of the beta distribution.
    #' @return The mode of the beta distribution.
    "mode" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      beta_mode(a, b, delta)
    },

    #' @description Standard deviation of the beta distribution.
    #' @return The standard deviation of the beta distribution.
    "sd" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      beta_sd(a, b, delta)
    },

    #' @description Variance of the beta distribution.
    #' @return The variance of the beta distribution.
    "variance" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      beta_variance(a, b, delta)
    },

    #' @description Skewness of the beta distribution.
    #' @return The skewness of the beta distribution.
    "skewness" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      if(delta > 0) {
        stop("Skewness of the non-central beta distribution is not available.")
      }
      beta_skewness(a, b)
    },

    #' @description Kurtosis of the beta distribution.
    #' @return The kurtosis of the beta distribution.
    "kurtosis" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      if(delta > 0) {
        stop("Kurtosis of the non-central beta distribution is not available.")
      }
      beta_kurtosis(a, b)
    },

    #' @description Kurtosis excess of the beta distribution.
    #' @return The kurtosis excess of the beta distribution.
    "kurtosisExcess" = function() {
      a     <- private[[".a"]]
      b     <- private[[".b"]]
      delta <- private[[".delta"]]
      if(delta > 0) {
        stop(
          "Kurtosis excess of the non-central beta distribution",
          " is not available."
        )
      }
      beta_kurtosis_excess(a, b)
    }
  )
)
