#' @importFrom stats rexp
#' @noRd
rhexp <- function(n, probs, rates) {
  x <- vapply(rates, function(lambda) {
    rexp(n, lambda)
  }, numeric(n))
  i <- sample.int(length(probs), size = n, replace = TRUE, prob = probs)
  x[cbind(seq_len(n), i)]
}

#' @title Hyperexponential distribution
#' @description A R6 class to represent a hyperexponential distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Hyperexponential_distribution}{Wikipedia}.
#'
#' @export
#' @importFrom R6 R6Class
Hyperexponential <- R6Class(
  "Hyperexponential",

  private = list(
    ".probs"    = NA_real_,
    ".rates" = NA_real_
  ),

  active = list(

    #' @field probs Get or set the value of \code{probs}.
    "probs" = function(value) {
      if(missing(value)) {
        return(private[[".probs"]])
      } else {
        private[[".probs"]] <- value
      }
    },

    #' @field rates Get or set the value of \code{rates}.
    "rates" = function(value) {
      if(missing(value)) {
        return(private[[".rates"]])
      } else {
        private[[".rates"]] <- value
      }
    }
  ),

  public = list(

    #' @description New hyperexponential distribution.
    #' @param probs probabilities (weights), a vector of positive numbers
    #' @param rates rate parameters, vector of positive numbers of the same
    #'   length as the \code{probs} vector
    #' @return A \code{Hyperexponential} object.
    "initialize" = function(probs, rates) {
      stopifnot(all(probs > 0))
      stopifnot(all(rates > 0))
      private[[".probs"]] <- probs
      private[[".rates"]] <- rates
    },

    #' @description Density function of the hyperexponential distribution.
    #' @param x vector of positive numbers
    #' @return The density evaluated at \code{x}.
    "d" = function(x) {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      rcpp_dhexp(x, probs, rates)
    },

    #' @description Cumulative distribution function of the hyperexponential
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      rcpp_phexp(q, probs, rates, lower)
    },

    #' @description Quantile function of the hyperexponential distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      rcpp_qhexp(p, probs, rates, lower)
    },

    #' @description Sampling from the hyperexponential distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      rhexp(n, probs/sum(probs), rates)
    },

    #' @description Mean of the hyperexponential distribution.
    #' @return The mean of the hyperexponential distribution.
    "mean" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_mean(probs, rates)
    },

    #' @description Mode of the hyperexponential distribution.
    #' @return The mode of the hyperexponential distribution.
    "mode" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_mode(probs, rates)
    },

    #' @description Standard deviation of the hyperexponential distribution.
    #' @return The standard deviation of the hyperexponential distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the hyperexponential distribution.
    #' @return The variance of the hyperexponential distribution.
    "variance" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_variance(probs, rates)
    },

    #' @description Skewness of the hyperexponential distribution.
    #' @return The skewness of the hyperexponential distribution.
    "skewness" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_skewness(probs, rates)
    },

    #' @description Kurtosis of the hyperexponential distribution.
    #' @return The kurtosis of the hyperexponential distribution.
    "kurtosis" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_kurtosis(probs, rates)
    },

    #' @description Kurtosis excess of the hyperexponential distribution.
    #' @return The kurtosis excess of the hyperexponential distribution.
    "kurtosisExcess" = function() {
      probs <- private[[".probs"]]
      rates <- private[[".rates"]]
      t_kurtosis_excess(probs, rates)
    }
  )
)

