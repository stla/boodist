#' @title Student distribution
#' @description A R6 class to represent a Student distribution.
#' @export
#' @importFrom R6 R6Class
Student <- R6Class(
  "Student",

  private = list(
    ".nu"    = NA_real_,
    ".delta" = NA_real_
  ),

  active = list(

    #' @field nu Get or set the value of \code{nu}.
    "nu" = function(value) {
      if(missing(value)) {
        return(private[[".nu"]])
      } else {
        private[[".nu"]] <- value
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

    #' @description New Student distribution.
    #' @param nu degrees of freedom parameter, \code{>0}
    #' @param delta non-centrality parameter
    #' @return A \code{Student} object.
    "initialize" = function(nu, delta) {
      stopifnot(nu > 0)
      private[[".nu"]]    <- nu
      private[[".delta"]] <- delta
    },

    #' @description Density function of the Student distribution.
    #' @param x numeric vector
    #' @return The density evaluated at \code{x}.
    "d" = function(x) {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      rcpp_dt(x, nu, delta)
    },

    #' @description Cumulative distribution function of the Student
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      rcpp_pt(q, nu, delta, lower)
    },

    #' @description Quantile function of the Student distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      rcpp_qt(p, nu, delta, lower)
    },

    #' @description Sampling from the Student distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      rt(n, df = nu, ncp = delta)
    },

    #' @description Mean of the Student distribution.
    #' @return The mean of the Student distribution.
    "mean" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_mean(nu, delta)
    },

    #' @description Median of the Student distribution.
    #' @return The median of the Student distribution.
    "median" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_median(nu, delta)
    },

    #' @description Mode of the Student distribution.
    #' @return The mode of the Student distribution.
    "mode" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_mode(nu, delta)
    },

    #' @description Standard deviation of the Student distribution.
    #' @return The standard deviation of the Student distribution.
    "sd" = function() {
      sqrt(self$sd())
    },

    #' @description Variance of the Student distribution.
    #' @return The variance of the Student distribution.
    "variance" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_variance(nu, delta)
    },

    #' @description Skewness of the Student distribution.
    #' @return The skewness of the Student distribution.
    "skewness" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_skewness(nu, delta)
    },

    #' @description Kurtosis of the Student distribution.
    #' @return The kurtosis of the Student distribution.
    "kurtosis" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_kurtosis(nu, delta)
    },

    #' @description Kurtosis excess of the Student distribution.
    #' @return The kurtosis excess of the Student distribution.
    "kurtosisExcess" = function() {
      nu    <- private[[".nu"]]
      delta <- private[[".delta"]]
      t_kurtosis_excess(nu, delta)
    }
  )
)
