#' @importFrom stats rnorm runif pnorm
rskewNormal <- function(n, xi, omega, alpha) {
  z <- rnorm(n)
  u <- runif(n)
  epsilon <- 2L*as.integer(u < pnorm(-alpha*z)) - 1L
  xi - omega * epsilon * z
}

#' @title Skew normal distribution
#' @description A R6 class to represent a skew normal distribution.
#' @export
#' @details
#' See \href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{Wikipedia}.
#' @importFrom R6 R6Class
SkewNormal <- R6Class(
  "SkewNormal",

  private = list(
    ".xi"    = NA_real_,
    ".omega" = NA_real_,
    ".alpha" = NA_real_
  ),

  active = list(

    #' @field xi Get or set the value of \code{xi}.
    "xi" = function(value) {
      if(missing(value)) {
        return(private[[".xi"]])
      } else {
        private[[".xi"]] <- value
      }
    },

    #' @field omega Get or set the value of \code{omega}.
    "omega" = function(value) {
      if(missing(value)) {
        return(private[[".omega"]])
      } else {
        private[[".omega"]] <- value
      }
    },

    #' @field alpha Get or set the value of \code{alpha}.
    "alpha" = function(value) {
      if(missing(value)) {
        return(private[[".alpha"]])
      } else {
        private[[".alpha"]] <- value
      }
    }
  ),

  public = list(

    #' @description New skew normal distribution.
    #' @param xi location parameter
    #' @param omega scale parameter, \code{>0}
    #' @param alpha shape parameter
    #' @return A \code{SkewNormal} object.
    "initialize" = function(xi, omega, alpha) {
      stopifnot(omega > 0)
      private[[".xi"]]    <- xi
      private[[".omega"]] <- omega
      private[[".alpha"]] <- alpha
    },

    #' @description Density function of the skew normal distribution.
    #' @param x numeric vector
    #' @return The density evaluated at \code{x}.
    "d" = function(x) {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      rcpp_dskewNormal(x, xi, omega, alpha)
    },

    #' @description Cumulative distribution function of the skew normal
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      rcpp_pskewNormal(q, xi, omega, alpha, lower)
    },

    #' @description Quantile function of the skew normal distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      rcpp_qskewNormal(p, xi, omega, alpha, lower)
    },

    #' @description Sampling from the skew normal distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      rskewNormal(n, xi, omega, alpha)
    },

    #' @description Mean of the skew normal distribution.
    #' @return The mean of the skew normal distribution.
    "mean" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_mean(xi, omega, alpha)
    },

    #' @description Mode of the skew normal distribution.
    #' @return The mode of the skew normal distribution.
    "mode" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_mode(xi, omega, alpha)
    },

    #' @description Standard deviation of the skew normal distribution.
    #' @return The standard deviation of the skew normal distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the skew normal distribution.
    #' @return The variance of the skew normal distribution.
    "variance" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_variance(xi, omega, alpha)
    },

    #' @description Skewness of the skew normal distribution.
    #' @return The skewness of the skew normal distribution.
    "skewness" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_skewness(xi, omega, alpha)
    },

    #' @description Kurtosis of the skew normal distribution.
    #' @return The kurtosis of the skew normal distribution.
    "kurtosis" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_kurtosis(xi, omega, alpha)
    },

    #' @description Kurtosis excess of the skew normal distribution.
    #' @return The kurtosis excess of the skew normal distribution.
    "kurtosisExcess" = function() {
      xi    <- private[[".xi"]]
      omega <- private[[".omega"]]
      alpha <- private[[".alpha"]]
      skewNormal_kurtosis_excess(xi, omega, alpha)
    }
  )
)
