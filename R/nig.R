#' @importFrom stats rnorm
#' @noRd
rnig <- function(n, mu, alpha, beta, delta) {
  stopifnot(alpha > 0, alpha > beta)
  stopifnot(delta > 0)
  z <- rig(n, 1/sqrt(alpha*alpha - beta*beta), delta)
  v <- delta * z
  rnorm(n, mu + beta * v, sqrt(v))
}

#' @title Normal-inverse Gaussian distribution
#' @description A R6 class to represent a normal-inverse Gaussian distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Normal-inverse_Gaussian_distribution}{Wikipedia}.
#' @note
#' The cumulative distribution function is evaluated by integrating the
#'   density function (in C++). Its returned value has two attributes: a
#'   numeric vector \code{"error_estimate"} and an integer vector
#'   \code{"error_code"}. The error code is 0 if no problem is detected. If an
#'   error code is not 0, a warning is thrown. The quantile function is
#'   evaluated by root-finding and then the user must provide some bounds
#'   enclosing the values of the quantiles. A maximum number of iterations is
#'   fixed in the root-finding algorithm. If it is reached, a warning is
#'   thrown.
#' @export
#' @importFrom R6 R6Class
NormalInverseGaussian <- R6Class(
  "NormalInverseGaussian",

  private = list(
    ".mu"    = NA_real_,
    ".alpha" = NA_real_,
    ".beta"  = NA_real_,
    ".delta" = NA_real_,
    ".gamma" = NA_real_
  ),

  active = list(

    #' @field mu Get or set the value of \code{mu}.
    "mu" = function(value) {
      if(missing(value)) {
        return(private[[".mu"]])
      } else {
        private[[".mu"]] <- value
      }
    },

    #' @field alpha Get or set the value of \code{alpha}.
    "alpha" = function(value) {
      if(missing(value)) {
        return(private[[".alpha"]])
      } else {
        private[[".alpha"]] <- value
      }
    },

    #' @field beta Get or set the value of \code{beta}.
    "beta" = function(value) {
      if(missing(value)) {
        return(private[[".beta"]])
      } else {
        private[[".beta"]] <- value
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

    #' @description New normal-inverse Gaussian distribution.
    #' @param mu location parameter
    #' @param alpha tail heaviness parameter, \code{>0}
    #' @param beta asymmetry parameter
    #' @param delta scale parameter, \code{>0}
    #' @return An \code{inverseGaussian} object.
    "initialize" = function(mu, alpha, beta, delta) {
      stopifnot(alpha > 0, alpha > beta)
      stopifnot(delta > 0)
      private[[".mu"]]    <- mu
      private[[".alpha"]] <- alpha
      private[[".beta"]]  <- beta
      private[[".delta"]] <- delta
      private[[".gamma"]] <- sqrt(alpha*alpha - beta*beta)
    },

    #' @description Density function of the normal-inverse
    #'   Gaussian distribution.
    #' @param x numeric vector
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x) {
      mu    <- private[[".mu"]]
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      gamma <- private[[".gamma"]]
      s <- sqrt(delta^2 + (x-mu)^2)
      alpha*delta*besselK(alpha*s, 1) * exp(delta*gamma+beta*(x-mu)) / (pi*s)
    },

    #' @description Cumulative distribution function of the normal-inverse
    #'   Gaussian distribution.
    #' @param q numeric vector of quantiles
    #' @return The cumulative probabilities corresponding to \code{q}, with two
    #'   attributes (see the \strong{Note}).
    "p" = function(q) {
      mu    <- private[[".mu"]]
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      pnig_rcpp(q, mu, alpha, beta, delta)
    },

    #' @description Quantile function of the normal-inverse
    #'   Gaussian distribution.
    #' @param p numeric vector of probabilities
    #' @param a,b bounds enclosing the quantiles to be found (see the
    #'   \strong{Note})
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, a = -30, b = 30) {
      if(any(self$p(a) - p >= 0)) {
        stop("The lower bound `a` is too large.")
      }
      if(any(self$p(b) - p <= 0)) {
        stop("The upper bound `b` is too small.")
      }
      mu    <- private[[".mu"]]
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      qnig_rcpp(p, a, b, mu, alpha, beta, delta)
    },

    #' @description Sampling from the normal-inverse Gaussian distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      mu    <- private[[".mu"]]
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      rnig(n, mu, alpha, beta, delta)
    },

    #' @description Mean of the normal-inverse Gaussian distribution.
    #' @return The mean of the normal-inverse Gaussian distribution.
    "mean" = function() {
      mu    <- private[[".mu"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      gamma <- private[[".gamma"]]
      mu + delta * beta / gamma
    },

    #' @description Standard deviation of the normal-inverse Gaussian
    #'  distribution.
    #' @return The standard deviation of the normal-inverse Gaussian
    #'   distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the normal-inverse Gaussian distribution.
    #' @return The variance of the normal-inverse Gaussian distribution.
    "variance" = function() {
      alpha <- private[[".alpha"]]
      delta <- private[[".delta"]]
      gamma <- private[[".gamma"]]
      delta * alpha * alpha / (gamma * gamma * gamma)
    },

    #' @description Skewness of the normal-inverse Gaussian distribution.
    #' @return The skewness of the normal-inverse Gaussian distribution.
    "skewness" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      gamma <- private[[".gamma"]]
      3 * beta / (alpha * sqrt(delta * gamma))
    },

    #' @description Kurtosis of the normal-inverse Gaussian distribution.
    #' @return The kurtosis of the normal-inverse Gaussian distribution.
    "kurtosis" = function() {
      self$kurtosisExcess() + 3
    },

    #' @description Kurtosis excess of the normal-inverse Gaussian distribution.
    #' @return The kurtosis excess of the normal-inverse Gaussian distribution.
    "kurtosisExcess" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      delta <- private[[".delta"]]
      gamma <- private[[".gamma"]]
      3 * (1 + 4 * beta*beta/(alpha*alpha)) / (delta*gamma)
    }
  )
)
