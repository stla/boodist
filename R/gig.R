#' @importFrom stats approxfun
#' @noRd
qbounds_distr2 <- function(distr, p) {
  mn  <- distr$mean()
  std <- distr$sd()
  q_ <- seq(mn - 6*std, mn + 6*std, length.out = 200L)
  prob_ <- distr$p(q_)
  f <- approxfun(prob_, q_)
  guess <- f(p)
  absguess <- abs(guess)
  c(guess - 0.1*absguess, guess + 0.1*absguess)
}

#' @title Generalized inverse Gaussian distribution
#' @description A R6 class to represent a generalized inverse Gaussian
#'   distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution}{Wikipedia}.
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
GeneralizedInverseGaussian <- R6Class(
  "GeneralizedInverseGaussian",

  private = list(
    ".theta"  = NA_real_,
    ".eta"    = NA_real_,
    ".lambda" = NA_real_
  ),

  active = list(

    #' @field theta Get or set the value of \code{theta}.
    "theta" = function(value) {
      if(missing(value)) {
        return(private[[".theta"]])
      } else {
        stopifnot(value > 0)
        private[[".theta"]] <- value
      }
    },

    #' @field eta Get or set the value of \code{eta}.
    "eta" = function(value) {
      if(missing(value)) {
        return(private[[".eta"]])
      } else {
        stopifnot(value > 0)
        private[[".eta"]] <- value
      }
    },

    #' @field lambda Get or set the value of \code{lambda}.
    "lambda" = function(value) {
      if(missing(value)) {
        return(private[[".lambda"]])
      } else {
        private[[".lambda"]] <- value
      }
    }
  ),

  public = list(

    #' @description New generalized inverse Gaussian distribution.
    #' @param theta concentration parameter, \code{>0}
    #' @param eta scale parameter, \code{>0}
    #' @param lambda parameter (denoted by \code{p} on Wikipedia)
    #' @return A \code{GeneralizedInverseGaussian} object.
    "initialize" = function(theta, eta, lambda) {
      stopifnot(theta > 0)
      stopifnot(eta > 0)
      private[[".theta"]]  <- theta
      private[[".eta"]]    <- eta
      private[[".lambda"]] <- lambda
    },

    #' @description Density function of the generalized inverse
    #'   Gaussian distribution.
    #' @param x numeric vector
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x) {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      1 / (2 * besselK(theta, lambda)) * (x/eta)^(lambda-1) *
        exp(-(x/eta + eta/x)/2) / eta
    },

    #' @description Cumulative distribution function of the generalized inverse
    #'   Gaussian distribution.
    #' @param q numeric vector of quantiles
    #' @return The cumulative probabilities corresponding to \code{q}, with two
    #'   attributes (see the \strong{Note}).
    "p" = function(q) {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      pgig_rcpp(q, theta, eta, lambda)
    },

    #' @description Quantile function of the generalized inverse
    #'   Gaussian distribution.
    #' @param p numeric vector of probabilities
    #' @param bounds bounds enclosing the quantiles to be found (see the
    #'   \strong{Note}), or \code{NULL} for automatic bounds
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, bounds = NULL) {
      if(!is.null(bounds)) {
        a <- bounds[, 1L]
        b <- bounds[, 2L]
        if(any(self$p(a) - p >= 0)) {
          stop("The lower bound is too large.")
        }
        if(any(self$p(b) - p <= 0)) {
          stop("The upper bound is too small.")
        }
      } else {
        bounds <- vapply(p, function(prob) {
          qbounds_distr2(self, prob)
        }, numeric(2L))
        a <- bounds[1L, ]
        b <- bounds[2L, ]
      }
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      qgig_rcpp(p, a, b, theta, eta, lambda)
    },

    #' @description Sampling from the generalized inverse Gaussian distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      rgig_rcpp(n, lambda, theta) * eta
    },

    #' @description Mean of the generalized inverse Gaussian distribution.
    #' @return The mean of the generalized inverse Gaussian distribution.
    "mean" = function() {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
    },

    #' @description Standard deviation of the generalized inverse Gaussian
    #'  distribution.
    #' @return The standard deviation of the generalized inverse Gaussian
    #'   distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the generalized inverse Gaussian distribution.
    #' @return The variance of the generalized inverse Gaussian distribution.
    "variance" = function() {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
    },

    #' @description Skewness of the generalized inverse Gaussian distribution.
    #' @return The skewness of the generalized inverse Gaussian distribution.
    "skewness" = function() {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
    },

    #' @description Kurtosis of the generalized inverse Gaussian distribution.
    #' @return The kurtosis of the generalized inverse Gaussian distribution.
    "kurtosis" = function() {
      self$kurtosisExcess() + 3
    },

    #' @description Kurtosis excess of the generalized inverse Gaussian distribution.
    #' @return The kurtosis excess of the generalized inverse Gaussian distribution.
    "kurtosisExcess" = function() {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
    }
  )
)
