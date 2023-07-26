#' @importFrom stats approxfun
#' @noRd
qbounds_distr2 <- function(distr, p) {
  mn  <- distr$mean()
  std <- distr$sd()
  q_ <- seq(max(0, mn - 6*std), mn + 6*std, length.out = 200L)
  prob_ <- distr$p(q_)
  f <- approxfun(prob_, q_)
  guess <- f(p)
  absguess <- abs(guess)
  c(max(0, guess - 0.1*absguess), guess + 0.1*absguess)
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
#'   enclosing the values of the quantiles or choose the automatic bounds.
#'   A maximum number of iterations is fixed in the root-finding algorithm.
#'   If it is reached, a warning is thrown.
#' @export
#' @importFrom R6 R6Class
#' @examples
#' if(require("plotly")) {
#' library(boodist)
#'
#' x_ <- seq(0, 3, length.out = 100L)
#' lambda_ <- seq(-1, 1, length.out = 100L)
#' dsty <- vapply(lambda_, function(lambda) {
#'   GeneralizedInverseGaussian$new(theta = 1, eta = 1, lambda)$d(x_)
#' }, numeric(length(x_)))
#' #
#' txt <- matrix(NA_character_, nrow = length(x_), ncol = length(lambda_))
#' for(i in 1L:nrow(txt)) {
#'   for(j in 1L:ncol(txt)) {
#'     txt[i, j] <- paste0(
#'       "x: ", formatC(x_[i]),
#'       "<br> lambda: ", formatC(lambda_[j]),
#'       "<br> density: ", formatC(dsty[i, j])
#'     )
#'   }
#' }
#' #
#' plot_ly(
#'   x = ~lambda_, y = ~x_, z = ~dsty, type = "surface",
#'   text = txt, hoverinfo = "text", showscale = FALSE
#' ) %>% layout(
#'   title = "Generalized inverse Gaussian distribution",
#'   margin = list(t = 40, r= 5, b = 5, l = 5),
#'   scene = list(
#'     xaxis = list(
#'       title = "lambda"
#'     ),
#'     yaxis = list(
#'       title = "x"
#'     ),
#'     zaxis = list(
#'       title = "density"
#'     )
#'   )
#' )
#' }
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
    #' @param x vector of positive numbers
    #' @param log Boolean, whether to return the log-density
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x, log = FALSE) {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      if(log) {
        -log(2 * besselK(theta, lambda)) + (lambda - 1)*log(x/eta) -
          theta*(x/eta + eta/x)/2 - log(eta)
      } else {
        ifelse(
          x == 0,
          0,
          1 / (2 * besselK(theta, lambda)) * (x/eta)^(lambda-1) *
            exp(-theta*(x/eta + eta/x)/2) / eta
        )
      }
    },

    #' @description Cumulative distribution function of the generalized inverse
    #'   Gaussian distribution.
    #' @param q numeric vector of quantiles (\code{>= 0})
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
      eta * besselK(theta, lambda + 1) / besselK(theta, lambda)
    },

    #' @description Mode of the generalized inverse Gaussian distribution.
    #' @return The mode of the generalized inverse Gaussian distribution.
    "mode" = function() {
      theta  <- private[[".theta"]]
      eta    <- private[[".eta"]]
      lambda <- private[[".lambda"]]
      eta * (lambda - 1 + sqrt((lambda-1)^2 + theta^2)) / theta
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
      K <- besselK(theta, lambda)
      eta^2 *
        (besselK(theta, lambda + 2) / K - (besselK(theta, lambda + 1) / K)^2)
    }

  )
)
