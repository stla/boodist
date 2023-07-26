#' @title Inverse Gamma distribution
#' @description A R6 class to represent an inverse Gamma distribution.
#' @details
#' See \href{https://en.wikipedia.org/wiki/Inverse-gamma_distribution}{Wikipedia}.
#' @export
#' @importFrom R6 R6Class
#' @importFrom stats pgamma qgamma rgamma
#' @examples
#' if(require("plotly")) {
#' x_ <- seq(0, 2, length.out = 100L)
#' alpha_ <- seq(0.5, 2.5, length.out = 100L)
#' dsty <- vapply(alpha_, function(alpha) {
#'   InverseGamma$new(alpha, beta = 1)$d(x_)
#' }, numeric(length(x_)))
#' #
#' txt <- matrix(NA_character_, nrow = length(x_), ncol = length(alpha_))
#' for(i in 1L:nrow(txt)) {
#'   for(j in 1L:ncol(txt)) {
#'     txt[i, j] <- paste0(
#'       "x: ", formatC(x_[i]),
#'       "<br> alpha: ", formatC(alpha_[j]),
#'       "<br> density: ", formatC(dsty[i, j])
#'     )
#'   }
#' }
#' #
#' plot_ly(
#'   x = ~alpha_, y = ~x_, z = ~dsty, type = "surface",
#'   text = txt, hoverinfo = "text", showscale = FALSE
#' ) %>% layout(
#'   title = "Inverse Gamma distribution",
#'   margin = list(t = 40, r= 5, b = 5, l = 5),
#'   scene = list(
#'     xaxis = list(
#'       title = "alpha"
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
InverseGamma <- R6Class(
  "InverseGamma",

  private = list(
    ".alpha" = NA_real_,
    ".beta"  = NA_real_
  ),

  active = list(

    #' @field alpha Get or set the value of \code{alpha}.
    "alpha" = function(value) {
      if(missing(value)) {
        return(private[[".alpha"]])
      } else {
        stopifnot(value > 0)
        private[[".alpha"]] <- value
      }
    },

    #' @field beta Get or set the value of \code{beta}.
    "beta" = function(value) {
      if(missing(value)) {
        return(private[[".beta"]])
      } else {
        stopifnot(value > 0)
        private[[".beta"]] <- value
      }
    }
  ),

  public = list(

    #' @description New inverse Gamma distribution.
    #' @param alpha shape parameter, \code{>0}
    #' @param beta scale parameter, \code{>0}
    #' @return An \code{inverseGamma} object.
    "initialize" = function(alpha, beta) {
      stopifnot(alpha > 0)
      stopifnot(beta > 0)
      private[[".alpha"]] <- alpha
      private[[".beta"]]  <- beta
    },

    #' @description Density function of the inverse Gamma distribution.
    #' @param x vector of positive numbers
    #' @param log Boolean, whether to return the logarithm of the density
    #' @return The density or the log-density evaluated at \code{x}.
    "d" = function(x, log = FALSE) {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      if(log) {
        alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(x) - beta/x
      } else {
        ifelse(
          x == 0, 0, beta^alpha / gamma(alpha) * x^(-alpha-1) * exp(-beta / x)
        )
      }
    },

    #' @description Cumulative distribution function of the inverse Gamma
    #'   distribution.
    #' @param q numeric vector of quantiles
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The cumulative probabilities corresponding to \code{q}.
    "p" = function(q, lower = TRUE) {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      pgamma(1/q, shape = alpha, rate = beta, lower.tail = !lower)
    },

    #' @description Quantile function of the inverse Gamma distribution.
    #' @param p numeric vector of probabilities
    #' @param lower Boolean, whether to deal with the lower tail
    #' @return The quantiles corresponding to \code{p}.
    "q" = function(p, lower = TRUE) {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      1 / qgamma(p, shape = alpha, rate = beta, lower.tail = !lower)
    },

    #' @description Sampling from the inverse Gamma distribution.
    #' @param n number of simulations
    #' @return A numeric vector of length \code{n}.
    "r" = function(n) {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      1 / rgamma(n, shape = alpha, rate = beta)
    },

    #' @description Mean of the inverse Gamma distribution.
    #' @return The mean of the inverse Gamma distribution.
    "mean" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      if(alpha <= 1) Inf else {beta / (alpha - 1)}
    },

    #' @description Median of the inverse Gamma distribution.
    #' @return The median of the inverse Gamma distribution.
    "median" = function() {
      self$q(0.5)
    },

    #' @description Mode of the inverse Gamma distribution.
    #' @return The mode of the inverse Gamma distribution.
    "mode" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      beta / (alpha + 1)
    },

    #' @description Standard deviation of the inverse Gamma distribution.
    #' @return The standard deviation of the inverse Gamma distribution.
    "sd" = function() {
      sqrt(self$variance())
    },

    #' @description Variance of the inverse Gamma distribution.
    #' @return The variance of the inverse Gamma distribution.
    "variance" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      if(alpha <= 2) {
        Inf
      } else {
        beta*beta / ((alpha-1)*(alpha-1)*(alpha-2))
      }
    },

    #' @description Skewness of the inverse Gamma distribution.
    #' @return The skewness of the inverse Gamma distribution.
    "skewness" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      if(alpha <= 3) {
        Inf
      } else {
        4 * sqrt(alpha - 2) / (alpha - 3)
      }
    },

    #' @description Kurtosis of the inverse Gamma distribution.
    #' @return The kurtosis of the inverse Gamma distribution.
    "kurtosis" = function() {
      self$kurtosisExcess() + 3
    },

    #' @description Kurtosis excess of the inverse Gamma distribution.
    #' @return The kurtosis excess of the inverse Gamma distribution.
    "kurtosisExcess" = function() {
      alpha <- private[[".alpha"]]
      beta  <- private[[".beta"]]
      if(alpha <= 4) {
        Inf
      } else {
        6 * (5*alpha - 11) / ((alpha - 3)*(alpha - 4))
      }

    }
  )
)
