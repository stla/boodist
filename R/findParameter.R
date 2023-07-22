#' @title Find non-centrality parameter
#' @description Find the non-centrality parameter of a Chi-squared distribution
#'   given a quantile, its corresponding probability, and the degrees of freedom.
#'
#' @param df degrees of freedom, a positive number
#' @param q a quantile
#' @param p probability corresponding to the quantile \code{q}
#'
#' @return The non-centrality parameter of the Chi-squared distribution with
#'   degrees of freedom parameter \code{df} and with cumulative probability
#'   \code{p} at the quantile \code{q}.
#' @export
#'
#' @examples
#' ncp <- findChi2ncp(1, 3, 0.1)
#' pchisq(3, df = 1, ncp = ncp) # should be 0.1
findChi2ncp <- function(df, q, p) {
  stopifnot(df > 0)
  stopifnot(q >= 0)
  stopifnot(p >= 0, p <= 1)
  find_chisq_ncp(nu, q, p)
}

#' @title Find degrees of freedom
#' @description Find the degrees of freedom parameter of a non-central
#'   Chi-squared distribution given a quantile, its corresponding
#'   probability, and the non-centrality parameter.
#'
#' @param ncp non-centrality parameter, a non-negative number
#' @param q a quantile
#' @param p probability corresponding to the quantile \code{q}
#'
#' @return The degrees of freedom parameter of the non-central Chi-squared
#'   distribution with non-centrality parameter \code{ncp} and with
#'   cumulative probability \code{p} at the quantile \code{q}.
#' @export
#'
#' @examples
#' nu <- findChi2df(ncp = 10, q = 3, p = 0.1)
#' pchisq(3, df = nu, ncp = 10) # should be 0.1
findChi2df <- function(ncp, q, p) {
  stopifnot(ncp >= 0)
  stopifnot(q >= 0)
  stopifnot(p >= 0, p <= 1)
  find_chisq_df(ncp, q, p)
}
