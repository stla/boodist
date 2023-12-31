% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gig.R
\name{GeneralizedInverseGaussian}
\alias{GeneralizedInverseGaussian}
\title{Generalized inverse Gaussian distribution}
\description{
A R6 class to represent a generalized inverse Gaussian
  distribution.
}
\details{
See \href{https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution}{Wikipedia}.
}
\note{
The cumulative distribution function is evaluated by integrating the
  density function (in C++). Its returned value has two attributes: a
  numeric vector \code{"error_estimate"} and an integer vector
  \code{"error_code"}. The error code is 0 if no problem is detected. If an
  error code is not 0, a warning is thrown. The quantile function is
  evaluated by root-finding and then the user must provide some bounds
  enclosing the values of the quantiles or choose the automatic bounds.
  A maximum number of iterations is fixed in the root-finding algorithm.
  If it is reached, a warning is thrown.
}
\examples{
if(require("plotly")) {
library(boodist)

x_ <- seq(0, 3, length.out = 100L)
lambda_ <- seq(-1, 1, length.out = 100L)
dsty <- vapply(lambda_, function(lambda) {
  GeneralizedInverseGaussian$new(theta = 1, eta = 1, lambda)$d(x_)
}, numeric(length(x_)))
#
txt <- matrix(NA_character_, nrow = length(x_), ncol = length(lambda_))
for(i in 1L:nrow(txt)) {
  for(j in 1L:ncol(txt)) {
    txt[i, j] <- paste0(
      "x: ", formatC(x_[i]),
      "<br> lambda: ", formatC(lambda_[j]),
      "<br> density: ", formatC(dsty[i, j])
    )
  }
}
#
plot_ly(
  x = ~lambda_, y = ~x_, z = ~dsty, type = "surface",
  text = txt, hoverinfo = "text", showscale = FALSE
) \%>\% layout(
  title = "Generalized inverse Gaussian distribution",
  margin = list(t = 40, r= 5, b = 5, l = 5),
  scene = list(
    xaxis = list(
      title = "lambda"
    ),
    yaxis = list(
      title = "x"
    ),
    zaxis = list(
      title = "density"
    )
  )
)
}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{theta}}{Get or set the value of \code{theta}.}

\item{\code{eta}}{Get or set the value of \code{eta}.}

\item{\code{lambda}}{Get or set the value of \code{lambda}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-GeneralizedInverseGaussian-new}{\code{GeneralizedInverseGaussian$new()}}
\item \href{#method-GeneralizedInverseGaussian-d}{\code{GeneralizedInverseGaussian$d()}}
\item \href{#method-GeneralizedInverseGaussian-p}{\code{GeneralizedInverseGaussian$p()}}
\item \href{#method-GeneralizedInverseGaussian-q}{\code{GeneralizedInverseGaussian$q()}}
\item \href{#method-GeneralizedInverseGaussian-r}{\code{GeneralizedInverseGaussian$r()}}
\item \href{#method-GeneralizedInverseGaussian-mean}{\code{GeneralizedInverseGaussian$mean()}}
\item \href{#method-GeneralizedInverseGaussian-mode}{\code{GeneralizedInverseGaussian$mode()}}
\item \href{#method-GeneralizedInverseGaussian-sd}{\code{GeneralizedInverseGaussian$sd()}}
\item \href{#method-GeneralizedInverseGaussian-variance}{\code{GeneralizedInverseGaussian$variance()}}
\item \href{#method-GeneralizedInverseGaussian-clone}{\code{GeneralizedInverseGaussian$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-new"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-new}{}}}
\subsection{Method \code{new()}}{
New generalized inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$new(theta, eta, lambda)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{theta}}{concentration parameter, \code{>0}}

\item{\code{eta}}{scale parameter, \code{>0}}

\item{\code{lambda}}{parameter (denoted by \code{p} on Wikipedia)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{GeneralizedInverseGaussian} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-d"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-d}{}}}
\subsection{Method \code{d()}}{
Density function of the generalized inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$d(x, log = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{vector of positive numbers}

\item{\code{log}}{Boolean, whether to return the log-density}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The density or the log-density evaluated at \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-p"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-p}{}}}
\subsection{Method \code{p()}}{
Cumulative distribution function of the generalized inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$p(q)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{q}}{numeric vector of quantiles (\code{>= 0})}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The cumulative probabilities corresponding to \code{q}, with two
  attributes (see the \strong{Note}).
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-q"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-q}{}}}
\subsection{Method \code{q()}}{
Quantile function of the generalized inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$q(p, bounds = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{p}}{numeric vector of probabilities}

\item{\code{bounds}}{bounds enclosing the quantiles to be found (see the
\strong{Note}), or \code{NULL} for automatic bounds}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The quantiles corresponding to \code{p}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-r"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-r}{}}}
\subsection{Method \code{r()}}{
Sampling from the generalized inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$r(n)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{n}}{number of simulations}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A numeric vector of length \code{n}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-mean"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-mean}{}}}
\subsection{Method \code{mean()}}{
Mean of the generalized inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$mean()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The mean of the generalized inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-mode"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-mode}{}}}
\subsection{Method \code{mode()}}{
Mode of the generalized inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$mode()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The mode of the generalized inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-sd"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-sd}{}}}
\subsection{Method \code{sd()}}{
Standard deviation of the generalized inverse Gaussian
 distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$sd()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The standard deviation of the generalized inverse Gaussian
  distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-variance"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-variance}{}}}
\subsection{Method \code{variance()}}{
Variance of the generalized inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$variance()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The variance of the generalized inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-GeneralizedInverseGaussian-clone"></a>}}
\if{latex}{\out{\hypertarget{method-GeneralizedInverseGaussian-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{GeneralizedInverseGaussian$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
