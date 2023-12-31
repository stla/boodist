% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nig.R
\name{NormalInverseGaussian}
\alias{NormalInverseGaussian}
\title{Normal-inverse Gaussian distribution}
\description{
A R6 class to represent a normal-inverse Gaussian distribution.
}
\details{
See \href{https://en.wikipedia.org/wiki/Normal-inverse_Gaussian_distribution}{Wikipedia}.
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
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{mu}}{Get or set the value of \code{mu}.}

\item{\code{alpha}}{Get or set the value of \code{alpha}.}

\item{\code{beta}}{Get or set the value of \code{beta}.}

\item{\code{delta}}{Get or set the value of \code{delta}.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-NormalInverseGaussian-new}{\code{NormalInverseGaussian$new()}}
\item \href{#method-NormalInverseGaussian-d}{\code{NormalInverseGaussian$d()}}
\item \href{#method-NormalInverseGaussian-p}{\code{NormalInverseGaussian$p()}}
\item \href{#method-NormalInverseGaussian-q}{\code{NormalInverseGaussian$q()}}
\item \href{#method-NormalInverseGaussian-r}{\code{NormalInverseGaussian$r()}}
\item \href{#method-NormalInverseGaussian-mean}{\code{NormalInverseGaussian$mean()}}
\item \href{#method-NormalInverseGaussian-sd}{\code{NormalInverseGaussian$sd()}}
\item \href{#method-NormalInverseGaussian-variance}{\code{NormalInverseGaussian$variance()}}
\item \href{#method-NormalInverseGaussian-skewness}{\code{NormalInverseGaussian$skewness()}}
\item \href{#method-NormalInverseGaussian-kurtosis}{\code{NormalInverseGaussian$kurtosis()}}
\item \href{#method-NormalInverseGaussian-kurtosisExcess}{\code{NormalInverseGaussian$kurtosisExcess()}}
\item \href{#method-NormalInverseGaussian-clone}{\code{NormalInverseGaussian$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-new"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-new}{}}}
\subsection{Method \code{new()}}{
New normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$new(mu, alpha, beta, delta)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mu}}{location parameter}

\item{\code{alpha}}{tail heaviness parameter, \code{>0}}

\item{\code{beta}}{asymmetry parameter}

\item{\code{delta}}{scale parameter, \code{>0}}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A \code{NormalInverseGaussian} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-d"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-d}{}}}
\subsection{Method \code{d()}}{
Density function of the normal-inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$d(x, log = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{numeric vector}

\item{\code{log}}{Boolean, whether to return the logarithm of the density}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The density or the log-density evaluated at \code{x}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-p"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-p}{}}}
\subsection{Method \code{p()}}{
Cumulative distribution function of the normal-inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$p(q)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{q}}{numeric vector of quantiles}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The cumulative probabilities corresponding to \code{q}, with two
  attributes (see the \strong{Note}).
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-q"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-q}{}}}
\subsection{Method \code{q()}}{
Quantile function of the normal-inverse
  Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$q(p, bounds = NULL)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-NormalInverseGaussian-r"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-r}{}}}
\subsection{Method \code{r()}}{
Sampling from the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$r(n)}\if{html}{\out{</div>}}
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
\if{html}{\out{<a id="method-NormalInverseGaussian-mean"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-mean}{}}}
\subsection{Method \code{mean()}}{
Mean of the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$mean()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The mean of the normal-inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-sd"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-sd}{}}}
\subsection{Method \code{sd()}}{
Standard deviation of the normal-inverse Gaussian
 distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$sd()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The standard deviation of the normal-inverse Gaussian
  distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-variance"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-variance}{}}}
\subsection{Method \code{variance()}}{
Variance of the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$variance()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The variance of the normal-inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-skewness"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-skewness}{}}}
\subsection{Method \code{skewness()}}{
Skewness of the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$skewness()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The skewness of the normal-inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-kurtosis"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-kurtosis}{}}}
\subsection{Method \code{kurtosis()}}{
Kurtosis of the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$kurtosis()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The kurtosis of the normal-inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-kurtosisExcess"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-kurtosisExcess}{}}}
\subsection{Method \code{kurtosisExcess()}}{
Kurtosis excess of the normal-inverse Gaussian distribution.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$kurtosisExcess()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The kurtosis excess of the normal-inverse Gaussian distribution.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-NormalInverseGaussian-clone"></a>}}
\if{latex}{\out{\hypertarget{method-NormalInverseGaussian-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{NormalInverseGaussian$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
