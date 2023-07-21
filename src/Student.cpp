#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dt(
    Rcpp::NumericVector x, double nu, double delta
){
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::non_central_t dist(nu, delta);
  for(int i = 0; i < n; i++) {
    out(i) = boost::math::pdf(dist, x(i));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pt(
  Rcpp::NumericVector q, double nu, double delta, bool lower
){
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::non_central_t dist(nu, delta);
  if(lower) {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::cdf(dist, q(i));
    }
  } else {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::cdf(boost::math::complement(dist, q(i)));
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_qt(
  Rcpp::NumericVector p, double nu, double delta, bool lower
){
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::non_central_t dist(nu, delta);
  if(lower) {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::quantile(dist, p(i));
    }
  } else {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::quantile(boost::math::complement(dist, p(i)));
    }
  }
  return out;
}

// [[Rcpp::export]]
double t_mean(double nu, double delta) {
  boost::math::non_central_t dist(nu, delta);
  return mean(dist);
}

// [[Rcpp::export]]
double t_variance(double nu, double delta) {
  boost::math::non_central_t dist(nu, delta);
  return variance(dist);
}

// [[Rcpp::export]]
double t_skewness(double nu, double delta) {
  boost::math::non_central_t dist(nu, delta);
  return skewness(dist);
}

// [[Rcpp::export]]
double t_kurtosis(double nu, double delta) {
  boost::math::non_central_t dist(nu, delta);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double t_kurtosis_excess(double nu, double delta) {
  boost::math::non_central_t dist(nu, delta);
  return kurtosis_excess(dist);
}
