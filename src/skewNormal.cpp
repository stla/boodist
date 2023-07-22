#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dskewNormal(
    Rcpp::NumericVector x, double xi, double omega, double alpha
){
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::skew_normal dist(xi, omega, alpha);
  for(int i = 0; i < n; i++) {
    out(i) = boost::math::pdf(dist, x(i));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pskewNormal(
  Rcpp::NumericVector q, double xi, double omega, double alpha, bool lower
){
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::skew_normal dist(xi, omega, alpha);
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
Rcpp::NumericVector rcpp_qskewNormal(
  Rcpp::NumericVector p, double xi, double omega, double alpha, bool lower
){
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::skew_normal dist(xi, omega, alpha);
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
double skewNormal_mean(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return mean(dist);
}

// [[Rcpp::export]]
double skewNormal_mode(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return mode(dist);
}

// [[Rcpp::export]]
double skewNormal_variance(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return variance(dist);
}

// [[Rcpp::export]]
double skewNormal_skewness(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return skewness(dist);
}

// [[Rcpp::export]]
double skewNormal_kurtosis(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double skewNormal_kurtosis_excess(double xi, double omega, double alpha) {
  boost::math::skew_normal dist(xi, omega, alpha);
  return kurtosis_excess(dist);
}
