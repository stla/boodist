#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dbeta(
    Rcpp::NumericVector x, double a, double b, double delta
){
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::beta_distribution dist(a, b);
  for(int i = 0; i < n; i++) {
    out(i) = boost::math::pdf(dist, x(i));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pbeta(
  Rcpp::NumericVector q, double a, double b, bool lower
){
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::beta_distribution dist(a, b);
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
Rcpp::NumericVector rcpp_qbeta(
  Rcpp::NumericVector p, double a, double b, bool lower
){
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::beta_distribution dist(a, b);
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
double beta_mean(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return mean(dist);
}

// [[Rcpp::export]]
double beta_variance(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return variance(dist);
}

// [[Rcpp::export]]
double beta_sd(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return standard_deviation(dist);
}

// [[Rcpp::export]]
double beta_mode(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return mode(dist);
}

// [[Rcpp::export]]
double beta_skewness(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return skewness(dist);
}

// [[Rcpp::export]]
double beta_kurtosis(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double beta_kurtosis_excess(double a, double b) {
  boost::math::beta_distribution dist(a, b);
  return kurtosis_excess(dist);
}
