#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dig(
    Rcpp::NumericVector x, double mu, double lambda, bool log
){
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::inverse_gaussian dist(mu, lambda);
  if(log) {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::logpdf(dist, x(i));
    }
  } else {
    for(int i = 0; i < n; i++) {
      out(i) = boost::math::pdf(dist, x(i));
    }
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_pig(
  Rcpp::NumericVector q, double mu, double lambda, bool lower
){
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::inverse_gaussian dist(mu, lambda);
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
Rcpp::NumericVector rcpp_qig(
  Rcpp::NumericVector p, double mu, double lambda, bool lower
){
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::inverse_gaussian dist(mu, lambda);
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
double ig_mean(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return mean(dist);
}

// [[Rcpp::export]]
double ig_median(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return median(dist);
}

// [[Rcpp::export]]
double ig_mode(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return mode(dist);
}

// [[Rcpp::export]]
double ig_variance(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return variance(dist);
}

// [[Rcpp::export]]
double ig_skewness(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return skewness(dist);
}

// [[Rcpp::export]]
double ig_kurtosis(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double ig_kurtosis_excess(double mu, double lambda) {
  boost::math::inverse_gaussian dist(mu, lambda);
  return kurtosis_excess(dist);
}
