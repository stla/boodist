#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dhexp(
    Rcpp::NumericVector x, Rcpp::NumericVector probs, Rcpp::NumericVector rates
){
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(probs, rts);
  for(int i = 0; i < n; i++) {
    out(i) = boost::math::pdf(dist, x(i));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_phexp(
  Rcpp::NumericVector q, Rcpp::NumericVector probs, Rcpp::NumericVector rates,
  bool lower
){
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(prbs, rts);
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
Rcpp::NumericVector rcpp_qhexp(
  Rcpp::NumericVector p, Rcpp::NumericVector probs, Rcpp::NumericVector rates,
  bool lower
){
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(prbs, rts);
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
double hexp_mean(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return mean(dist);
}

// [[Rcpp::export]]
double hexp_mode(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return mode(dist);
}

// [[Rcpp::export]]
double hexp_variance(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return variance(dist);
}

// [[Rcpp::export]]
double hexp_skewness(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return skewness(dist);
}

// [[Rcpp::export]]
double hexp_kurtosis(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double hexp_kurtosis_excess(
  Rcpp::NumericVector probs, Rcpp::NumericVector rates
) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  boost::math::hyperexponential dist(prbs, rts);
  return kurtosis_excess(dist);
}
