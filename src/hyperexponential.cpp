#include "boodist.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dhexp(
    Rcpp::NumericVector x, Rcpp::NumericVector probs, Rcpp::NumericVector rates
){
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double* p = prbs.data();
  const double* lambda = rts.data();
  int n = x.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(p, lambda);
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
  const double* p = prbs.data();
  const double* lambda = rts.data();
  int n = q.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(p, lambda);
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
  const double* pr = prbs.data();
  const double* lambda = rts.data();
  int n = p.size();
  Rcpp::NumericVector out(n);
  boost::math::hyperexponential dist(pr, lambda);
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
  const double* p = prbs.data();
  const double* lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return mean(dist);
}

// [[Rcpp::export]]
double hexp_mode(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double *p = prbs.data();
  const double *lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return mode(dist);
}

// [[Rcpp::export]]
double hexp_variance(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double* p = prbs.data();
  const double* lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return variance(dist);
}

// [[Rcpp::export]]
double hexp_skewness(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double* p = prbs.data();
  const double* lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return skewness(dist);
}

// [[Rcpp::export]]
double hexp_kurtosis(Rcpp::NumericVector probs, Rcpp::NumericVector rates) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double* p = prbs.data();
  const double* lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return kurtosis(dist);
}

// [[Rcpp::export]]
double hexp_kurtosis_excess(
  Rcpp::NumericVector probs, Rcpp::NumericVector rates
) {
  std::vector<double> prbs(probs.begin(), probs.end());
  std::vector<double> rts(rates.begin(), rates.end());
  const double* p = prbs.data();
  const double* lambda = rts.data();
  boost::math::hyperexponential dist(p, lambda);
  return kurtosis_excess(dist);
}
