#include "boodist.h"

// [[Rcpp::export]]
double find_chisq_ncp(double nu, double q, double p){
  // gives delta such that pchisq(q, nu, delta) = p
  return non_central_chi_squared::find_non_centrality(nu, q, p);
}

// [[Rcpp::export]]
double find_chisq_df(double delta, double q, double p){
  // gives nu such that pchisq(q, nu, delta) = p
  return non_central_chi_squared::find_degrees_of_freedom(delta, q, p);
}
