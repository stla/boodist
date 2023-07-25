#include "boodist.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

double psi(const double x, const double alpha, const double lambda) {
  //const double alpha = std::sqrt(omega*omega + lambda*lambda) - lambda;
  return -alpha*(std::cosh(x) - 1) - lambda*(std::expm1(x) - x);
}

std::pair<double,double> psipsiprime(
  const double x, const double alpha, const double lambda
) {
  const double expxm1 = std::expm1(x);
  const double psix = -alpha*(std::cosh(x) - 1) - lambda*(expxm1 - x);
  const double psiprimex = -alpha*std::sinh(x) - lambda*expxm1;
  return std::make_pair(psix, psiprimex);
}

double chi(
    const double x,
    const double s, const double sprime, const double t, const double tprime,
    const double eta, const double zeta, const double theta, const double xi
) {
  double out;
  if(x < -sprime) {
    out = std::exp(-theta + xi*(x + s));
  } else if(x <= tprime) {
    out = 1.0;
  } else {
    out = std::exp(-eta - zeta*(x - t));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector rgig_rcpp(
  const unsigned n, const double lambda, const double omega
) {
  const double alpha = std::sqrt(omega*omega + lambda*lambda) - lambda;

  const double mpsi1  = - psi(1.0, alpha, lambda);
  const double mpsim1 = - psi(-1.0, alpha, lambda);
  double s, t;
  if(mpsi1 < 0.5) {
    t = std::log(4 / (alpha + lambda + lambda));
  } else if(mpsi1 <= 2) {
    t = 1;
  } else {
    t = std::sqrt(2 / (alpha + lambda));
  }
  if(mpsim1 < 0.5) {
    const double beta = 1.0 / alpha;
    s = std::min(1/lambda, std::log1p(beta + std::sqrt(beta*(beta + 2.0))));
  } else if(mpsim1 <= 2) {
    s = 1;
  } else {
    s = std::sqrt(2 / (alpha + lambda));
  }

  std::pair<double,double> ppps = psipsiprime(-s, alpha, lambda);
  std::pair<double,double> pppt = psipsiprime(t, alpha, lambda);
  const double eta   = -pppt.first;
  const double zeta  = -pppt.second;
  const double theta = -ppps.first;
  const double xi    = ppps.second;
  const double p = 1 / xi;
  const double r = 1 / zeta;
  const double tprime = t - r * eta;
  const double sprime = s - p * theta;
  const double q = tprime + sprime;
  const double invsumpqr = 1.0 / (p + q + r);
  const double u1 = q * invsumpqr;
  const double u2 = (q + r) * invsumpqr;
  const double lambda_over_omega = lambda / omega;
  const double outfactor =
    lambda_over_omega + std::sqrt(1 + lambda_over_omega*lambda_over_omega);

  Rcpp::NumericVector out(n);
  boost::mt19937 gen;
  boost::random::uniform_real_distribution<double> runif(0.0, 1.0);
  double x;
  double u, v, w;
  for(unsigned i = 0; i < n; i++) {
    do {
      u = runif(gen);
      v = runif(gen);
      w = runif(gen);
      if(u < u1) {
        x = -sprime + q*v;
      } else if(u < u2) {
        x = tprime - r*std::log(v);
      } else {
        x = -sprime + p*std::log(v);
      }
    } while(
        w*chi(x, s, sprime, t, tprime, eta, zeta, theta, xi) >
      std::exp(psi(x, alpha, lambda))
    );
    out(i) = outfactor * std::exp(x);
  }
  return out;
}

