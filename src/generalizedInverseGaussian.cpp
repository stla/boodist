#include "boodist.h"

double psi(const double x, const double alpha, const double lambda) {
  //const double alpha = std::sqrt(omega*omega + lambda*lambda) - lambda;
  return -alpha*(std::cosh(x) - 1) - lambda*(std::exp(x) - x - 1);
}

std::pair<double,double> psipsiprime(
  const double x, const double alpha, const double lambda
) {
  const double expxm1 = std::expm1(x);
  const double psix = -alpha*(std::cosh(x) - 1) - lambda*(expxm1 - x);
  const double psiprimex = -alpha*std::sinh(x) - lambda*expxm1;
  return std::make_pair(psix, psiprimex);
}

double rgig_rcpp(
  const double x, const double omega, const double lambda
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
  if(mpsi1 < 0.5) {
    const double beta = 1.0 / alpha;
    s = std::min(1/lambda, std::log1p(beta + std::sqrt(beta*(beta + 2.0))));
  } else if(mpsi1 <= 2) {
    t = 1;
  } else {
    t = std::sqrt(2 / (alpha + lambda));
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
  //const double q = tprime + sprime;
  double out;
  if(x < sprime) {
    out = std::exp(-theta + xi*(x + s));
  } else if(x <= tprime) {
    out = 1.0;
  } else {
    out = std::exp(-eta - zeta*(x - t));
  }
  return out;
}
