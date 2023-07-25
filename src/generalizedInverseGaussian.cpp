#include "boodist.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

double psi(const double x, const double alpha, const double lambda) {
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
  const unsigned n, const double lambda, const double omega // omega <=> theta
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

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
double dgig(double x, double theta, double eta, double lambda) {
  double y = x / eta;
  return 1.0 / (2.0 * eta * cyl_bessel_k(lambda, theta)) *
    std::pow(y, lambda-1) * std::exp(-theta * (y + 1.0/y) / 2.0);
}

class GIGpdf : public Func {
private:
  double theta;
  double eta;
  double lambda;

public:
  GIGpdf(double theta_, double eta_, double lambda_)
    : theta(theta_), eta(eta_), lambda(lambda_) {}

  double operator()(const double& x) const {
    return dgig(x, theta, eta, lambda);
  }
};


// [[Rcpp::export]]
Rcpp::NumericVector pgig_rcpp(Rcpp::NumericVector q,
                              const double theta,
                              const double eta,
                              const double lambda) {
  GIGpdf f(theta, eta, lambda);
  const int subdiv = 150;
  const double eps_abs = 1e-8;
  const double eps_rel = 1e-8;
  int n = q.size();
  Rcpp::NumericVector out(n);
  Rcpp::NumericVector error_estimate(n);
  Rcpp::IntegerVector error_code(n);
  for(int i = 0; i < n; i++) {
    double err_est;
    int err_code;
    const double upper = q(i);
    const double res =
      integrate(f, 0.0, upper, err_est, err_code, subdiv, eps_abs, eps_rel,
                Integrator<double>::GaussKronrod201);
    out(i) = res;
    error_estimate(i) = err_est;
    error_code(i) = err_code;
    if(err_code != 0) {
      Rcpp::warning("An anomaly occured (see the error codes).");
    }
  }
  out.attr("error_estimate") = error_estimate;
  out.attr("error_code") = error_code;
  return out;
}



// [[Rcpp::export]]
Rcpp::NumericVector qgig_rcpp(
    Rcpp::NumericVector p, Rcpp::NumericVector g_a, Rcpp::NumericVector g_b,
    const double theta, const double eta, const double lambda
) {

  auto pdf = [theta, eta, lambda](double x) {
    return dgig(x, theta, eta, lambda);
  };

  int n = p.size();
  Rcpp::NumericVector out(n);

  for(int i = 0; i < n; i++) {
    double prob = p(i);
    auto integral = [pdf, prob](double f_q) {
      double error;
      return gauss_kronrod<double, 61>::integrate(
          pdf, 0.0, -std::log1p(-f_q), 15, 1e-6, &error
      ) - prob;
    };
    const double a = -std::expm1(-g_a(i));
    const double b = -std::expm1(-g_b(i));
    std::uintmax_t max_iter = 300;
    std::pair<double, double> interval = toms748_solve(
      integral, a, b,
      [](double l, double r){return fabs(l-r) < 1e-6;},
      max_iter
    );
    if(max_iter >= 300) {
      Rcpp::warning("Reached maximum number of iterations.");
    }
    out(i) = -(std::log1p(-interval.first) + std::log1p(-interval.second)) / 2;
  }

  return out;
}
