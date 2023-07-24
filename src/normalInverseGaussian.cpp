#include <RcppNumerical.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
using namespace boost::math::quadrature;
#include <boost/math/tools/roots.hpp>
using namespace boost::math::tools;
#include "boodist.h"
using namespace Numer;
#include <cmath>

double dnig(double x, double mu, double alpha, double beta, double delta) {
  double gamm = std::sqrt(alpha * alpha - beta * beta);
  double s = std::sqrt(delta * delta + (x - mu) * (x - mu));
  return alpha * delta * cyl_bessel_k(1, alpha * s) *
         exp(delta * gamm + beta * (x - mu)) / (M_PI * s);
}

class NIGpdf : public Func {
 private:
  double mu;
  double alpha;
  double beta;
  double delta;

 public:
  NIGpdf(double mu_, double alpha_, double beta_, double delta_)
      : mu(mu_), alpha(alpha_), beta(beta_), delta(delta_) {}

  double operator()(const double& x) const {
    double gamm = std::sqrt(alpha * alpha - beta * beta);
    double s = std::sqrt(delta * delta + (x - mu) * (x - mu));
    return alpha * delta * cyl_bessel_k(1, alpha * s) *
           exp(delta * gamm + beta * (x - mu)) / (M_PI * s);
  }
};

// [[Rcpp::export]]
Rcpp::NumericVector pnig_rcpp(Rcpp::NumericVector q,
                              const double mu,
                              const double alpha,
                              const double beta,
                              const double delta) {
  const double lower = -std::numeric_limits<double>::infinity();
  NIGpdf f(mu, alpha, beta, delta);
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
        integrate(f, lower, upper, err_est, err_code, subdiv, eps_abs, eps_rel,
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
Rcpp::NumericVector qnig_rcpp(double p,
                              const double mu,
                              const double alpha,
                              const double beta,
                              const double delta) {

  auto pdf = [mu, alpha, beta, delta](double x) {
    return dnig(x, mu, alpha, beta, delta)
  };

  const double lower = -std::numeric_limits<double>::infinity();

  auto integral = [pdf, lower, p](double atanq) {
    double error;
    return gauss_kronrod<double, 61>::integrate(
                  pdf, lower, std::tan(atanq), 0, 0, &error
                ) - p;
  };

  const double a = -M_PI + 0.01;
  const double b = M_PI - 0.01;
  std::uintmax_t max_iter = 300;
  std::pair<double, double> interval = toms748_solve(
    integral,
    a, b,
    [](double l, double r){return abs(l-r) < 1e-6;},
    max_iter
  );
  // if(err_code != 0) {
  //   Rcpp::warning("An anomaly occured (see the error codes).");
  // }
  return Rcpp::NumericVector::create(
    std::tan(interval.first), std::tan(interval.second)
  );
}
