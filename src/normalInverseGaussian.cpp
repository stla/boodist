#include "boodist.h"
#include <boost/math/special_functions/bessel.hpp>
#include <RcppNumerical.h>
using namespace Numer;
#include <cmath>

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
Rcpp::List pnig_rcpp(const double q,
                     const double mu,
                     const double alpha,
                     const double beta,
                     const double delta) {
  const double lower = -std::numeric_limits<double>::infinity();
  NIGpdf f(mu, alpha, beta, delta);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, q, err_est, err_code);
  return Rcpp::List::create(Rcpp::Named("result") = res,
                            Rcpp::Named("error_estimate") = err_est,
                            Rcpp::Named("error_code") = err_code);
}
