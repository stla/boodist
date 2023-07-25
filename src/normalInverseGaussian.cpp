#include "boodist.h"

double dnig(double x, double mu, double alpha, double beta, double delta) {
  double gamm = std::sqrt(alpha * alpha - beta * beta);
  double s = std::sqrt(delta * delta + (x - mu) * (x - mu));
  return alpha * delta * cyl_bessel_k(1, alpha * s) *
         std::exp(delta * gamm + beta * (x - mu)) / (M_PI * s);
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
    return dnig(x, mu, alpha, beta, delta);
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
Rcpp::NumericVector qnig_rcpp(
  Rcpp::NumericVector p, Rcpp::NumericVector tan_a, Rcpp::NumericVector tan_b,
  const double mu, const double alpha, const double beta, const double delta
) {

  auto pdf = [mu, alpha, beta, delta](double x) {
    return dnig(x, mu, alpha, beta, delta);
  };

  const double lower = -std::numeric_limits<double>::infinity();

  int n = p.size();
  Rcpp::NumericVector out(n);

  for(int i = 0; i < n; i++) {
    double prob = p(i);
    auto integral = [pdf, lower, prob](double atanq) {
      double error;
      return gauss_kronrod<double, 61>::integrate(
          pdf, lower, std::tan(atanq), 15, 1e-6, &error
      ) - prob;
    };
    const double a = std::atan(tan_a(i));
    const double b = std::atan(tan_b(i));
    std::uintmax_t max_iter = 300;
    std::pair<double, double> interval = toms748_solve(
      integral, a, b,
      [](double l, double r){return fabs(l-r) < 1e-6;},
      max_iter
    );
    if(max_iter >= 300) {
      Rcpp::warning("Reached maximum number of iterations.");
    }
    out(i) = (std::tan(interval.first) + std::tan(interval.second)) / 2;
  }

  return out;
}
