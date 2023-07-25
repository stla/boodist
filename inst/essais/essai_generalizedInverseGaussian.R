lambda <- 3
omega <- 1
dgig <- function(x) {
 1 / 2 / besselK(1, lambda) * x^(lambda-1) * exp(-(x+1/x)/2)
}
sims_gig <- boodist:::rgig_rcpp(40000, lambda, omega)

x <- seq(0, 20, length.out = 100L)
y <- dgig(x)

plot(x, y, type = "l")
lines(density(sims_gig, n = 1000), col ="red")

#### p & q ####
theta <- 2; eta <- 3; lambda <- 2
q <- 3
( p <- boodist:::pgig_rcpp(q, theta, eta, lambda) )
boodist:::qgig_rcpp(p, 2, 4, theta, eta, lambda)
