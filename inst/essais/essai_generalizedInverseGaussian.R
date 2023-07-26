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

# general rgig ####
theta <- 2; eta <- 3; lambda <- 2
sims_gig <- boodist:::rgig_rcpp(10000, lambda, theta) * eta
plot(ecdf(sims_gig))
x <- seq(0, 25, length.out = 100L)
y <- boodist:::pgig_rcpp(x, theta, eta, lambda)
lines(x, y, lwd = 3, lty = "dashed", col = "red")


