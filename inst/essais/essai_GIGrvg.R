library(boodist)
library(GIGrvg)

lambda <- 2
chi <- 3
psi <- 4
x <- 2

dgig(x, lambda, chi, psi)

theta <- sqrt(chi * psi)
eta <- sqrt(chi / psi)
gig <- GeneralizedInverseGaussian$new(theta, eta, lambda)
gig$d(x)
