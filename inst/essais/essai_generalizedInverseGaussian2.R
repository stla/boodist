library(boodist)

theta <- 1
eta <- 2
lambda <- 2

gig <- GeneralizedInverseGaussian$new(theta, eta, lambda)

sims_gig <- gig$r(20000)

x <- seq(0, 20, length.out = 100L)
y <- gig$d(x)

plot(x, y, type = "l")
lines(density(sims_gig, n = 1000), col ="red")

mean(sims_gig)
gig$mean()
var(sims_gig)
gig$variance()

#### p & q ####
theta <- 2; eta <- 3; lambda <- 2
q <- 3
( p <- gig$p(q) )
gig$q(p)

