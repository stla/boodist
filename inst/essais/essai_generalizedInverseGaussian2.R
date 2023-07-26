library(boodist)

theta <- 2
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

# picture ####
theta <- 1
eta <- 2
lambda <- 2

gig <- GeneralizedInverseGaussian$new(theta, eta, lambda)

x <- seq(0, 25, length.out = 250L)
y <- gig$d(x)

opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Generalized inverse Gaussian distribution"
)
axis(1L)

svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Generalized inverse Gaussian distribution"
)
axis(1L)
dev.off()
rsvg::rsvg_png("x.svg", "generalizedInverseGaussian.png", width = 512, height = 256)

