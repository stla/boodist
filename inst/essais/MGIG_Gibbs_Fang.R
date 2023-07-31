library(boodist)
library(EigenR)

q <- 3
z <- c(2, 3)
A <- toeplitz(c(3, 1))
d <- 2
q <- 3
p <- q + (1-d)/2
b <- 2
H <- tcrossprod(z)


theta <- sqrt(c(t(z) %*% A %*% z * b))
eta <- sqrt(b) / sqrt(c(t(z) %*% A %*% z))

gig <- GeneralizedInverseGaussian$new(
  theta = theta, eta = eta, lambda = -p)
gig$r(1L)


nsims <- 30000L
sims <- array(NA_real_, dim = c(d, d, nsims))
x <- gig$r(nsims)
W <- rWishart(nsims, q, A)
for(i in 1L:nsims) {
  # x <- gig$r(1L)
  # W <- rWishart(1L, q, A)[, , 1L]
  sims[, , i] <- chol2inv(chol(x[i]*H + W[, , i]))
}
sims <- sims[, , -(1L:2000L)]

q*A + gig$mean() * H
apply(sims, 1:2, mean)

inv <- apply(sims, 3, solve, simplify = TRUE)
apply(inv, 1, mean)


dets <- numeric(dim(sims)[3L])
for(i in seq_along(dets)) {
  dets[i] <- det(sims[, , i])
}
plot(dets)
mean(dets)

Psi <- b*H
Ga <- A
la <- -q

Bessel2(Psi %*% Ga / 4, la + 1 + 3/2) / Bessel2(Psi %*% Ga / 4, la +3/2) / det(Ga / 2)
Bessel2(Psi %*% Ga / 4, la + 1) / Bessel2(Psi %*% Ga / 4, la) / det(Ga / 2)

q*A + gig$mean() * H

inv <- apply(sims, 3, solve, simplify = TRUE)
apply(inv, 1, mean)



Bessel1 <- function(t, nu) {
  2 * besselK(2*sqrt(t), nu) / t^(nu/2)
}
Bessel2 <- function(Z, delta) {
  trZ <- sum(diag(Z))
  detZ <- det(Z)
  f <- function(x) {
    sqrt(pi)*exp(-x)*exp(-trZ/x)*Bessel1(detZ/x^2, delta)/x^(delta+3/2)
  }
  integrate(f, 0, Inf)$value
}



rMGIG <- function(q, z, A, b, d) {
  p <- q + (1-d)/2
  H <- tcrossprod(z)
  theta <- sqrt(c(t(z) %*% A %*% z * b))
  eta <- sqrt(b) / sqrt(c(t(z) %*% A %*% z))
  gig <- GeneralizedInverseGaussian$new(
    theta = theta, eta = eta, lambda = -p
  )
  x <- gig$r(1L)
  W <- rWishart(1L, q, A)[, , 1L]
  chol2inv(chol(x[i]*H + W[, , i]))
}
