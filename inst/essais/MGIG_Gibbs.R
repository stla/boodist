###----------------------------------------------------###
###----------------------------------------------------###
###       Gibbs sampler for Matrix Generalized         ###
###       Inverse Gaussian (MGIG) distribution         ###
###----------------------------------------------------###
###----------------------------------------------------###


###---------------------------------------###
###     (proposed)   Gibbs sampler        ###
###---------------------------------------###
## INPUT
# la: shape parameter (scalar)
# Psi: shape parameter (matrix)
# Ga: shape parameter (matrix)
# mc: length of Gibbs sampler
# burn: burn-in period

## OUTPUT
# random samples from MGIG

library(boodist)

rMGIG_GS <- function(la, Psi, Ga, mc = 1000L, burn = mc/2L, print = TRUE){

  # dimension
  p <- nrow(Psi)

  # initial values
  a <- rep(1, p)
  B <- diag(p)

  # arrays to store random samples
  Si_rn <- array(NA_real_, dim = c(p, p, mc))

  # Gibbs sampler (when p>1)

  for(item in 1L:mc) {
    # a
    inv_B <- solve(B)
    Ga_tilde <- inv_B %*% Ga %*% t(inv_B)
    Chi_a <- diag(Ga_tilde)
    Psi_a <- diag(t(B) %*% Psi %*% B)
    Lambda_a <- la + p - (1L:p) + 1
    Theta <- sqrt(Chi_a * Psi_a)
    Eta <- sqrt(Chi_a / Psi_a)
    a <- numeric(p)
    for(i in 1L:p) {
      gig <- GeneralizedInverseGaussian$new(
        theta = Theta[i], eta = Eta[i], lambda = Lambda_a[i])
      a[i] <- gig$r(1L)
    }

    # b (i = 1)
    P <- crossprod(inv_B, inv_B / a)
    M1 <- Psi
    R1 <- B
    R1[2L:p, 1L] <- 0
    R1 <- t(sqrt(a) * t(R1))
    M2 <- Ga
    R2 <- t(inv_B)
    R2[1L, ] <-
      R2[1L, , drop = FALSE] + t(B[2L:p, 1L]) %*% R2[2L:p, , drop = FALSE]
    R2 <- t(t(R2) / sqrt(a))
    vn <- - M1[2L:p, 1L:p, drop = FALSE] %*% (R1 %*% R1[1L, ]) +
      R2[2L:p, , drop = FALSE] %*% t(M2[1L, , drop = FALSE] %*% R2)
    N <- a[1L] * Psi[2L:p, 2L:p, drop = FALSE] +
      M2[1L, 1L] * P[2L:p, 2L:p, drop = FALSE]
    chol_N <- chol(N)
    z <- rnorm(p - 1L) + backsolve(r = t(chol_N), x = vn, upper.tri = FALSE)
    B[2L:p, 1L] <- backsolve(r = chol_N, x = z, upper.tri = TRUE)
    #  b (i > 1)
    for(i in seq(from = 2L, length.out = p - 2L)) {
      M1[, i - 1L] <- M1[, i - 1L] + M1[, i:p, drop = FALSE] %*% B[i:p, i - 1L]
      M1[i - 1L, ] <- M1[i - 1L, , drop = FALSE] +
        t(B[i:p, i - 1L]) %*% M1[i:p, , drop = FALSE]
      R1[(i + 1L):p, ] <- R1[(i + 1L):p, , drop = FALSE] -
        outer(B[(i + 1L):p, i], R1[i, ])
      M2[, i:p] <- M2[, i:p, drop = FALSE] - outer(M2[, i - 1L], B[i:p, i - 1L])
      M2[i:p, ] <- M2[i:p, , drop = FALSE] - outer(B[i:p, i - 1L], M2[i - 1L, ])
      R2[i, ] <- R2[i, , drop = FALSE] +
        t(B[(i + 1L):p, i]) %*% R2[(i + 1L):p, , drop = FALSE]
      vn <- - M1[(i + 1L):p, 1L:p, drop = FALSE] %*% (R1 %*% R1[i, ]) +
        R2[(i + 1L):p, , drop = FALSE] %*% t(M2[i, , drop = FALSE] %*% R2)
      N <- a[i] * Psi[(i + 1L):p, (i + 1L):p, drop = FALSE] +
        M2[i, i] * P[(i + 1L):p, (i + 1L):p, drop = FALSE]
      chol_N <- chol(N)
      z <- rnorm(p - i) + backsolve(r = t(chol_N), x = vn, upper.tri = FALSE)
      B[(i + 1L):p, i] <- backsolve(r = chol_N, x = z, upper.tri = TRUE)
    }

    # Si
    t_sqrt_Si <- sqrt(a) * t(B)
    Si_rn[, , item] <- crossprod(t_sqrt_Si)

    # print
    if(print && item %% 1000L == 0L){ print(item) }
  }

  # output
  return( Si_rn[, , -(1L:burn), drop = FALSE] )
}


la <- 3
Psi <- diag(c(4,5))#rWishart(1L, df = 6, Sigma = diag(2L))[, , 1L]
Ga <- toeplitz(4:3)
sims <- rMGIG_GS(la, Psi, Ga, mc = 20000L, burn = 1000L)

sims11 <- sims[1L, 1L, ]
plot(sims11, type = "l")

plot(density(sims11))
mean(sims11)

dets <- numeric(dim(sims)[3L])
for(i in seq_along(dets)) {
  dets[i] <- det(sims[, , i])
}
plot(dets)
mean(dets)

Bessel2(Psi %*% Ga / 4, la + 1 + 3/2) / Bessel2(Psi %*% Ga / 4, la +3/2) / det(Ga / 2)
Bessel2(Psi %*% Ga / 4, la + 1) / Bessel2(Psi %*% Ga / 4, la) / det(Ga / 2)

library(HypergeoMat)

BesselA(20L, Psi %*% Ga / 4, la + 1) / BesselA(20L, Psi %*% Ga / 4, la) / sqrt(det(Psi / 2))


Bessel1 <- function(t, nu) {
  2 * besselK(2*sqrt(t), nu) / t^(nu/2) # je trouve K(2*sqrt(t))
}

t <- 4
nu <- 3
Bessel1(t, nu)
Bessel1(t, -nu) / t^nu

lambda <- 0.1
delta <- 1/2
BesselA(20L, lambda, delta) - BesselA(20L, -lambda, -delta) /lambda^delta
-sin(pi*delta)/pi * Bessel1(lambda, delta)


Bessel2 <- function(Z, delta) {
  trZ <- sum(diag(Z))
  detZ <- det(Z)
  f <- function(x) {
    sqrt(pi)*exp(-x)*exp(-trZ/x)*Bessel1(detZ/x^2, delta)/x^(delta+3/2)
  }
  integrate(f, 0, Inf)$value
}

Z <- toeplitz(2:1)
Bessel2(Z, 3)
