library(boodist)
library(TruncatedNormal)
library(mvtnorm)
library(gtools)
library(ggplot2)

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
  chol2inv(chol(x*H + W))
}


dat <- iris[, 3:5]
d <- 2L
N <- nrow(dat)
names(dat) <- c("x1", "x2", "Group")
# dat$x1 <- log(dat$x1)
# dat$x2 <- log(dat$x2)
dat$id <- as.character(1L:N)
Z <- model.matrix(~ 0 + Group, data = dat)
names(Z) <- c("z1", "z2", "z3")
dat <- cbind(dat, as.data.frame(Z))
dat$u <- rep(0.011, N)
groups <- split(dat, ~ Group)
Levels <- names(groups)
G <- length(Levels)

set.seed(666)
for(g in 1L:G) {
  df <- groups[[g]]
  sims <- rmvnorm(50L, mean = c(5*g, 5*g))
  df$x1 <- sims[, 1L]
  df$x2 <- sims[, 2L]
  groups[[g]] <- df
}
dat <- do.call(rbind, groups)
#ggplot(dat) + geom_point(aes(x = x1, y = x2, color = Group))

# prior parameters
a0 <- rep(1, G)
a1 <- matrix(1, nrow = G, ncol = d)
a2 <- matrix(1, nrow = G, ncol = d)
a3 <- rep(2, G)
a4 <- rep(1, G)
a5 <- array(1, dim = c(d, d, G))
nu0 <- d
Lambda0 <- diag(d)

# initialization
gamma <- delta <- rep(1, G)
Mu <- aggregate(cbind(x1, x2) ~ Group, FUN = mean, data = dat)
Mu <- cbind(Mu$x1, Mu$x2)
Beta <- matrix(0.05, nrow = G, ncol = 2L)
Delta <- vapply(groups, function(df) {
  Sigma <- cov(cbind(df$x1, df$x2))
  D <- det(Sigma)
  Sigma / D^(1/2)
}, matrix(NA_real_, nrow = 2L, ncol = 2L))
Lambda <- apply(Delta, 3, solve, simplify = FALSE)
rho <- vapply(groups, nrow, integer(1L)) / nrow(dat)


for(k in 1L:10L) {

  ig <- InverseGaussian$new(gamma[g], delta[g])
  # Gibbs step 1
  f <- function(g) {
    df <- groups[[g]]
    n <- nrow(df)
    y <- cbind(df$x1, df$x2)
    mu <- Mu[g, ]
    Sigma <- Delta[, , g]
    iSigma <- Lambda[[g]]
    beta <- Beta[g, ]
    u <- df$u
    Z <- vapply(1L:n, function(i) {
      v <- y[i, ] - mu - c(u[i] * Sigma %*% beta)
      #- t(v) %*% iSigma %*% v / (2*u[i])
      dmvnorm(t(v), sigma = 2*u[i]*Sigma, log = TRUE) +
        ig$d(u[i], log = TRUE)
    }, numeric(1L))
    #1/det(Sigma)^(n/2) * sum(Z)
    sum(Z)
  }
  logPrZ <- vapply(1L:G, f, numeric(1L))
  logPrZ <- logPrZ - max(logPrZ)
  PrZ <- rho * exp(logPrZ)
  Z <- rmultinom(N, 1L, PrZ)
  grp <- apply(Z, 2L, function(z) {
    Levels[which(z == 1L)]
  })
  gindices <- lapply(Levels, function(lvl) {
    which(grp == lvl)
  })

  # Gibbs step 2
  al <- function(g) {
    Sigma <- Delta[, , g]
    sqrt(gamma[g]^2 + t(Beta[g, ] %*% Sigma %*% Beta[g, ]))
  }
  alpha <- vapply(1:G, al, numeric(1L))

  # newu <- function(i) {
  #   y <- as.numeric(dat[i, c("x1", "x2")])
  #   g <- which(Levels == dat[i, "Group"])
  #   iSigma <- Lambda[[g]]
  #   mu <- Mu[g, ]
  #   q <- sqrt(delta[g]^2 + t(y-mu) %*% iSigma %*% (y-mu))
  #   alph <- alpha[g]
  #   GeneralizedInverseGaussian$new(q*alph, alph/q, (d+1)/2)$r(1L)
  # }
  # u <- vapply(1L:N, newu, numeric(1L))
  for(g in 1L:G) {
    df <- groups[[g]]
    iSigma <- Lambda[[g]]
    mu <- Mu[g, ]
    alph <- alpha[g]
    for(i in 1L:nrow(df)) {
      y <- as.numeric(df[i, c("x1", "x2")])
      q <- sqrt(delta[g]^2 + t(y-mu) %*% iSigma %*% (y-mu))
      df$u[i] <-
        GeneralizedInverseGaussian$new(q*alph, alph/q, (d+1)/2)$r(1L)
    }
    groups[[g]] <- df
  }

  # gindices...
  dat <- do.call(rbind, groups)
  groups <- split(dat, grp)
  #G <- length(groups)


  t0 <- rowSums(Z)
  t1 <- matrix(NA_real_, nrow = G, ncol = d)
  t2 <- matrix(NA_real_, nrow = G, ncol = d)
  t3 <- rep(NA_real_, G)
  t4 <- rep(NA_real_, G)
  t5 <- array(0, dim = c(d, d, G))
  S0 <- array(0, dim = c(d, d, G))
  for(g in 1L:G) {
    df <- groups[[g]]
    y <- as.matrix(df[, c("x1", "x2")])
    u <- df$u
    t1[g, ] <- colSums(y)
    t2[g, ] <- colSums(y / u)
    t3[g] <- sum(u) / 2
    t4[g] <- sum(1/u) / 2
    for(j in 1L:nrow(df)) {
      t5[, , g] <- t5[, , g] + crossprod(y[j, , drop = FALSE]) / u[j]
      S0[, , g] <- S0[, , g] + crossprod(y[j, , drop = FALSE] - Mu[g, ]) / u[j]
    }
    # i <- gindices[[g]]
    # t1[g, ] <- colSums(y[i, , drop = FALSE])
    # t2[g, ] <- colSums(y[i, , drop = FALSE] / u[i])
    # t3[g] <- sum(u[i]) / 2
    # t4[g] <- sum(1/u[i]) / 2
    # for(j in i) {
    #   t5[, , g] <- t5[, , g] + crossprod(y[j, , drop = FALSE]) / u[j]
    #   S0[, , g] <- S0[, , g] + crossprod(y[j, , drop = FALSE] - Mu[g, ]) / u[j]
    # }
  }

  a0 <- a0 + t0
  a1 <- a1 + t1
  a2 <- a2 + t2
  a3prev <- a3
  a3 <- a3 + t3
  a4 <- a4 + t4
  a5 <- a5 + t5 # not used!..

  delta <- sqrt(rgamma(G, shape = a0/2+1, rate = a4 - a0^2/(4*a3)))
  gamma <- vapply(1:G, function(g) {
    rtnorm(1, a0[g]*delta[g]/(2*a3[g]), 1/(2*a3[g]), lb = 0, ub = Inf)
  }, numeric(1L))

  h <- 4*a3*a4 - a0^2
  f <- function(g) {
    D <- Delta[, , g]
    Sigma11 <- 2*a3[g] * D / h[g]
    iD <- Lambda[[g]]
    Sigma22 <- 2*a4[g] * iD / h[g]
    Sigma12 <- - a0[g] * diag(d) / h[g]
    Sigma <- rbind(cbind(Sigma11, Sigma12), cbind(Sigma12, Sigma22))
    m1 <- (2*a2[g, ]*a3[g] - a0[g]*a1[g, ]) / h[g]
    m2 <- iD %*% (2*a1[g, ]*a4[g] - a0[g]*a2[g, ]) / h[g]
    m <- c(m1, m2)
    c(rmvnorm(1L, mean = m, sigma = Sigma))
  }
  MuBeta <- t(vapply(1L:G, f, numeric(2L*d)))
  Mu <- MuBeta[, 1L:d]
  Beta <- MuBeta[, -(1L:d)]

  Delta <- array(NA_real_, dim = c(d, d, G))
  for(g in 1L:G) {
    Delta[, , g] <- rMGIG(
      q = (nu0+t0[g])/2, z = Beta[g, ], A = Lambda0 + S0[, , g],
      b = 2*t3[g], d = d
    )
  }
  Lambda <- apply(Delta, 3, solve, simplify = FALSE)

  rho <- rdirichlet(1L, a0)[1L, ]
  o <- order(rho)
  rho <- rho[o]
  groups <- groups[o]
  Levels <- names(groups)
}
