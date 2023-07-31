library(boodist)

dat <- iris[, 3:5]
d <- 2L
N <- nrow(dat)
names(dat) <- c("x1", "x2", "Group")
Z <- model.matrix(~ 0 + Group, data=dat)
names(Z) <- c("z1", "z2", "z3")
dat <- cbind(dat, as.data.frame(Z))
dat$u <- rep(1, N)
groups <- split(dat, ~ Group)
G <- 3L

# prior parameters
a0 <- rep(1, G)
a1 <- matrix(1, nrow = G, ncol = d)
a2 <- matrix(1, nrow = G, ncol = d)
a3 <- rep(1, G)
a4 <- rep(1, G)
a5 <- array(1, dim = c(d, d, G))

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
rho <- vapply(groups, nrow, integer(1L)) / nrow(dat)


Lambda <- apply(Delta, 3, solve, simplify = FALSE)

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
    v <- y[i, ] - mu - u[i] * Sigma %*% beta
    - t(v) %*% iSigma %*% v / (2*u[i])
  }, numeric(1L))
  1/det(Sigma)^(n/2) * exp(sum(Z))
}

PrZ <- rho * vapply(1L:G, f, numeric(1L))
Z <- rmultinom(N, 1L, PrZ)

Levels <- c("setosa", "versicolor", "virginica")
grp <- apply(Z, 2L, function(z) {
  Levels[which(z == 1L)]
})
gindices <- lapply(Levels, function(lvl) {
  which(grp == lvl)
})

#
al <- function(g) {
  Sigma <- Delta[, , g]
  sqrt(gamma[g]^2 + t(Beta[g, ] %*% Sigma %*% Beta[g, ]))
}
alpha <- vapply(1:G, al, numeric(1L))

newu <- function(i) {
  y <- as.numeric(dat[i, c("x1", "x2")])
  g <- which(Levels == dat[i, "Group"])
  iSigma <- Lambda[[g]]
  mu <- Mu[g, ]
  q <- sqrt(delta[g]^2 + t(y-mu) %*% iSigma %*% (y-mu))
  alpha <- alpha[g]
  GeneralizedInverseGaussian$new(q*alpha, alpha/q, (d+1)/2)$r(1L)
}

u <- vapply(1L:N, newu, numeric(1L))

t0 <- rowSums(Z)
t1 <- matrix(NA_real_, nrow = G, ncol = d)
t2 <- matrix(NA_real_, nrow = G, ncol = d)
t3 <- rep(NA_real_, G)
t4 <- rep(NA_real_, G)
t5 <- array(0, dim = c(d, d, G))

y <- cbind(dat$x1, dat$x2)

for(g in 1L:G) {
  i <- gindices[[g]]
  t1[g, ] <- colSums(y[i, , drop = FALSE])
  t2[g, ] <- colSums(y[i, , drop = FALSE] / u[i])
  t3[g] <- sum(u[i]) / 2
  t4[g] <- sum(1/u[i]) / 2
  for(j in i) {
    t5[, , g] <- t5[, , g] + crossprod(y[j, , drop = FALSE]) / u[j]
  }
}
