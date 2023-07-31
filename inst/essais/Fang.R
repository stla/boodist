library(boodist)
library(TruncatedNormal)
library(mvtnorm)

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


a0 <- a0 + t0
a1 <- a1 + t1
a2 <- a2 + t2
a3 <- a3 + t3
a4 <- a4 + t4
a5 <- a5 + t5


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


