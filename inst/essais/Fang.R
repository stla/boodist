library(boodist)

dat <- iris[, 3:5]
N <- nrow(dat)
names(dat) <- c("x1", "x2", "Group")
Z <- model.matrix(~ 0 + Group, data=dat)
names(Z) <- c("z1", "z2", "z3")
dat <- cbind(dat, as.data.frame(Z))
dat$u <- rep(1, N)
groups <- split(dat, ~ Group)

# initialization
G <- 3L
gamma <- deltag <- rep(1, G)
Mu <- aggregate(cbind(x1, x2) ~ Group, FUN = mean, data = dat)
Mu <- cbind(Mu$x1, Mu$x2)
Beta <- matrix(0.05, nrow = G, ncol = 2L)
Delta <- vapply(groups, function(df) {
  Sigma <- cov(cbind(df$x1, df$x2))
  D <- det(Sigma)
  Sigma / D^(1/2)
}, matrix(NA_real_, nrow = 2L, ncol = 2L))
rho <- vapply(groups, nrow, integer(1L)) / nrow(dat)


f <- function(g) {
  df <- groups[[g]]
  n <- nrow(df)
  y <- cbind(df$x1, df$x2)
  mu <- Mu[g, ]
  Sigma <- Delta[, , g]
  Lambda <- solve(Sigma)
  beta <- Beta[g, ]
  u <- df$u
  Z <- vapply(1L:n, function(i) {
    v <- y[i, ] - mu - u[i] * Sigma %*% beta
    - t(v) %*% Lambda %*% v / (2*u[i])
  }, numeric(1L))
  1/det(Sigma)^(n/2) * exp(sum(Z))
}

PrZ <- rho * vapply(1L:G, f, numeric(1L))
Z <- rmultinom(N, 1L, PrZ)

Levels <- c("setosa", "versicolor", "virginica")
apply(Z, 2L, function(z) {
  Levels[which(z == 1L)]
})

