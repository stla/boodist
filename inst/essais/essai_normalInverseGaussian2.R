library(boodist)


nig <- NormalInverseGaussian$new(0, 2, 1, 2)

( quants <- nig$q(c(0.2, 0.6), a = 0, b = 2) )
nig$p(quants)

sims <- nig$r(55000L)

plot(ecdf(sims))
q_ <- seq(-2, 6, length.out = 100)
y <- nig$p(q_)
lines(q_, y, col = "red", lwd = 3, lty = "dashed")


q <- nig$mean() + 6*nig$sd()
nig$p(q)

mn <- nig$mean()
std <- nig$sd()
q_ <- seq(mn - 6*std, mn + 6*std, length.out = 200L)
prob_ <- nig$p(q_)
f <- approxfun(prob_, q_)
f(0.6)


p <- 0.6

f <- function(aq) {
  nig$p(tan(aq)) - 0.6
}

uroot <- uniroot(f, lower = -pi/2+0.01, upper = pi/2-0.01,
                 tol = .Machine$double.eps^0.5)
root <- uroot$root
q <- tan(root)

nig$p(q)

boodist:::qnig_rcpp(0.6, 0, 2, 1, 2)

