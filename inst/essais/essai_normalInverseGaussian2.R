library(boodist)


nig <- NormalInverseGaussian$new(0, 2, 1, 2)

sims <- nig$r(55000L)

plot(ecdf(sims))
q_ <- seq(-2, 6, length.out = 100)
y <- nig$p(q_)
lines(q_, y, col = "red", lwd = 3, lty = "dashed")


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

