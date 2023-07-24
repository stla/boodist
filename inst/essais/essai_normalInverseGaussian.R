library(boodist)


nig <- NormalInverseGaussian$new(0, 2, 1, 2)

x <- seq(-2, 6, length.out = 200L)
y <- nig$d(x)


opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.39),
  main = "Normal-inverse Gaussian distribution"
)
axis(1L)

sims <- nig$r(55000L)
lines(density(sims), col = "red", lwd = 3)

plot(ecdf(sims))
q_ <- seq(-2, 6, length.out = 100)
yy <- vapply(q_, function(q) {
  ll <- boodist:::pnig_rcpp(q, 0, 2, 1, 2)
  ll$result
}, numeric(1L))
lines(q_, yy, col = "red", lwd = 3)


svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.39),
  main = "Normal-inverse Gaussian distribution"
)
axis(1L)
dev.off()

rsvg::rsvg_png("x.svg", "normalInverseGaussian.png", width = 512, height = 256)

