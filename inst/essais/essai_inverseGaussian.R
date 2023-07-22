library(boodist)

ig <- InverseGaussian$new(2, 1)

x <- seq(0, 5, length.out = 200L)
y <- ig$d(x)


opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.76),
  main = "Inverse Gaussian distribution"
)
axis(1L)

sims <- ig$r(45000L)
lines(density(sims))

svg("x.svg")
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.8),
  main = "Hyperexponential distribution"
)
axis(1L)
dev.off()

rsvg::rsvg_png("x.svg", "hyperexponential.png", width = 512, height = 512)

