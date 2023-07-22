library(boodist)

sn <- SkewNormal$new(2, 1, -5)

x <- seq(-2, 3, length.out = 200L)
y <- sn$d(x)


opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.75),
  main = "Skew normal distribution"
)
axis(1L)

sims <- sn$r(45000L)
lines(density(-sims+4), col="red", lwd = 3)

svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.75),
  main = "Skew normal distribution"
)
axis(1L)
dev.off()

rsvg::rsvg_png("x.svg", "skewNormal.png", width = 512, height = 256)

