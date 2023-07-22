library(boodist)

he <- Hyperexponential$new(c(0.3, 0.7), c(0.2, 1))

x <- seq(0, 6, length.out = 200L)
y <- he$d(x)


opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.8),
  main = "Hyperexponential distribution"
)
axis(1L)


svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.8),
  main = "Hyperexponential distribution"
)
axis(1L)
dev.off()

rsvg::rsvg_png("x.svg", "hyperexponential.png", width = 512, height = 256)

