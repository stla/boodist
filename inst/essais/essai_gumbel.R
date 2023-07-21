library(boodist)

gb <- Gumbel$new(3, 4)

x <- seq(-5, 25, length.out = 400L)
y <- gb$d(x)


opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Gumbel(3,4)"
)
axis(1L)


svg("x.svg")
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Gumbel(3,4)"
)
axis(1L)
dev.off()

rsvg::rsvg_png("x.svg", "Gumbel.png", width = 512, height = 512)

