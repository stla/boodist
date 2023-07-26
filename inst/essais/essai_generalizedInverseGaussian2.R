library(boodist)

theta <- 2
eta <- 2
lambda <- 2

gig <- GeneralizedInverseGaussian$new(theta, eta, lambda)

sims_gig <- gig$r(20000)

x <- seq(0, 20, length.out = 100L)
y <- gig$d(x)

plot(x, y, type = "l")
lines(density(sims_gig, n = 1000), col ="red")

mean(sims_gig)
gig$mean()
var(sims_gig)
gig$variance()

#### p & q ####
theta <- 2; eta <- 3; lambda <- 2
q <- 3
( p <- gig$p(q) )
gig$q(p)

# picture ####
theta <- 1
eta <- 2
lambda <- 2

gig <- GeneralizedInverseGaussian$new(theta, eta, lambda)

x <- seq(0, 25, length.out = 250L)
y <- gig$d(x)

opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Generalized inverse Gaussian distribution"
)
axis(1L)

svg("x.svg", width = 8, height = 4)
opar <- par(mar = c(3, 1, 3, 1))
plot(
  x, y, type = "l", xlab = NA, ylab = NA, axes = FALSE,
  lwd = 4, col = "blue", xaxs = "i", yaxs = "i", ylim = c(0, 0.1),
  main = "Generalized inverse Gaussian distribution"
)
axis(1L)
dev.off()
rsvg::rsvg_png("x.svg", "generalizedInverseGaussian.png", width = 512, height = 256)

# plotly ####
library(plotly)
library(boodist)

x_ <- seq(0, 3, length.out = 100L)
lambda_ <- seq(-1, 1, length.out = 100L)
dsty <- vapply(lambda_, function(lambda) {
  GeneralizedInverseGaussian$new(theta = 1, eta = 1, lambda)$d(x_)
}, numeric(length(x_)))
#
txt <- matrix(NA_character_, nrow = length(x_), ncol = length(lambda_))
for(i in 1L:nrow(txt)) {
  for(j in 1L:ncol(txt)) {
    txt[i, j] <- paste0(
      "x: ", formatC(x_[i]),
      "<br> lambda: ", formatC(lambda_[j]),
      "<br> density: ", formatC(dsty[i, j])
    )
  }
}
#
plot_ly(
  x = ~lambda_, y = ~x_, z = ~dsty, type = "surface",
  text = txt, hoverinfo = "text", showscale = FALSE
) %>% layout(
  title = "Generalized inverse Gaussian distribution",
  margin = list(t = 40, r= 5, b = 5, l = 5),
  scene = list(
    xaxis = list(
      title = "lambda"
    ),
    yaxis = list(
      title = "x"
    ),
    zaxis = list(
      title = "density"
    )
  )
)

