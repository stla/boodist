library(plotly)
library(boodist)

dinvgamma <- function(x_, alpha_, beta = 1) {
  vapply(alpha_, function(alpha) {
    InverseGamma$new(alpha, beta)$d(x_)
  }, numeric(length(x_)))
}

x_ <- seq(0, 2, length.out = 100L)
alpha_ <- seq(0.5, 2.5, length.out = 100L)
dsty <- vapply(alpha_, function(alpha) {
  InverseGamma$new(alpha, beta = 1)$d(x_)
}, numeric(length(x_)))
dsty[1, ] <- 0
#
txt <- matrix(NA_character_, nrow = length(x_), ncol = length(alpha_))
for(i in 1L:nrow(txt)) {
  for(j in 1L:ncol(txt)) {
    txt[i, j] <- paste0(
      "x: ", formatC(x_[i]),
      "<br> alpha: ", formatC(alpha_[j]),
      "<br> density: ", formatC(dsty[i, j])
    )
  }
}
#
plot_ly(
  x = ~alpha_, y = ~x_, z = ~dsty, type = "surface",
  text = txt, hoverinfo = "text", showscale = FALSE
) %>% layout(
  title = "Inverse Gamma distribution",
  margin = list(t = 40, r= 5, b = 5, l = 5),
  scene = list(
    xaxis = list(
      title = "alpha"
    ),
    yaxis = list(
      title = "x"
    ),
    zaxis = list(
      title = "density"
    )
  )
)



dat$p <- 100*dat$p
gg <- ggplot(data=dat, aes(x=Delta, y=p)) +
  geom_line(size=1.5) +
  geom_hline(yintercept=p0, colour="blue", linetype = "dashed", size=1.2) +
  #xlab(GreekLetter("Delta")) + ylab("<i>p</i>") +
  xlab(sprintf("acceptable cumulative uncertainty %s", GreekLetter("Delta"))) +
  ylab("probability (%)") +
  theme(axis.title.x = element_text(size = rel(2))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(plot.title = element_text(size=21, face="italic")) +
  # ggtitle(sprintf("%s=%s, %s=%s",
  #                 GreekLetter("sigma"), formatC(sigma),
  #                 GreekLetter("delta"), formatC(delta)))
  ggtitle(sprintf("precision measure %s=%s, accuracy measure %s=%s",
                  GreekLetter("sigma"), formatC(sigma),
                  GreekLetter("delta"), formatC(delta)))
pgg <- ggplotly(gg, width=850, height=500)
pgg$x$data[[1]]$text <- paste(paste0(GreekLetter("Delta"), ":"), formatC(Delta),
                              '<br>p: ', formatC(100*p))
