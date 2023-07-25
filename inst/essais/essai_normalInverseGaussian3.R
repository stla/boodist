library(boodist)

nig <- NormalInverseGaussian$new(0, 2, 1, 2)
p_ <- seq(0.1, 0.9, length.out = 10L)
q_ <- nig$q(p_)
nig$p(q_) - p_


