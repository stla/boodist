###---------------------------------------###
###        Hit-and-Run sampler            ###
###---------------------------------------###
rMGIG_HR <- function(ka, Psi, Phi, Si0=NULL, mc=1000, burn=mc/2, print=T){
  # function
  fun <- function(d){
    p <- length(d)
    d <- sort(d, decreasing=T)
    Mat_1 <- matrix(d, p, p, byrow=F)
    Mat_2 <- matrix(d, p, p, byrow=T)
    output <- sum( log((Mat_1-Mat_2)[upper.tri(Mat_1, diag=F)]) )
    return(output)
  }

  # preparation
  p <- dim(Psi)[1]
  L <- matrix(0, p, p)
  Lower <- lower.tri(L, diag = TRUE)

  # initial values
  if(is.null(Si0)){
    Si <- diag(seq(0, 2, length.out=p+1)[-1], nrow=p)
  }else{
    Si <- Si0
  }
  eigen_Si <- eigen(Si)
  eigen_Si_vectors <- eigen_Si$vectors
  eigen_Si_values <- eigen_Si$values
  inv_Si <- eigen_Si_vectors%*%diag(1/eigen_Si_values, nrow=p)%*%t(eigen_Si_vectors)
  eigen_log_Si_values <- log(eigen_Si_values)

  # arrays to store random samples
  Si_rn <- array(NA, dim = c(p, p, mc))

  # iterations
  for (item in 1:mc) {
    l <- rnorm(p*(p+1)/2, mean=0, sd=1)
    L[Lower] <- l
    s <- sqrt( sum(l^2) )
    tilde_L <- L/s
    D <- tilde_L + t(tilde_L) - diag(diag(tilde_L), nrow=p)
    la <- rnorm(1, mean=0, sd=1)
    log_proposal <- eigen_Si_vectors%*%diag(eigen_log_Si_values, nrow=p)%*%t(eigen_Si_vectors) + la*D
    eigen_log_proposal <- eigen(log_proposal)
    eigen_log_proposal_vectors <- eigen_log_proposal$vectors
    eigen_log_proposal_values <- eigen_log_proposal$values
    eigen_proposal_values <- exp(eigen_log_proposal_values)
    proposal <- eigen_log_proposal_vectors%*%diag(eigen_proposal_values, nrow=p)%*%t(eigen_log_proposal_vectors)
    inv_proposal <- eigen_log_proposal_vectors%*%diag(1/eigen_proposal_values, nrow=p)%*%t(eigen_log_proposal_vectors)
    log_ratio <- sum(ka*(eigen_log_proposal_values - eigen_log_Si_values) - (0.5)*diag((inv_proposal-inv_Si)%*%Psi + (proposal-Si)%*% Phi)) + sum(eigen_log_proposal_values - eigen_log_Si_values) + (fun(eigen_proposal_values) - fun(eigen_log_proposal_values)) - (fun(eigen_Si_values) - fun(eigen_log_Si_values))
    log_uniform <- log(runif(1, min=0, max=1))
    if (log_uniform <= log_ratio) {
      Si <- proposal
      eigen_Si_vectors <- eigen_log_proposal_vectors
      eigen_Si_values <- eigen_proposal_values
      inv_Si <- inv_proposal
      eigen_log_Si_values <- eigen_log_proposal_values
    }
    Si_rn[,,item] <- Si

    # print
    if(print & item%%1000==0){ print(item) }
  }

  # output
  return( Si_rn[,,-(1:burn), drop=FALSE] )
}



la <- 3
Psi <- diag(2)#rWishart(1L, df = 6, Sigma = diag(2L))[, , 1L]
Ga <- diag(2)
sims <- rMGIG_HR(la, Psi, Ga, mc = 20000L, burn = 1000L)


dets <- numeric(dim(sims)[3L])
for(i in seq_along(dets)) {
  dets[i] <- det(sims[, , i])
}
plot(dets)
mean(dets)

Bessel2(Psi %*% Ga / 4, la + 1 - 3/2) / Bessel2(Psi %*% Ga / 4, la - 3/2) * det(Ga / 2)
Bessel2(Psi %*% Ga / 4, la + 1) / Bessel2(Psi %*% Ga / 4, la) * det(Ga / 2)
