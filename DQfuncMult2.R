
#=============== The DQ Function for Model II of Multiplicative Joint Model ===============#

DQfuncMult2 <- function (ptheta, theta) { # ptheta means "theta prime"
  
  pgamma <- ptheta$gamma
  pphi <- ptheta$phi
  palpha <- ptheta$alpha
  pYsigma2  <- (ptheta$Ysigma) ^ 2
  pBsigma2 <- (ptheta$Bsigma) ^ 2
  plamb <- ptheta$lamb
  gamma <- theta$gamma
  phi <- theta$phi
  alpha <- theta$alpha
  Ysigma2  <- (theta$Ysigma) ^ 2
  Bsigma2 <- (theta$Bsigma) ^ 2
  lamb <- theta$lamb
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) as.matrix(Bsigma2 * BTg[[i]] %*% t(BTg[[i]]) + Ysigma2 * diag(1, ni[i])))
  #VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * t(BTg[[i]]) %*% solve(VY[[i]]) %*% BTg[[i]]))
  VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * sum( forwardsolve(t(chol(VY[[i]])), BTg[[i]])^2))) 
  #muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - BTg[[i]])))
  muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]], as.vector(Y.st[[i]] - BTg[[i]]))))

  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2 * VB[[i]]) * t(b)))
  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Ztime %*% phi) + alpha * bi # n*nknot matrix #
  eta.s <- as.vector(Ztime2 %*% phi) + alpha * bi[Index, ] # M*nknot matrix #
  exp.es <- exp(eta.s) # M*nknot matrix #
  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index) # n*nknot matrix # 
  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
  log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*nknot matrix #
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*nknot matrix #
  
  len <- ncz + ncb + 3
  Q <- rep(0, len)
  
  Q[len] <- - n / sqrt(pBsigma2) + sum((Integral * (bi - 1) ^ 2) %*% wGQ) / (pBsigma2 ^ (3 / 2))
  Q[len - 1] <- - N / sqrt(pYsigma2) + sum(((Y - as.vector(B %*% pgamma) * bi[ID, ]) ^ 2 * 
                                            Integral[ID, ]) %*% wGQ) / (pYsigma2 ^ (3 / 2))
  
  eta.sp <- as.vector(Ztime2 %*% pphi) + palpha * bi[Index, ] # M*nknot matrix #
  exp.esp <- exp(eta.sp) # M*nknot matrix #
  
  temp1 <- as.vector((CondExp[Index, ] * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  temp2 <- as.vector((CondExp[Index, ] * bi[Index, ] * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  temp3 <- Ztime2 * temp1 # M*ncz matrix #
  temp4 <- lapply(1:ncb, function(i) (Y - as.vector(B %*% pgamma) * bi[ID, ]) * B[, i] * bi[ID, ]) # N*nknot matrices #
  
  post1 <- as.vector(tapply(temp1, Index1, sum)) # vector of length n_u #
  post2 <- as.vector(tapply(temp2, Index1, sum)) # vector of length n_u #
  post3 <- as.matrix(apply(temp3, 2, function(x) tapply(x, Index1, sum))) # n_u*ncz matrix #
  post4 <- unlist(lapply(temp4, function(x) sum((x * Integral[ID, ]) %*% wGQ))) # vector of length ncb #
  post.bi <- as.vector((Integral * bi) %*% wGQ) # vector of length n #
  
  Q[1 : ncb] <- post4 / pYsigma2
  Q[(ncb + 1) : (ncb + ncz)] <- colSums(d * Ztime) - colSums(Index2 * post3 / post1) # vector of length ncz #
  Q[ncb + ncz + 1] <- sum(d * post.bi) - sum(Index2 * post2 / post1)
  
  return(Q)
}
