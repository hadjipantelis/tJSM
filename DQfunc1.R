
#=============== The DQ Function for Model I ===============#
#=============== Transformation model is fitted for the survival part ===============#

DQfunc1 <- function (ptheta, theta) { # ptheta means "theta prime"

  pbeta <- ptheta$beta
  beta <- theta$beta
  pYsigma2 <- (ptheta$Ysigma) ^ 2
  Ysigma2 <- (theta$Ysigma) ^ 2
  pBSigma <- ptheta$BSigma
  BSigma <- theta$BSigma
  pphi <- ptheta$phi
  phi <- theta$phi
  palpha <- ptheta$alpha
  alpha <- theta$alpha
  plamb <- ptheta$lamb
  lamb <- theta$lamb
  

  # VY <- lapply(1:n, function(i) as.matrix(Z.st[[i]] %*% BSigma %*% t(Z.st[[i]]) + Ysigma2 * diag(1, ni[i])))
  # VB <- lapply(1:n, function(i) BSigma - BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma)
  # muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta)))
  # bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b)))
  
  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma, Ysigma2)) 
  VB <-  lapply(1:n, function(i) calc_VB(y_i = Z.st[[i]], a_i = BSigma, VY[[i]]))
  muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]],Y.st[[i]] - X.st[[i]] %*% beta)))
  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * backsolve(chol(solve(VB[[i]])), t(b))))

  bi <- do.call(rbind, bi.st) # (n*ncz)*GQ matrix #
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <- do.call(rbind, lapply((1:n)[nk != 0], function(i) Ztime2.st[[i]] %*% bi.st[[i]])) # M*GQ matrix #
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix #
  eta.s <- as.vector(Wtime2 %*% phi + alpha * Xtime2 %*% beta) + alpha * Ztime2.b # M*GQ matrix #
  # exp.es <- exp(eta.s) # M*GQ matrix #
  n_ <- ncol(eta.s)
  m_ <- nrow(eta.s)
  exp.es <- matrix(example_expM(A1),m_,n_)
  const <- matrix(0, n, GQ) # n*GQ matrix #
  const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index) 
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix #
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  len <- ncx + ncw + p + 2
  Q <- rep(0, len)
  
  post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- if(ncz > 1) t(sapply(1:n, function(i) post.bi[i, ((i - 1) * ncz + 1):(i * ncz)])) else 
             matrix(diag(post.bi), nrow = n) # n*ncz matrix #
  
  if (ncz>1) {
    tempB <- do.call(rbind, lapply(1:n, function(i) apply(t(bi.st[[i]]), 1, function(x) x %o% x)))
    # (n*ncz^2)*GQ matrix #     
  }else{
    tempB <- bi ^ 2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1):(i * ncz2)])) else 
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix #
  
  pYmu <- as.vector(X %*% pbeta) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) 
  post.resid <- ((Y - pYmu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Q[ncx + ncw + 2] <- - N / sqrt(pYsigma2) + sum(post.resid) / (pYsigma2 ^ (3 / 2))
  
  #tempB <- - n * solve(pBSigma) / 2 + solve(pBSigma) %*% matrix(colSums(post.bi2), ncz, ncz) %*% solve(pBSigma) / 2
  pBSigmaInv = solve(pBSigma);
  tempB <- - n * pBSigmaInv / 2 + pBSigmaInv %*% matrix(colSums(post.bi2), ncz, ncz) %*% pBSigmaInv / 2
  ind <- Indexing(ncz)
  Q[(ncx + ncw + 3):len] <- as.vector(tapply(c(tempB), ind, sum))
  
  eta.sp <- as.vector(Wtime2 %*% pphi + palpha * Xtime2 %*% pbeta) + palpha * Ztime2.b # M*GQ matrix #
  exp.esp <- exp(eta.sp) # M*GQ matrix #
  XZb2 <- as.vector(Xtime2 %*% pbeta) + Ztime2.b # M*GQ matrix #
  temp1 <- as.vector((CondExp[Index, ] * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  temp2 <- as.vector((CondExp[Index, ] * XZb2 * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  temp3 <- Wtime2 * temp1 # M*ncw matrix # 

  post1 <- as.vector(tapply(temp1, Index1, sum)) # vector of length n_u #
  post2 <- as.vector(tapply(temp2, Index1, sum)) # vector of length n_u #
  post3 <- as.matrix(apply(temp3, 2, function(x) tapply(x, Index1, sum))) # n_u*ncw matrix #
  
  Q[(ncx + 1):(ncx + ncw)] <- colSums(d * Wtime) - colSums(Index2 * post3 / post1) # vector of length ncw #
  Q[ncx + ncw + 1] <- sum(d * Xtime %*% pbeta) + sum(d * rowSums(Ztime * post.bi)) - sum(Index2 * post2 / post1)
  
  temp4 <- Xtime2 * temp1 # M*ncx matrix #
  post4 <- palpha * as.matrix(apply(temp4, 2, function(x) tapply(x, Index1, sum))) # n_u*ncx matrix #
  pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
  Q[1:ncx] <- colSums(X * pResid) / pYsigma2 + palpha * colSums(d * Xtime) - colSums(Index2 * post4 / post1) 
  
  return(Q)
}
