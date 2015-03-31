
#=============== The DQ Function for Model II ===============#
#=============== Transformation model is fitted for the survival part ===============#

DQfunc2 <- function (ptheta, theta) { # ptheta means "theta prime"
  ## This might have significant (e-14) numerical difference with the original code. Something odd is happening with VB/muB.
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
  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma, Ysigma2))
  # VB <- lapply(1:n, function(i) BSigma - BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma) 
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = Z.st[[i]], M_i1 = BSigma, M_i3 = VY[[i]]))
  # muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta)))  
  muB <- lapply(1:n, function(i) calc_muB( BSigma, M_i3=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta, M_i1=VY[[i]], M_i2=X.st[[i]]))
  # bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b)))   
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # (n*ncz)*GQ matrix #
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix # 

  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk !=      0] - 1)# M*GQ matrix # 
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi) + alpha*Ztime.b # n*GQ matrix #

  # eta.s <- as.vector(Wtime2 %*% phi) + alpha * Ztime2.b # M*GQ matrix #
  # exp.es <- exp(eta.s) # M*GQ matrix #

  calc_y_a( Ztime2.b,alpha); # Ztime2.b gets altered because we do in-place multiplication - Ztime2.b *= alpha
  exp.es <- as.numeric(Wtime2 %*% phi) + Ztime2.b  
  calc_expM2(exp.es)

  const <- matrix(0, n, GQ) # n*GQ matrix #

  #   const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index)  
  const[nk != 0, ] <- calc_mult0_rowsum((Index), lamb[Index1], exp.es)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 

  #  log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix #
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix #
  
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
   # tempB <- do.call(rbind, lapply(1:n, function(i) apply(t(bi.st[[i]]), 1, function(x) x %o% x)))
   # tempB <- do.call(rbind, lapply(1:n, function(i) apply(t(bi.st[[i]]), 1, function(x) tcrossprod(x)))) 
    tempB <-  fast_rbind_lapply( bi.st )
    # (n*ncz^2)*GQ matrix #     
  } else {
    tempB <- bi ^ 2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1):(i * ncz2)])) else 
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix #
  
  pYmu <- as.vector(X %*% pbeta) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) 
  post.resid <- ((Y - pYmu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Q[ncx + ncw + 2] <- - N / sqrt(pYsigma2) + sum(post.resid) / (pYsigma2 ^ (3 / 2))
  
  # tempB <- - n * solve(pBSigma) / 2 + solve(pBSigma) %*% matrix(colSums(post.bi2), ncz, ncz) %*% solve(pBSigma) / 2
  pBSigmaInv = solve(pBSigma);
  tempB <- - n * pBSigmaInv / 2 + pBSigmaInv %*% matrix(colSums(post.bi2), ncz, ncz) %*% pBSigmaInv / 2

  ind <- Indexing(ncz)
  Q[(ncx + ncw + 3):len] <- as.vector(tapply(c(tempB), ind, sum))
  
  # eta.sp <- as.vector(Wtime2 %*% pphi) + palpha * Ztime2.b # M*GQ matrix #
  calc_y_a( Ztime2.b,palpha/ alpha) # Ztime2.b gets altered because we do in-place multiplication - Ztime2.b *= (palpha/ alpha)
  # eta.sp <- as.numeric(Wtime2 %*% pphi) +  Ztime2.b # M*GQ matrix #
  # exp.esp <- exp(eta.sp) # M*GQ matrix #
  
  exp.esp <-as.numeric(Wtime2 %*% pphi) +  Ztime2.b # M*GQ matrix #
  calc_expM2(exp.esp)


  # temp1 <- as.vector((CondExp[Index, ] * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  # temp2 <- as.vector((CondExp[Index, ] * Ztime2.b * exp.esp * Integral[Index, ]) %*% wGQ) # vector of length M #
  
  temp0 <- exp.esp; temp0[1] = temp0[1] +0 # "touch the variable"
  calc_M1_M2_M3_Hadamard(temp0, CondExp ,  Integral, as.integer(Index-1))
  temp1 <- calc_M_y(y_i = wGQ, M_i = temp0) 
  calc_y_a( Ztime2.b,1/palpha) # Ztime2.b gets altered because we do in-place multiplication - Ztime2.b *= (1/ alpha)
  calc_M1_M2_Hadamard(temp0,Ztime2.b)
  temp2 <- calc_M_y(y_i =wGQ, M_i= temp0)

  temp3 <- Wtime2 * temp1 # M*ncw matrix #
   
  # post1 <- as.vector(tapply(temp1, Index1, sum)) # vector of length n_u #
  # post2 <- as.vector(tapply(temp2, Index1, sum)) # vector of length n_u #.
  # post3 <- as.matrix(apply(temp3, 2, function(x) tapply(x, Index1, sum))) # n_u*ncw matrix #
  post1 <- calc_tapply_vect_sum( temp1, as.integer(Index1-1));
  post2 <- calc_tapply_vect_sum( temp2, as.integer(Index1-1));  
  post3 <- as.matrix(apply(temp3, 2, function(x) calc_tapply_vect_sum( x, as.integer(Index1-1))))

  Q[(ncx + 1):(ncx + ncw)] <- colSums(d * Wtime) - colSums(Index2 * post3 / post1) # vector of length ncw #
  Q[ncx + ncw + 1] <- sum(d * rowSums(Ztime * post.bi)) - sum(Index2 * post2 / post1)
  
  pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
  Q[1:ncx] <- colSums(X * pResid) / pYsigma2 # vector of length ncx # 

return(Q)
}
