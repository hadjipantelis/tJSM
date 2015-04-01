
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model I ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterTM1 <- function (theta.old) { # Use apply instead of matrix calculation #
  
  # Get Old Estimates #
  beta.old <- theta.old$beta
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  BSigma.old <- theta.old$BSigma
  phi.old <- theta.old$phi
  alpha.old <- theta.old$alpha
  lamb.old <- theta.old$lamb
 
  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma.old, Ysigma2.old)) 
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = Z.st[[i]], M_i1 = BSigma.old, M_i3 = VY[[i]]))
  muB <- lapply(1:n, function(i) calc_muB( BSigma.old, M_i3=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta.old, M_i1=VY[[i]], M_i2=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 
  
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk != 0] - 1)# M*GQ matrix # 
    
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi.old + alpha.old * Xtime %*% beta.old) + alpha.old * Ztime.b # n*GQ matrix #
  eta.s <- as.vector(Wtime2 %*% phi.old + alpha.old * Xtime2 %*% beta.old) + alpha.old * Ztime2.b
  const <- matrix(0, n, GQ) # n*GQ matrix #
   
  calc_expM2(eta.s)
  temp0a <- eta.s * lamb.old[Index1]; 
  const[nk != 0, ] <- calc_rowsum( (Index), temp0a)
  log.density2 <- -log(1 + rho * const) # n*GQ matrix #  
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix # 

  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix #
    
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta.old), VY[[i]]))

  lgLik <- sum(log(f.long * deno / (pi ^ (ncz / 2))))
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- if(ncz > 1) {t(sapply(1:n, function(i) post.bi[i, ((i - 1) * ncz + 1) : (i * ncz)])) } else {
             matrix(diag(post.bi), nrow = n) } # n*ncz matrix Ehat(bi) #
  
  #========== Update BSigma ==========#
  if (ncz > 1) {      
     tempB <-  fast_rbind_lapply( bi.st )     # (n*ncz^2)*GQ matrix #      
  } else {
    tempB <- bi^2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1) : (i * ncz2)])) else
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix Ehat(bibi^T) #
  BSigma.new <- if (ncz > 1) matrix(colMeans(post.bi2), ncz, ncz) else mean(post.bi2) # ncz*ncz matrix #
  
  #========== Update Ysigma ==========#
  Ymu <- as.vector(X %*% beta.old) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) # N*GQ matrix #
  post.resid <- ((Y - Ymu)^2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and alpha ==========# 
  XZb2 <- as.vector(Xtime2 %*% beta.old) + Ztime2.b # M*GQ matrix #
  CondExp2 <- CondExp[nk != 0, ]
  temp0b <- XZb2 * temp0a; 

  temp1 <- CondExp2 * calc_rowsum( (Index), temp0b)
  temp2 <- calc_mult_rowsum2(y_i = Index, M_i2 = CondExp2, M_i1 = temp0b,      XZb2)
  temp3 <- lapply(1:(ncw), function(i) calc_mult_rowsum(y_i = Index, y_i2 = Wtime2[, i], M_i2 = CondExp2, temp0a))
  temp4 <- lapply(1:(ncw^2), function(i) calc_mult_rowsum(y_i = Index, y_i2 = Wtime22[, i], M_i2 = CondExp2, temp0a)) 
  temp5 <- lapply(1:(ncw), function(i) calc_mult_rowsum(y_i = Index, y_i2 = Wtime2[, i], M_i2 = CondExp2, XZb2 *temp0a)) 

  Integral2 <- Integral[nk != 0, ]
  post1 <- sum((temp1 * Integral2) %*% wGQ)
  post2 <- sum((temp2 * Integral2) %*% wGQ)
  post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
  post4 <- unlist(lapply(temp4, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw^2 #
  post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw # 
   
  phiScore <- colSums(d * Wtime) - post3 # vector of length ncw #
  alphaScore <- sum(d * Xtime %*% beta.old) + sum(d * rowSums(Ztime * post.bi)) - post1
  pa.score <- c(phiScore, alphaScore)
  pa.info <- matrix(0, (ncw + 1), (ncw + 1)) # (ncw+1)*(ncw+1) matrix #
  pa.info[1:ncw, 1:ncw] <- -post4
  pa.info[(ncw + 1), 1:ncw] <- -post5
  pa.info[1:ncw, (ncw + 1)] <- -post5
  pa.info[(ncw + 1), (ncw + 1)] <- -post2
  
  #=============== Update phi and alpha ===============#
  pa.old <- c(phi.old, alpha.old) # vector of length (ncw+1) #
  paSVD <- svd(pa.info)
  pa.info.inv <- paSVD$v %*% diag(1 / paSVD$d) %*% t(paSVD$u)
  pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
  phi.new <- pa.new[1:ncw]
  alpha.new <- pa.new[ncw + 1]
  
  #========== calculate the score of beta ==========#
  newZtime2.b = alpha.new * Ztime2.b
  eta.s.n1 <- as.vector(Wtime2 %*% phi.new + alpha.new * Xtime2 %*% beta.old) + newZtime2.b # M*GQ matrix # 
  calc_expM2(eta.s.n1)
  temp0c <- alpha.new * eta.s.n1 * lamb.old[Index1]; 
  temp6 <- lapply(1:(ncx), function(i) calc_mult_rowsum(y_i = Index, y_i2 = Xtime2[, i], M_i2 = CondExp2, temp0c))
  temp0d <- alpha.new*temp0c 
  temp7 <- lapply(1:(ncx^2), function(i) calc_mult_rowsum(y_i = Index, y_i2 = Xtime22[, i], M_i2 = CondExp2, temp0d))


  post6 <- unlist(lapply(temp6, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncx #
  post7 <- unlist(lapply(temp7, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncx^2 #
  post8 <- as.vector(Y - X %*% beta.old - rowSums(Z * post.bi[ID, ])) * X # N*ncx matrix #
  
  betaScore <- colSums(post8) / Ysigma2.new + alpha.new * colSums(d * Xtime) - post6
  betaInfo <- - X2.sum / Ysigma2.new - post7
  
  #========== Update beta ==========#
  betaSVD <- svd(betaInfo)
  betaInfo.inv <- betaSVD$v %*% diag(1 / betaSVD$d) %*% t(betaSVD$u)
  beta.new <- as.vector(beta.old - betaInfo.inv %*% betaScore) # vector of length ncx #
  
  #========== Calculate the new lambda with new parameters ==========#
  eta.s.n2 <- as.vector(Wtime2 %*% phi.new + alpha.new * Xtime2 %*% beta.new) + newZtime2.b  # M*GQ matrix # 
  calc_expM2(eta.s.n2) # This is exponentiated 
  calc_M1_M2_M3_Hadamard(eta.s.n2, CondExp ,  Integral, as.integer(Index-1))
  tempLamb <- calc_M_y(y_i =wGQ, M_i=eta.s.n2)
  postLamb <- calc_tapply_vect_sum( tempLamb, as.integer(Index1-1)); ## Check this!

  lamb.new <- Index2 / postLamb
  
  result <- list(beta = beta.new, Ysigma = sqrt(Ysigma2.new), BSigma = BSigma.new, phi = phi.new, 
                 alpha = alpha.new, lamb = lamb.new, lgLik = lgLik, est.bi = post.bi)
  return(result)
}
