
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model II ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterTM2 <- function (theta.old) { # Use apply instead of matrix calculation #
  
  # Get Old Estimates #
  beta.old <- theta.old$beta
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  BSigma.old <- theta.old$BSigma
  phi.old <- theta.old$phi
  alpha.old <- theta.old$alpha
  lamb.old <- theta.old$lamb
  
  # VY <- lapply(1:n, function(i) as.matrix(Z.st[[i]] %*% BSigma.old %*% t(Z.st[[i]]) + Ysigma2.old * diag(1, ni[i])))
  # VB <- lapply(1:n, function(i) BSigma.old - BSigma.old %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma.old)
  # muB <- lapply(1:n, function(i) as.vector(BSigma.old %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta.old))) 
  # bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b))) 

  #VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma.old, Ysigma2.old))  
  #VB <-  lapply(1:n, function(i) calc_VB(y_i = Z.st[[i]], a_i = BSigma.old, VY[[i]])) 
  #muB <- lapply(1:n, function(i) calc_muB( BSigma.old,  y_i0=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta.old, M_i1=VY[[i]], M_i2=X.st[[i]])) 
  #bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) )



  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma.old, Ysigma2.old)) 
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = Z.st[[i]], M_i1 = BSigma.old, M_i3 = VY[[i]]))
  muB <- lapply(1:n, function(i) calc_muB( BSigma.old,  y_i0=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta.old, M_i1=VY[[i]], M_i2=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) )
 
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <- do.call(rbind, lapply((1:n)[nk != 0], function(i) Ztime2.st[[i]] %*% bi.st[[i]])) # M*GQ matrix #
  
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi.old) + alpha.old * Ztime.b # n*GQ matrix #
  eta.s <- as.vector(Wtime2 %*% phi.old) + alpha.old * Ztime2.b # M*GQ matrix #
  # exp.es <- exp(eta.s) # M*GQ matrix #
  n_ <- ncol(eta.s)
  m_ <- nrow(eta.s)
  exp.es <- matrix(calc_expM(eta.s),m_,n_) # This faster for rectangular matrices
  const <- matrix(0, n, GQ) # n*GQ matrix #
  temp0a <- exp.es * lamb.old[Index1]
  #const[nk != 0, ] <- rowsum(lamb.old[Index1] * exp.es, Index)
  const[nk != 0, ] <- calc_rowsum( (Index), temp0a)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  # log.survival <- if (rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix # 
  log.survival <- if (rho > 0) log.density2 / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix f(bi|Oi) #
  
  # f.long <- sapply(1:n, function(i) dmvnorm(Y.st[[i]], as.vector(X.st[[i]] %*% beta.old), VY[[i]]))
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta.old), VY[[i]]))
  lgLik <- sum(log(f.long * deno / (pi ^ (ncz / 2))))
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- if(ncz > 1) t(sapply(1:n, function(i) post.bi[i, ((i - 1) * ncz + 1) : (i * ncz)])) else 
             matrix(diag(post.bi), nrow = n) # n*ncz matrix Ehat(bi) #
  
  #========== Update BSigma ==========#
  if (ncz>1) {
    tempB <- do.call(rbind, lapply(1:n, function(i) apply(t(bi.st[[i]]), 1, function(x) x %o% x)))
    # (n*ncz^2)*GQ matrix #      
  } else {
    tempB <- bi ^ 2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1) : (i * ncz2)])) else 
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix Ehat(bibi^T) #
  BSigma.new <- if (ncz > 1) matrix(colMeans(post.bi2), ncz, ncz) else mean(post.bi2) # ncz*ncz matrix #
  
  #========== Update beta: the linear regresion coefficents of regression Yi-E(Zi*bi) on X_i ==========#
  tempX <- Y - rowSums(Z * post.bi[ID, ]) # vector of length N #
beta.new <- as.vector(solve(t(X) %*% X) %*% (t(X) %*% tempX)) # vector of length ncx #
   #  beta.new <- as.vector(qr.solve(X,tempX));
  #========== Update Ysigma ==========#
  Ymu.new <- as.vector(X %*% beta.new) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) 
  # N*GQ matrix #
  post.resid <- ((Y - Ymu.new) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and alpha ==========#
  CondExp2 <- CondExp[nk != 0, ]
  # temp1 <- CondExp2 * rowsum(Ztime2.b * exp.es * lamb.old[Index1], Index) # n*GQ matrix #
  # temp2 <- CondExp2 * rowsum(Ztime2.b ^ 2 * exp.es * lamb.old[Index1], Index) # n*GQ matrix #
  # temp3 <- lapply(1:ncw, function(i) CondExp2 * rowsum(Wtime2[, i] * exp.es * lamb.old[Index1], Index))
  # temp4 <- lapply(1:(ncw ^ 2), function(i) CondExp2 * rowsum(Wtime22[, i] * exp.es * lamb.old[Index1], Index)) 
  # temp5 <- lapply(1:ncw, function(i) CondExp2 * rowsum(Ztime2.b * Wtime2[, i] * exp.es * lamb.old[Index1], Index))

  
  temp0b <- Ztime2.b * temp0a
  temp1 <- CondExp2 * rowsum(temp0b, Index) # n*GQ matrix #
  temp2 <- CondExp2 * rowsum(Ztime2.b * temp0b, Index) # n*GQ matrix #
  temp3 <- lapply(1:ncw, function(i) CondExp2 * rowsum(Wtime2[, i] * temp0a, Index))
  temp4 <- lapply(1:(ncw ^ 2), function(i) CondExp2 * rowsum(Wtime22[, i] * temp0a, Index)) 
  temp5 <- lapply(1:ncw, function(i) CondExp2 * rowsum(  Wtime2[, i] * temp0b, Index))

  Integral2 <- Integral[nk != 0,]
  post1 <- sum((temp1 * Integral2) %*% wGQ)
  post2 <- sum((temp2 * Integral2) %*% wGQ)
  post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
  post4 <- unlist(lapply(temp4, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw^2 #
  post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw # 
  
  phiScore <- colSums(d * Wtime) - post3 # vector of length ncw #
  alphaScore <- sum(d * rowSums(Ztime * post.bi)) - post1
  pa.score <- c(phiScore, alphaScore)
  pa.info <- matrix(0, (ncw + 1), (ncw + 1)) # (ncw+1)*(ncw+1) matrix #
  pa.info[1:ncw, 1:ncw] <- -post4
  pa.info[(ncw + 1), 1:ncw] <- -post5
  pa.info[1:ncw, (ncw + 1)] <- -post5
  pa.info[(ncw + 1), (ncw + 1)] <- -post2
  
  #=============== Update phi and alpha ===============#
  pa.old <- c(phi.old, alpha.old) # vector of length (ncw+1) #
  paSVD <- svd(pa.info)
  pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
  pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
  phi.new <- pa.new[1:ncw]
  alpha.new <- pa.new[ncw + 1]
  
  #========== Calculate the new lambda with new parameters ==========#
  eta.sn <- as.vector(Wtime2 %*% phi.new) + alpha.new * Ztime2.b # M*GQ matrix #  ## CHECK IF YOU CAN INPLACE THIS
  # tempLamb <- (CondExp[Index, ] * exp.eta.sn * Integral[Index, ]) %*% wGQ # vector of length M #
  # postLamb <- as.vector(tapply(tempLamb, Index1, sum)) # vector of length n_u #

  n_ <- ncol(eta.sn)
  m_ <- nrow(eta.sn)
  exp.eta.es <- matrix(calc_expM(eta.sn),m_,n_) # This faster for rectangular matrices
  calc_M1_M2_M3_Hadamard(exp.eta.es, CondExp ,  Integral, as.integer(Index-1))
  tempLamb <- calc_M_y(y_i =wGQ, M_i=exp.eta.es)
  postLamb <- calc_tapply_vect_sum( tempLamb, as.integer(Index1-1)); ## Check this!

  lamb.new <- Index2 / postLamb
  
  result <- list(beta = beta.new, Ysigma = sqrt(Ysigma2.new), BSigma = BSigma.new, phi = phi.new, 
                 alpha = alpha.new, lamb = lamb.new, lgLik = lgLik, est.bi = post.bi)
  return(result)
}
