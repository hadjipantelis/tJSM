
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model II with NMRE ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterMult2 <- function (theta.old) { # Use apply instead of matrix calculation #
  
  # Get Old Estimates #
  gamma.old <- theta.old$gamma
  phi.old <- theta.old$phi
  alpha.old <- theta.old$alpha
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  Bsigma2.old <- (theta.old$Bsigma) ^ 2
  lamb.old <- theta.old$lamb
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma.old))
  VY <- lapply(1:n, function(i) calc_VY( BTg[[i]], Bsigma2.old, Ysigma2.old) )
  VB <-  lapply(1:n, function(i) calc_VB(M1 = Bsigma2.old, M2 = BTg[[i]], M3 = VY[[i]]))
  muB <- lapply(1:n, function(i) calc_muBMult(  Bsigma2.old,VY[[i]],BTg[[i]],Y.st[[i]] )+1 )
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Ztime %*% phi.old) + alpha.old * bi # n*nknot matrix #
  eta.s <- as.vector(Ztime2 %*% phi.old) + alpha.old * bi[Index, ] # M*nknot matrix #
  exp.es <- exp(eta.s) # M*nknot matrix #
  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- rowsum(lamb.old[Index1] * exp.es, Index) # n*nknot matrix # 
  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
  log.survival <- if (rho > 0) - log(1 + rho * const) / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*nknot matrix #
  
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]])) # vector of length n #
  lgLik <- sum(log(f.long * deno / sqrt(pi)))
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*nknot matrix #
  
  #========== Update Bsigma ==========#
  Bsigma2.new <- mean((Integral * (bi - 1) ^ 2) %*% wGQ) 
  
  #========== Update gamma: the linear regresion coefficents of regression Yi on E(bi)*B ==========#
  post.bi <- as.vector((Integral * bi) %*% wGQ) # vector of length n #
  post.bi2 <- as.vector((Integral * bi ^ 2) %*% wGQ) # vector of length n #
  gamma.new <- as.vector(solve(matrix(colSums(post.bi2[ID] * B2), ncb)) %*% colSums(post.bi[ID] * B * Y)) 
  # vector of length ncb #
  
  #========== Update Ysigma ==========#
  post.resid <- ((Y - as.vector(B %*% gamma.new) * bi[ID, ]) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and alpha ==========#
  CondExp2 <- CondExp[nk != 0, ]
  temp1 <- lapply(1:ncz, function(i) CondExp2 * rowsum(Ztime2[, i] * exp.es * lamb.old[Index1], Index)) 
  # n*nknot matrices #
  temp2 <- CondExp2 * rowsum(bi[Index, ] * exp.es * lamb.old[Index1], Index) # n*nknot matrix #
  temp3 <- lapply(1:(ncz ^ 2), function(i) CondExp2 * rowsum(Ztime22[, i] * exp.es * lamb.old[Index1], Index)) 
  # n*nknot matrices #
  temp4 <- CondExp2 * rowsum(bi[Index, ] ^ 2 * exp.es * lamb.old[Index1], Index) # n*nknot matrix #
  temp5 <- lapply(1:ncz, function(i) CondExp2 * rowsum(bi[Index, ] * Ztime2[, i] * 
                                     exp.es * lamb.old[Index1], Index)) # n*nknot matrices #
  Integral2 <- Integral[nk != 0, ]
  post1 <- unlist(lapply(temp1, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz #
  post2 <- sum((temp2 * Integral2) %*% wGQ)
  post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz^2 #
  post4 <- sum((temp4 * Integral2) %*% wGQ)
  post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz #
  
  phiScore <- colSums(d * Ztime) - post1 # vector of length ncz #
  alphaScore <- sum(d * post.bi) - post2
  pa.score <- c(phiScore, alphaScore)
  pa.info <- matrix(0, (ncz + 1), (ncz + 1)) # (ncz+1)*(ncz+1) matrix #
  pa.info[1:ncz, 1:ncz] <- - post3
  pa.info[(ncz + 1), (ncz + 1)] <- - post4
  pa.info[(ncz + 1), 1:ncz] <- - post5
  pa.info[1:ncz, (ncz + 1)] <- - post5
  
  #=============== Update phi and alpha ===============#
  pa.old <- c(phi.old, alpha.old) # vector of length (ncz+1) #
  paSVD <- svd(pa.info)
  pa.info.inv <- paSVD$v %*% diag(1 / paSVD$d) %*% t(paSVD$u)
  pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncz+1) #
  phi.new <- pa.new[1 : ncz]
  alpha.new <- pa.new[ncz + 1]
  
  #========== Calculate the new lambda with new parameters ==========#
  eta.s.n <- as.vector(Ztime2 %*% phi.new) + alpha.new * bi[Index, ] # M*nknot matrix #
  tempLamb <- (CondExp[Index, ] * exp(eta.s.n) * Integral[Index, ]) %*% wGQ # vector of length M #
  postLamb <- as.vector(tapply(tempLamb, Index1, sum)) # vector of length n_u #
  lamb.new <- Index2 / postLamb
  
  result <- list(gamma = gamma.new, phi = phi.new, alpha = alpha.new, Ysigma = sqrt(Ysigma2.new), 
                 Bsigma = sqrt(Bsigma2.new), lamb = lamb.new, lgLik = lgLik, est.bi = post.bi)
  return(result)
}
