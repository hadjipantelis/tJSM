
#========== Function to Obtain Lamb Given Other Finite Dimensional Parameters for Model I ==========#
#=============== Transformation model is fitted for the survival part ===============#

Lamb1 <- function (para, lamb.init, tol, iter) {
  
  para.list <- Vec2List(para, ncx, ncz, ncw)
  beta <- para.list$beta
  Ysigma2 <- (para.list$Ysigma) ^ 2
  BSigma <- para.list$BSigma
  phi <- para.list$phi
  alpha <- para.list$alpha
  
  #  VY <- lapply(1:n, function(i) as.matrix(Z.st[[i]] %*% BSigma %*% t(Z.st[[i]]) + Ysigma2 * diag(1, ni[i])))
  #  VB <- lapply(1:n, function(i) BSigma - BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma) 
  #  muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta)))
  #  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b)))   
 
  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma, Ysigma2))  
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = Z.st[[i]], M_i1 = BSigma, M_i3 = VY[[i]]))
  muB <- lapply(1:n, function(i) calc_muB( BSigma, M_i3=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta, M_i1=VY[[i]], M_i2=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 

  # each element is ncz*GQ matrix #
  bi <- do.call(rbind, bi.st) # (n*ncz)*GQ matrix #
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <- do.call(rbind, lapply((1:n)[nk != 0], function(i) Ztime2.st[[i]] %*% bi.st[[i]])) # M*GQ matrix #
  
  eta.h <- as.vector(Wtime %*% phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix #

  # eta.s <- as.vector(Wtime2 %*% phi + alpha * Xtime2 %*% beta) + alpha * Ztime2.b # M*GQ matrix # 
  calc_y_a( Ztime2.b,alpha); # Ztime2.b gets altered
  
#eta.s <- as.numeric(Wtime2 %*% phi + alpha * Xtime2 %*% beta) + Ztime2.b  
  
  # exp.es <- exp(eta.s) # M*GQ matrix #
 # n_ <- ncol(eta.s)
 # m_ <- nrow(eta.s)
 # exp.es <- matrix(calc_expM(eta.s),m_,n_)
  



exp.es<- as.numeric(Wtime2 %*% phi + alpha * Xtime2 %*% beta) + Ztime2.b  
  calc_expM2(exp.es)



  lamb.old <- lamb.init
  err <- 1
  
  for (step in 1:iter) {
    Lamb.old <- cumsum(lamb.old)
    
    log.lamb <- log(lamb.old[Index0])
    log.lamb[is.na(log.lamb)] <- 0
    log.density1 <- log.lamb + eta.h # n*GQ matrix #
    const <- matrix(0, n, GQ) # n*GQ matrix #
    
    # const[nk != 0, ] <- rowsum(lamb.old[Index1] * exp.es, Index)
    const[nk != 0, ] <- calc_mult0_rowsum((Index), lamb.old[Index1], exp.es)
    log.density2 <- - log(1 + rho * const) # n*GQ matrix # 

    # log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix # 
    log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix # 
    
    f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
    deno <- as.vector(f.surv %*% wGQ) # vector of length n #
    Integral <- f.surv / deno # n*GQ matrix #
    CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
    
    # tempLamb <- (CondExp[Index, ] * exp.es * Integral[Index, ]) %*% wGQ # vector of length M # 
    # postLamb <- as.vector(tapply(tempLamb, Index1, sum)) # vector of length n_u #
    # WE DO IN PLACE MULTIPLICATION / variable 'exp.es' is hold the results
    tempLamb0 <- exp.es; tempLamb0[1] = tempLamb0[1] +0 # "touch the variable"
    calc_M1_M2_M3_Hadamard(tempLamb0, CondExp ,  Integral, as.integer(Index-1))
    tempLamb <- calc_M_y(y_i =wGQ, M_i=tempLamb0)
    postLamb <- calc_tapply_vect_sum( tempLamb, as.integer(Index1-1)); ## Check this!
    lamb.new <- Index2 / postLamb
    
    Lamb.new <- cumsum(lamb.new)
    err <- max(abs((Lamb.new - Lamb.old) / Lamb.old))
    if (step > 3 & err < tol) break
    lamb.old <- lamb.new
  }
  converge <- as.numeric(step <= iter & err < tol)
  return(list(lamb = lamb.new, converge = converge))
}
