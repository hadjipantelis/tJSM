
#========== Function to Obtain Lamb Given Other Finite Dimensional Parameters ==========#
#=============== Model II for Multiplicative Joint Modeling ===============#

LambMult2 <- function (para, lamb.init, tol, iter) {
  
  para.list <- Vec2ListMult(para, ncz, ncb)
  gamma <- para.list$gamma
  phi <- para.list$phi
  alpha <- para.list$alpha
  Ysigma2 <- (para.list$Ysigma) ^ 2
  Bsigma2 <- (para.list$Bsigma) ^ 2
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) as.matrix(Bsigma2 * BTg[[i]] %*% t(BTg[[i]]) + Ysigma2 * diag(1, ni[i])))
  VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * t(BTg[[i]]) %*% solve(VY[[i]]) %*% BTg[[i]]))
 #  VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * sum( forwardsolve(t(chol(VY[[i])), BTg[[i]])^2)))
  muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - BTg[[i]])))
  # muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]],as.vector(Y.st[[i]] - BTg[[i]]))))

  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2 * VB[[i]]) * t(b)))
  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  
  eta.h <- as.vector(Ztime %*% phi) + alpha * bi # n*nknot matrix #
  eta.s <- as.vector(Ztime2 %*% phi) + alpha * bi[Index, ] # M*nknot matrix #
  exp.es <- exp(eta.s) # M*nknot matrix #
  
  lamb.old <- lamb.init
  err <- 1
  
  for (step in 1:iter) {
    Lamb.old <- cumsum(lamb.old)
    
    log.lamb <- log(lamb.old[Index0])
    log.lamb[is.na(log.lamb)] <- 0
    log.density1 <- log.lamb + eta.h # n*nknot matrix #
    const <- matrix(0, n, nknot) # n*nknot matrix #
    const[nk != 0, ] <- rowsum(lamb.old[Index1] * exp.es, Index) # n*nknot matrix # 
    log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
    log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*nknot matrix # 
    
    f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
    deno <- as.vector(f.surv %*% wGQ) # vector of length n #
    Integral <- f.surv / deno # n*nknot matrix #
    CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*nknot matrix #
    
    tempLamb <- (CondExp[Index, ] * exp.es * Integral[Index, ]) %*% wGQ # vector of length M #
    postLamb <- as.vector(tapply(tempLamb, Index1, sum)) # vector of length n_u #
    lamb.new <- Index2 / postLamb
    
    Lamb.new <- cumsum(lamb.new)
    err <- max(abs((Lamb.new - Lamb.old) / Lamb.old))
    if (step > 3 & err < tol) break
    lamb.old <- lamb.new
  }
  converge <- as.numeric(step <= iter & err < tol)
  return(list(lamb = lamb.new, converge = converge))
}
