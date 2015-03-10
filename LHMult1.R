
#=============== Function to Calculate the Likelihood Value for Model I ===============#
#================ Multiplicative Joint Modeling ===============#

LHMult1 <- function (theta) {
  
  gamma <- theta$gamma
  phi <- theta$phi
  alpha <- theta$alpha
  Ysigma2 <- (theta$Ysigma) ^ 2
  Bsigma2 <- (theta$Bsigma) ^ 2
  lamb <- theta$lamb
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) as.matrix(Bsigma2 * BTg[[i]] %*% t(BTg[[i]]) + Ysigma2 * diag(1, ni[i])))
  #VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * t(BTg[[i]]) %*% solve(VY[[i]]) %*% BTg[[i]])) 
  VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) *  sum( forwardsolve(t(chol(VY[[i]])), BTg[[i]])^2)))
  #muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - BTg[[i]])))
  muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]],as.vector(Y.st[[i]] - BTg[[i]]))))

  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2 * VB[[i]]) * t(b)))
  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  Btime.b <- as.vector(Btime %*% gamma) * bi # n*nknot matrix #
  Btime2.b <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Ztime %*% phi) + alpha * Btime.b # n*nknot matrix #
  eta.s <- as.vector(Ztime2 %*% phi) + alpha * Btime2.b # M*nknot matrix #
  exp.es <- exp(eta.s) # M*nknot matrix #
  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index) # n*nknot matrix # 
  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
  log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  
  f.long <- sapply(1:n, function(i) dmvnorm(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]])) # vector of length n #
  return(sum(log(f.long * deno / sqrt(pi))))
}
