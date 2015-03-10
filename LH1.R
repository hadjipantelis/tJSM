
#=============== Function to Calculate the Likelihood Value for Model I ===============#

LH1 <- function (theta) {
  
  beta <- theta$beta
  Ysigma2 <- (theta$Ysigma) ^ 2
  BSigma <- theta$BSigma
  phi <- theta$phi
  alpha <- theta$alpha
  lamb <- theta$lamb
  
  VY <- lapply(1:n, function(i) as.matrix(Z.st[[i]] %*% BSigma %*% t(Z.st[[i]]) + Ysigma2 * diag(1,ni[i])))
  # VB <- lapply(1:n, function(i) BSigma - BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma)
  VB <- lapply(1:n, function(i) BSigma - sum( forwardsolve(t(chol(VY[[i]])), Z.st[[i]])^2)*BSigma*BSigma)
  # muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta))) 
  muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]], as.vector(Y.st[[i]] - X.st[[i]] %*% beta))))
  # bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b)))
  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * backsolve(chol(solve(VB[[i]])), t(b))))

  # each element is ncz*GQ matrix #
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <- do.call(rbind, lapply((1:n)[nk != 0], function(i) Ztime2.st[[i]] %*% bi.st[[i]])) # M*GQ matrix #
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix #
  eta.s <- as.vector(Wtime2 %*% phi + alpha * Xtime2 %*% beta) + alpha * Ztime2.b # M*GQ matrix #
  exp.es <- exp(eta.s) # M*GQ matrix #
  const <- matrix(0, n, GQ) # n*GQ matrix #
  const[nk != 0,] <- rowsum(lamb[Index1] * exp.es, Index)  
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  
  f.long <- sapply(1:n, function(i) dmvnorm(Y.st[[i]], as.vector(X.st[[i]] %*% beta), VY[[i]])) 
  # vector of length n #
  return(sum(log(f.long * deno / (pi ^ (ncz / 2)))))
}
