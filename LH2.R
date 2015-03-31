
#=============== Function to Calculate the Likelihood Value for Model II ===============#

LH2 <- function (theta) {
  
  beta <- theta$beta
  Ysigma2 <- (theta$Ysigma) ^ 2
  BSigma <- theta$BSigma
  phi <- theta$phi
  alpha <- theta$alpha
  lamb <- theta$lamb
  
  # VY <- lapply(1:n, function(i) as.matrix(Z.st[[i]] %*% BSigma %*% t(Z.st[[i]]) + Ysigma2 * diag(1, ni[i])))
  # VB <- lapply(1:n, function(i) BSigma - BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% Z.st[[i]] %*% BSigma)
  # muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta)))
  # bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2) * solve(chol(solve(VB[[i]]))) %*% t(b)))

  VY <- lapply(1:n, function(i) calc_VY(Z.st[[i]], BSigma, Ysigma2)) 
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = Z.st[[i]], M_i1 = BSigma, M_i3 = VY[[i]]))
  muB <- lapply(1:n, function(i) calc_muB( BSigma, M_i3=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta, M_i1=VY[[i]], M_i2=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 

  # each element is ncz*GQ matrix #
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  # Ztime2.b <- do.call(rbind, lapply((1:n)[nk != 0], function(i) Ztime2.st[[i]] %*% bi.st[[i]])) # M*GQ matrix #
  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk !=      0] - 1)# M*GQ matrix # 

  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Wtime %*% phi) + alpha * Ztime.b # n*GQ matrix #

  # eta.s <- as.vector(Wtime2 %*% phi) + alpha * Ztime2.b # M*GQ matrix #
  # exp.es <- exp(eta.s) # M*GQ matrix #
  calc_y_a( Ztime2.b,alpha); # Ztime2.b gets altered
  exp.es<- as.numeric(Wtime2 %*% phi ) + Ztime2.b  
  calc_expM2(exp.es) 

  const <- matrix(0, n, GQ) # n*GQ matrix #

  #const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index)  
  const[nk != 0, ] <- calc_rowsum( (Index),  exp.es * lamb[Index1])

  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  
  # log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix #   
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  
  # f.long <- sapply(1:n, function(i) dmvnorm(Y.st[[i]], as.vector(X.st[[i]] %*% beta), VY[[i]])) 
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta), VY[[i]]))

  # vector of length n #
  return(sum(log(f.long * deno / (pi ^ (ncz / 2)))))
}
