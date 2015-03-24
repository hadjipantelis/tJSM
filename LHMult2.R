
#=============== Function to Calculate the Likelihood Value for Model II ===============#
#=============== Multiplicative Joint Modeling ===============#

LHMult2 <- function (theta) {
  
  gamma <- theta$gamma
  phi <- theta$phi
  alpha <- theta$alpha
  Ysigma2 <- (theta$Ysigma) ^ 2
  Bsigma2 <- (theta$Bsigma) ^ 2
  lamb <- theta$lamb
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))

  # VY <- lapply(1:n, function(i) as.matrix(Bsigma2 * BTg[[i]] %*% t(BTg[[i]]) + Ysigma2 * diag(1, ni[i]))) 
  VY <- lapply(1:n, function(i) calc_VY(BTg[[i]], BSigma2, Ysigma2)) 

  #VB <- lapply(1:n, function(i) as.numeric(Bsigma2 - (Bsigma2 ^ 2) * t(BTg[[i]]) %*% solve(VY[[i]]) %*% BTg[[i]])) 
  VB <-  lapply(1:n, function(i) calc_VB(M_i2 = BTg[[i]], M_i1 = BSigma2, M_i3 = VY[[i]]))

  muB <- lapply(1:n, function(i) as.numeric(1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - BTg[[i]])))
  ## For reference purposes only
  # muB <- lapply(1:n, function(i) as.vector(BSigma %*% t(Z.st[[i]]) %*% solve(VY[[i]]) %*% as.vector(Y.st[[i]] - X.st[[i]] %*% beta))) 
  # muB <- lapply(1:n, function(i) calc_muB( BSigma, M_i3=Z.st[[i]], y_i1=Y.st[[i]],  y_i2=beta, M_i1=VY[[i]], M_i2=X.st[[i]])) 

  bi.st <- lapply(1:n, function(i) as.matrix(muB[[i]] + sqrt(2 * VB[[i]]) * t(b)))
  ## For reference purposes only
  # bi.st <- lapply(1:n, function(i) calc_bi_st(muB[[i]], b ,VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  log.density1 <- log.lamb + as.vector(Ztime %*% phi) + alpha * bi # n*nknot matrix #

  # eta.s <- as.vector(Ztime2 %*% phi) + alpha * bi[Index, ] # M*nknot matrix #
  # exp.es <- exp(eta.s) # M*nknot matrix #
  calc_y_a( Ztime2,phi); # Ztime2 gets altered
  eta.s <- alpha * bi[Index, ] + Ztime2  
  n_ <- ncol(eta.s)
  m_ <- nrow(eta.s)
  exp.es <- matrix(calc_expM(eta.s), m_, n_)


  const <- matrix(0, n, nknot) # n*nknot matrix #
  
  # const[nk != 0, ] <- rowsum(lamb[Index1] * exp.es, Index) # n*nknot matrix # 
  const[nk != 0, ] <- calc_rowsum( (Index),  exp.es * lamb[Index1])

  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 

# log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*GQ matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*nknot matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  
  # f.long <- sapply(1:n, function(i) dmvnorm(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]])) # vector of length n #
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]]))

  return(sum(log(f.long * deno / sqrt(pi))))
}
