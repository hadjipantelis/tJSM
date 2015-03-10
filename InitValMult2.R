
#=============== Initial Value Calculation for Model II with NMRE ===============#

InitValMult2 <- function (gamma) {
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  G <- unlist(lapply(1:n, function(i) as.vector(Y.st[[i]] - BTg[[i]]) %x% as.vector(Y.st[[i]] - BTg[[i]])))
  tempCov1 <- unlist(lapply(BTg, function(x) x %x% x))
  tempCov2 <- unlist(lapply(1:n, function(i) c(diag(1, ni[i]))))
  tempCov <- cbind(tempCov1, tempCov2)
  sigmas <- solve(t(tempCov) %*% tempCov) %*% t(tempCov) %*% G
  Bsigma2 <- sigmas[1]
  Ysigma2 <- sigmas[2]
  if(Bsigma2 < 0) Bsigma2 <- 0.1
  if(Ysigma2 < 0) Ysigma2 <- 0.1
  
  VY <- lapply(1:n, function(i) as.matrix(Bsigma2 * BTg[[i]] %*% t(BTg[[i]]) + Ysigma2 * diag(1, ni[i])))
  bBLUP <- unlist(lapply(1:n, function(i) 1 + Bsigma2 * t(BTg[[i]]) %*% solve(VY[[i]]) %*% 
                                          as.vector(Y.st[[i]] - BTg[[i]])))
  
  #========== fit the Cox model ==========#
  rand <- bBLUP[ID]
  data.init <- data.frame(start = start, stop = stop, event = event, Z = Z, rand = rand)
  fit <- coxph(Surv(start, stop, event) ~ Z + rand, data = data.init)
  phi.old <- fit$coefficients[1 : ncz]
  alpha.old <- fit$coefficients[ncz + 1]
  temp <- as.vector(exp(Ztime2 %*% phi.old + alpha.old * bBLUP[Index])) # M*1 vector #
  lamb.old <- Index2 / as.vector(tapply(temp, Index1, sum)) # vector of length n_u #
  
  if (rho == 0) {
    phi.new <- phi.old
    alpha.new <- alpha.old
    lamb.new <- lamb.old
  } else {
    for (it in 1:iter) {
      exp.es <- exp(as.vector(Ztime2 %*% phi.old + alpha.old * bBLUP[Index]))
      const <- rep(0, n)
      const[nk != 0] <- as.vector(tapply(lamb.old[Index1] * exp.es, Index, sum)) # vector of length n #
      CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|Oi), vector of length n #
      CondExp2 <- CondExp[nk != 0]
      temp1 <- lapply(1:ncz, function(i) CondExp2 * as.vector(tapply(Ztime2[, i] * exp.es * lamb.old[Index1], 
                                         Index, sum)))
      temp1 <- sapply(temp1, sum) # vector of length ncz #
      temp2 <- sum(CondExp2 * as.vector(tapply(bBLUP[Index] * exp.es * lamb.old[Index1], Index, sum))) # scalar #
      temp3 <- lapply(1:(ncz ^ 2), function(i) CondExp2 * as.vector(tapply(Ztime22[, i] * exp.es * 
                                               lamb.old[Index1], Index, sum)))
      temp3 <- sapply(temp3, sum) # vector of length ncz^2 #
      temp4 <- sum(CondExp2 * as.vector(tapply(bBLUP[Index] ^ 2 * exp.es * lamb.old[Index1], Index, sum))) # scalar #
      temp5 <- lapply(1:ncz, function(i) CondExp2 * as.vector(tapply(Ztime2[, i] * bBLUP[Index] * exp.es * 
                                         lamb.old[Index1], Index, sum)))
      temp5 <- sapply(temp5, sum) # vector of length ncz #
      
      phiScore <- colSums(d * Ztime) - temp1 # vector of length ncz #
      alphaScore <- sum(d * bBLUP) - temp2
      pa.score <- c(phiScore, alphaScore)
      pa.info <- matrix(0, (ncz + 1), (ncz + 1)) # (ncz+1)*(ncz+1) matrix #
      pa.info[1:ncz, 1:ncz] <- - temp3
      pa.info[(ncz + 1), (ncz + 1)] <- - temp4
      pa.info[(ncz + 1), 1:ncz] <- - temp5
      pa.info[1:ncz, (ncz + 1)] <- - temp5
      
      #=============== Update phi and alpha ===============#
      pa.old <- c(phi.old, alpha.old) # vector of length (ncz+1) #
      paSVD <- svd(pa.info)
      pa.info.inv <- paSVD$v %*% diag(1 / paSVD$d) %*% t(paSVD$u)
      pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncz+1) #
      phi.new <- pa.new[1 : ncz]
      alpha.new <- pa.new[ncz + 1]
      
      #========== Calculate the new lambda with new parameters ==========#
      exp.esn <- exp(as.vector(Ztime2 %*% phi.new + alpha.new * bBLUP[Index]))
      tempLamb <- as.vector(tapply(CondExp[Index] * exp.esn, Index1, sum)) # vector of length M #
      lamb.new <- Index2 / tempLamb
      
      #========== Check Convergence ==========#
      err <- max(abs(pa.new - pa.old) / (abs(pa.old) + tol.P))
      if (err <= tol.P) break
      else {
        phi.old <- phi.new
        alpha.old <- alpha.new
        lamb.old <- lamb.new
      }
    }
  }
  
  result <- list(phi = phi.new, alpha = alpha.new, lamb = lamb.new, Ysigma = sqrt(Ysigma2), 
                 Bsigma = sqrt(Bsigma2))
  return(result)
}