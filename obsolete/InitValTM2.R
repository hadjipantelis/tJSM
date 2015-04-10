
#=============== Initial Value Calculation for Transformation Model II ===============#

InitValTM2 <- function (beta) {
  
  rand <- rowSums(Z * bBLUP[ID, ]) # vector of length N #
  rand.time <- rowSums(Ztime * bBLUP) # vector of length n #
  rand.time2 <- rowSums(Ztime2 * bBLUP[Index, ]) # vector of length M #
  
  #========== first fit the Cox model ==========#
  data.init <- data.frame(start = start, stop = stop, event = event, W = W, rand = rand)
  fit <- coxph(Surv(start, stop, event) ~ W + rand, data = data.init) 
  phi.old <- fit$coefficients[1:ncw]
  alpha.old <- fit$coefficients[ncw + 1]
  temp <- as.vector(exp(Wtime2 %*% phi.old + alpha.old * rand.time2)) # M*1 vector #
  lamb.old <- Index2 / calc_tapply_vect_sum(  v1=temp, v2=  as.integer(Index1-1)) # vector of length n_u #
  
  if (rho == 0) {
    phi.new <- phi.old
    alpha.new <- alpha.old
    lamb.new <- lamb.old
  } else {
    for (it in 1:iter) {
      exp.es <- exp(as.vector(Wtime2 %*% phi.old + alpha.old * rand.time2))
      const <- rep(0, n)
      temp0a <- exp.es * lamb.old[Index1];
      const[nk != 0] <- calc_tapply_vect_sum(  v1=temp0a, v2=  as.integer(Index-1)) # vector of length n #

      CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|Oi), vector of length n #
      CondExp2 <- CondExp[nk != 0]

      temp0b <- rand.time2 *temp0a;

      temp1 <- sum(CondExp2 *calc_tapply_vect_sum(  v1= temp0b, v2=  as.integer(Index-1)));
      temp2 <- sum(CondExp2 *calc_tapply_vect_sum(  v1=rand.time2 * temp0b, v2=  as.integer(Index-1)));
      temp3 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1=Wtime2[, i] * temp0a, v2=  as.integer(Index-1)))
      temp3 <- sapply(temp3, sum) # vector of length ncw #
      temp4 <-  lapply(1:ncw^2, function(i) CondExp2 * calc_tapply_vect_sum(  v1=Wtime22[, i] * temp0a, v2=  as.integer(Index-1)))
      temp4 <- sapply(temp4, sum) # vector of length ncw^2 #
      temp5 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1=Wtime2[, i] * temp0b, v2=  as.integer(Index-1)))
      temp5 <- sapply(temp5, sum) # vector of length ncw #  


      phiScore <- colSums(d * Wtime) - temp3 # vector of length ncw #
      alphaScore <- sum(d * rand.time) - temp1
      pa.score <- c(phiScore, alphaScore)
      pa.info <- matrix(0, (ncw + 1), (ncw + 1)) # (ncw+1)*(ncw+1) matrix #
      pa.info[1:ncw, 1:ncw] <- - temp4
      pa.info[(ncw + 1), 1:ncw] <- - temp5
      pa.info[1:ncw, (ncw + 1)] <- - temp5
      pa.info[(ncw + 1), (ncw + 1)] <- - temp2
      
      #=============== Update phi and alpha ===============#
      pa.old <- c(phi.old, alpha.old) # vector of length (ncw+1) #
      paSVD <- svd(pa.info)
      pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
      pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
      phi.new <- pa.new[1:ncw]
      alpha.new <- pa.new[ncw + 1]
      
      #========== Calculate the new lambda with new parameters ==========#
      exp.esn <- exp(as.vector(Wtime2 %*% phi.new + alpha.new * rand.time2))
      tempLamb <- calc_tapply_vect_sum(  v1=CondExp[Index] * exp.esn, v2=  as.integer(Index1-1))
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
  
  result <- list(phi = phi.new, alpha = alpha.new, lamb = lamb.new)
  return(result)
}
