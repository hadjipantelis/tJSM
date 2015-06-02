result <- MakeResultObj(className, theta.new, call, Vcov, converge, controlvals, time.SE, N, n, d, rho, U,  ncb, alpha.name, phi.names){ 
 theta.new$lamb <- cbind("time" = U, "bashaz" = theta.new$lamb)
  if(className == 'jmodelMult'){ 
    names(theta.new$gamma) <- paste("gamma.", 1:ncb, sep = "")
    names(theta.new$Bsigma) <- "sigma.b"  
  } else {
    names(theta.new$beta) <- beta.names
    if(ncz > 1){
      dimnames(theta.new$Bsigma) <- dimnames(Bsigma)
    } else { 
      names(theta.new$Bsigma) <- "sigma.b"
    }
  }
  names(theta.new$phi) <- phi.names
  names(theta.new$alpha) <- if(model == 1) alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"

  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if(converge == 1) "success" else "failure"
  result$control <- controlvals
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  class(result) <- className
  
  return(result)
  
}
