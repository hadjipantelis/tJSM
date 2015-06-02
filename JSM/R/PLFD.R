
#================== Profile Likelihood Method with Forward Difference for Model I & II ==================#

PLFD <- function (model, theta, tol, iter, delta, n, ncx, ncz, ncw, alpha.name, beta.names, phi.names, p) {
  
  pl <- LHGeneric(theta) / n
  
  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  PLs <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      para1 <- para
      para1[i] <- para[i] + delta
      para1[j] <- para1[j] + delta
      result <- LambGeneric(para1, lamb.init, tol, iter)
      para1.list <- Vec2List(para1, ncx, ncz, ncw)
      theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                           Ysigma = para1.list$Ysigma, BSigma = para1.list$BSigma, lamb = result$lamb)
      PLs[i, j] <- LHGeneric(theta.input1) / n
    }
  }
  
  pls <- rep(0, len)
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- LambGeneric(para1, lamb.init, tol, iter)
    para1.list <- Vec2List(para1, ncx, ncz, ncw)
    theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, BSigma = para1.list$BSigma, lamb = result$lamb)
    pls[i] <- LHGeneric(theta.input1) / n
  }
  
  I <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      I[i, j] <- - (PLs[i, j] - pls[i] - pls[j] + pl) / (delta ^ 2)
    }
  }
  I <- I + t(I) - diag(diag(I))
  # svd_I <- svd(I)
  # V <- svd_I$v %*% diag(1 / svd_I$d) %*% t(svd_I$u) / n ##CHECK THIS TWICE
  V <- solve(I) / n;
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Valpha.name, "sigma.e", paste("BSigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
