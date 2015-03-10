
#========== Differentiate the S function with Forward Difference for Model I & II ==========#

PFDS <- function (model, theta, tol, iter, delta) {
  
  S <- Sfunc(model, theta)
  
  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- if (model == 1) Lamb1(para1, lamb.init, tol, iter) else Lamb2(para1, lamb.init, tol, iter)
    para1.list <- Vec2List(para1, ncx, ncz, ncw)
    theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, BSigma = para1.list$BSigma, lamb = result$lamb)
    S1 <- Sfunc(model, theta.input1)
    DS[i, ] <- (S1 - S) / delta
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2
  svd_DS <- svd(DS)
  V <- - svd_DS$v %*% diag(1/svd_DS$d) %*% t(svd_DS$u)
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Valpha.name, "sigma.e", paste("BSigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
