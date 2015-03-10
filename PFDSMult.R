
#========== Differentiate the S function with Forward Difference for Model I & II ==========#
#========== Multipicative Joint Modeling ==========#

PFDSMult <- function (model, theta, tol, iter, delta) {
  
  S <- SfuncMult(model, theta)
  
  para <- List2VecMult(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- if(model == 1) LambMult1(para1, lamb.init, tol, iter) else LambMult2(para1, lamb.init, tol, iter)
    para1.list <- Vec2ListMult(para1, ncz, ncb)
    theta.input1 <- list(gamma = para1.list$gamma, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    S1 <- SfuncMult(model, theta.input1)
    DS[i, ] <- (S1 - S) / delta
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2
  svd_DS <- svd(DS)
  V <- - svd_DS$v %*% diag(1 / svd_DS$d) %*% t(svd_DS$u)
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste("gamma.", 1:ncb, sep = ""), paste(rep("phi:", ncz), phi.names, sep = ""), Valpha.name, 
              "sigma.e", "sigma.b")
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
