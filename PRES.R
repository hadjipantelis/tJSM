
#========== Differentiate the S function with Richardson Extrapolation for Model I & II ==========#

PRES <- function (model, theta, tol, iter, delta) {     
  
  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para2 <- para3 <- para4 <- para
    para1[i] <- para[i] - 2 * delta
    para2[i] <- para[i] - delta
    para3[i] <- para[i] + delta
    para4[i] <- para[i] + 2 * delta
    result1 <- if (model == 1) Lamb1(para1, lamb.init, tol, iter) else Lamb2(para1, lamb.init, tol, iter)
    result2 <- if (model == 1) Lamb1(para2, lamb.init, tol, iter) else Lamb2(para2, lamb.init, tol, iter)
    result3 <- if (model == 1) Lamb1(para3, lamb.init, tol, iter) else Lamb2(para3, lamb.init, tol, iter)
    result4 <- if (model == 1) Lamb1(para4, lamb.init, tol, iter) else Lamb2(para4, lamb.init, tol, iter)
    list1 <- Vec2List(para1, ncx, ncz, ncw)
    list2 <- Vec2List(para2, ncx, ncz, ncw)
    list3 <- Vec2List(para3, ncx, ncz, ncw)
    list4 <- Vec2List(para4, ncx, ncz, ncw)
    theta.input1 <- list(beta = list1$beta, phi = list1$phi, alpha = list1$alpha, Ysigma = list1$Ysigma, BSigma = list1$BSigma, lamb = result1$lamb)
    theta.input2 <- list(beta = list2$beta, phi = list2$phi, alpha = list2$alpha, Ysigma = list2$Ysigma, BSigma = list2$BSigma, lamb = result2$lamb)
    theta.input3 <- list(beta = list3$beta, phi = list3$phi, alpha = list3$alpha, Ysigma = list3$Ysigma, BSigma = list3$BSigma, lamb = result3$lamb)
    theta.input4 <- list(beta = list4$beta, phi = list4$phi, alpha = list4$alpha, Ysigma = list4$Ysigma, BSigma = list4$BSigma, lamb = result4$lamb)
    S1 <- Sfunc(model, theta.input1)
    S2 <- Sfunc(model, theta.input2)
    S3 <- Sfunc(model, theta.input3)
    S4 <- Sfunc(model, theta.input4)
    DS[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta)
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2
  V <- -solve(DS);
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Valpha.name, "sigma.e", paste("BSigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
