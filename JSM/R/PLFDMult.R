
#================== Profile Likelihood Method with Forward Difference for Model I & II ==================#

<<<<<<< HEAD
PLFDMult <- function (model, theta, tol, iter, delta, ncz, ncb, alpha.name, phi.names, B.st, n, Y.st,  b, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2) {
  
  pl <- LHMultGeneric(theta = theta, B.st =B.st , n = n, Y.st = Y.st,  b = b, model=model, Btime = Btime, Btime2 = Btime2, Index = Index, Ztime = Ztime, Ztime2 = Ztime2, Index0 = Index0, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ) / n
=======
PLFDMult <- function (model, theta, tol, iter, delta, B.st, n, Y.st, b, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ, ncz, ncb, Index2, alpha.name, phi.names) {
  
  pl <- LHMultGeneric(theta, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ) / n
>>>>>>> Just_Before_Merge
  
  para <- List2VecMult(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  PLs <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      para1 <- para
      para1[i] <- para[i] + delta
      para1[j] <- para1[j] + delta
      result <- LambMultGeneric(para1, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2)
      para1.list <- Vec2ListMult(para1, ncz, ncb)
      theta.input1 <- list(gamma = para1.list$gamma, phi = para1.list$phi, alpha = para1.list$alpha, 
                           Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
<<<<<<< HEAD
      PLs[i, j] <- LHMultGeneric(theta.input1, B.st =B.st , n = n, Y.st = Y.st,  b = b, model=model, Btime = Btime, Btime2 = Btime2, Index = Index, Ztime = Ztime, Ztime2 = Ztime2, Index0 = Index0, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ)/n
=======
      PLs[i, j] <- LHMultGeneric(theta.input1, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ)/n
>>>>>>> Just_Before_Merge
    }
  }
  
  pls <- rep(0, len)
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
<<<<<<< HEAD
    result <- LambMultGeneric(para1, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st,  b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2)  
    para1.list <- Vec2ListMult(para1, ncz, ncb)
    theta.input1 <- list(gamma = para1.list$gamma, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    pls[i] <- LHMultGeneric(theta.input1, B.st =B.st , n = n, Y.st = Y.st,  b = b, model=model, Btime = Btime, Btime2 = Btime2, Index = Index, Ztime = Ztime, Ztime2 = Ztime2, Index0 = Index0, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ) / n  
=======
    result <- LambMultGeneric(para1, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2)  
    para1.list <- Vec2ListMult(para1, ncz, ncb)
    theta.input1 <- list(gamma = para1.list$gamma, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    pls[i] <- LHMultGeneric(theta.input1, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ) / n  
>>>>>>> Just_Before_Merge
  }
  
  I <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      I[i, j] <- - (PLs[i,j] - pls[i] - pls[j] + pl) / (delta ^ 2)
    }
  }
  I <- I + t(I) - diag(diag(I)) 
  V <- solve(I) / n;
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste("gamma.", 1:ncb, sep = ""), paste(rep("phi:", ncz), phi.names, sep = ""), Valpha.name, 
              "sigma.e", "sigma.b")
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
