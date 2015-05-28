
List2Vec <- function (theta) {

  BSigma <- theta$BSigma
  BSigma <- if (is.matrix(BSigma)) BSigma[lower.tri(BSigma, diag = TRUE)] else BSigma
  
  para <- c(theta$beta, theta$phi, theta$alpha, theta$Ysigma, BSigma)
  return(para)
}
