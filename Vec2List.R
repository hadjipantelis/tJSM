
Vec2List <- function (para, ncx, ncz, ncw) {
  
  p <- ncz * (ncz + 1) / 2
  beta <- para[1 : ncx]
  phi <- para[(ncx + 1) : (ncx + ncw)]
  alpha <- para[ncx + ncw + 1]
  Ysigma <- para[ncx + ncw + 2]
  BSigma <- para[(ncx + ncw + 3):(ncx + ncw + p + 2)]
  if (ncz > 1) {
    BSigma.new <- matrix(0, ncz, ncz)
    BSigma.new[lower.tri(BSigma.new, diag = TRUE)] <- BSigma
    BSigma.new <- BSigma.new + t(BSigma.new) - diag(diag(BSigma.new))
    BSigma <- BSigma.new
  }
  result <- list(beta = beta, phi = phi, alpha = alpha, Ysigma = Ysigma, BSigma = BSigma)
  return(result)
}

