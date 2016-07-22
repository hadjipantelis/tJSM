
Vec2ListMult <- function (para, ncz, ncb) {
  
  gamma <- para[1 : ncb]
  phi <- para[(ncb + 1) : (ncb + ncz)]
  alpha <- para[ncb + ncz + 1]
  Ysigma <- para[ncz + ncb + 2]
  Bsigma <- para[ncz + ncb + 3]

  result <- list(gamma = gamma, phi = phi, alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma)
  return(result)
}

