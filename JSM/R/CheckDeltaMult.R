#' Check whether the delta is too large when Richardson Extrapolation is used
#'
#' Check whether the delta is too large when Richardson Extrapolation is used
#'
#' @param theta List of parameters returned by the EM iteration procedure  
#' @param delta increment used by (forward/Richardson Extrapolation) difference algorithm
#'

CheckDeltaMult <- function (theta, delta) {
  
  Ysigma <- theta$Ysigma
  check1 <- (Ysigma - 2 * delta > 0)
  
  Bsigma <- theta$Bsigma
  check2 <- (Bsigma - 2 * delta > 0)

  return(all(c(check1, check2))) 
}
