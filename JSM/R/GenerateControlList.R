#' Generate the control variables used by the model 
#' 
#' @param tol.P : tolerance for parameters (numeric)
#' @param tol.L  : tolerance for log-likelihood  (numeric)
#' @param max.iter : maximum number of EM iterations (numeric)
#' @param SE.method : standard error estimation method  (char)
#' @param delta : increment used by forward difference  (numeric)
#' @param nknot : number of Gauss-Hermite quadrature knots  (numeric)
#' @return controlVals list

GenerateControlList <- function( control ){

 controlVals <- list(tol.P = 10 ^ (-4), tol.L = 10 ^ (-8), max.iter = 200, SE.method = 'PRES', delta = 10 ^ (- 5), 
                      nknot = 12)
  control <- c(control)
  namec <- names(control)
  if(length(uname <- namec[!namec %in% names(controlVals)]) > 0){
    warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
  }
  controlVals[namec] <- control
  if(controlVals$SE.method == 'PLFD' | controlVals$SE.method == 'PFDS'){
    controlVals$delta <- 10 ^ (- 3)
  }
  controlVals[namec] <- control

  return(controlVals)
  
}
