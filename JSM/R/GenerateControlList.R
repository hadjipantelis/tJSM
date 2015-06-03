#' Generate the control variables used by the model 
#' 
#' @param tol.P :  
#' @param tol.L  : 
#' @param max.iter :  
#' @param SE.method :  
#' @param delta :
#' @param nknot :
#' @return VOID


GenerateControlList <- function( control ){

 controlvals <- list(tol.P = 10 ^ (-4), tol.L = 10 ^ (-8), max.iter = 200, SE.method = 'PRES', delta = 10 ^ (- 5), 
                      nknot = 12)
  control <- c(control)
  namec <- names(control)
  if(length(uname <- namec[!namec %in% names(controlvals)]) > 0){
    warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
  }
  controlvals[namec] <- control
  if(controlvals$SE.method == 'PLFD' | controlvals$SE.method == 'PFDS'){
    controlvals$delta <- 10 ^ (- 3)
  }
  controlvals[namec] <- control

  return(controlvals)
  
}
