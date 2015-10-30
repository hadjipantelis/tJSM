#' Return the AIC score associated with a jmodelTM object
#'
#' Return the AIC score associated with a jmodelTM object
#'
#' @param object A jmodelTM object as this is produced by function jmodelTM
#' @param ... Not used.
#'
#' @export
AIC.jmodelTM <-  function (object, ...) {
print('blah blah!')
  - 2 * object$logLik + 2 * nrow(object$Vcov)
}
  
