#' Return the BIC score associated with a jmodelTM object
#'
#' Return the BIC score associated with a jmodelTM object
#'
#' @param object A jmodelTM object as this is produced by function jmodelTM
#' @param ... Not used.
#'
#' @export
BIC.jmodelTM <-  function (object, ...) {
  - 2 * object$logLik + log(object$n) * nrow(object$Vcov)
}
