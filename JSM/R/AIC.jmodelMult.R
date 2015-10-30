#' Return the AIC score associated with a jmodelMult object (a joint model with NMRE)
#'
#' Return the AIC score associated with a jmodelMult object (a joint model with NMRE)
#'
#' @param object A jmodelMult object as this is produced by function jmodelMult
#' @param ... Not used.
#'
#' @export
AIC.jmodelMult <-  function (object, ...) {
  - 2 * object$logLik + 2 * nrow(object$Vcov)
}
