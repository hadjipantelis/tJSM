#' Return the BIC score associated with a jmodelMult object (a joint model with NMRE)
#'
#' Return the BIC score associated with a jmodelMult object (a joint model with NMRE)
#'
#' @param object A jmodelMult object as this is produced by function jmodelMult
#' @param ... Not used.
#'
#' @export
BIC.jmodelMult <-  function (object, ...) {
  - 2 * object$logLik + log(object$n) * nrow(object$Vcov)
}
