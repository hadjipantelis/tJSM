#' Return the variance covariance matrix of a jmodelMult object (a joint model with NMRE)
#'
#' Return the variance covariance matrix of a jmodelMult object (a joint model with NMRE)
#'
#' @param object A jmodelMult object as this is produced by function jmodelMult
#' @param ... Not used.
#'
#' @export      
vcov.jmodelMult <-  function (object, ...) {
  object$Vcov
}
