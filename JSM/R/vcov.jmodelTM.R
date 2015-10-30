#' Return the variance covariance matrix of a jmodelTM object
#'
#' Return the variance covariance matrix of a jmodelTM object
#'
#' @param object A jmodelTM object as this is produced by function jmodelTM
#' @param ... Not used.
#'
#' @export
vcov.jmodelTM <-  function (object, ...) {
    object$Vcov
  }
