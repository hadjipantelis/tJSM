#' Print the logLikelihood of joint model model with NMRE
#'
#' Print the logLikelihood of joint model model with NMRE
#'
#' @param object A jmodelMult object
#' @param ... Not used.
#'
#' @export

logLik.jmodelMult <-  function (object, ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Only used for 'jmodelMult' objects.\n")
  out <- object$logLik
  attr(out, "df") <- nrow(object$Vcov)
  attr(out, "n") <- object$n
  class(out) <- "logLik"
  out
}
