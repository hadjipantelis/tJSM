#' Print the logLikelihood of joint model model
#'
#' Print the logLikelihood of joint model model
#'
#' @param object A jmodelTM object.
#' @param ... Not used.
#'
#' @export

logLik.jmodelTM <-  function (object, ...) {
    if (!inherits(object, "jmodelTM"))
      stop("Only used for 'jmodelTM' objects.\n")
    out <- object$logLik
    attr(out, "df") <- nrow(object$Vcov)
    attr(out, "n") <- object$n
    class(out) <- "logLik"
    out
  }
