
#=============== The S Function for Model I & II ===============#

Sfunc <- function (model, theta) {
  S <-  DQfuncGeneric(theta, theta)
  return(S)
}
