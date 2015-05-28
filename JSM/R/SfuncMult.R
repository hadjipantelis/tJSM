
#=============== The S Function for Model I & II of Multiplicative Joint Model ===============#

SfuncMult <- function (model, theta) {
  S <-  DQfuncMultGeneric(theta, theta)
  return(S)
}
