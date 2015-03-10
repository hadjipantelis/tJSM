
#=============== The S Function for Model I & II of Multiplicative Joint Model ===============#

SfuncMult <- function (model, theta) {
  S <- if(model == 1) DQfuncMult1(theta, theta) else DQfuncMult2(theta, theta)
  return(S)
}