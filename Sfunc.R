
#=============== The S Function for Model I & II ===============#

Sfunc <- function (model, theta) {
  S <- if (model == 1) DQfunc1(theta, theta) else DQfunc2(theta, theta)
  return(S)
}
