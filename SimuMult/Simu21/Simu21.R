
#=============== Simulation Study for Model II with rho = 1 ===============#

#=============== Steup for running on Gauss ===============#

"ppaste" <- function(...){paste(...,sep="")}

args <- commandArgs(TRUE)

cat(ppaste("Command-line arguments:\n"))
print(args)

###################
sim_start <- 0
###################

if (length(args)==0){
  sim_num <- sim_start + 1
  set.seed(121231)
  sinkit <- FALSE
} else {
  # SLURM can use either 0- or 1-indexing...
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
  set.seed(762*(sim_num-1) + 121231)
}

i <- sim_num
sinkfile <- paste("Output",i,".txt",sep="")

cat(paste("\nAnalyzing dataset number ",i,"...\n\n",sep=""))


#============================== Run the simulation study ==============================#

gauss <- TRUE

if (sinkit){
  cat(paste("Sinking output to: ",sinkfile,"\n",sep=""))
  sink(sinkfile)
}  

library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(splines)

source('jmodelMult.R')
source('InitValMult1.R')
source('InitValMult2.R')
source('EMiterMult1.R')
source('EMiterMult2.R')
source('DQfuncMult1.R')
source('DQfuncMult2.R')
# source('CheckEtaMult.R')
source('CheckDeltaMult.R')
source('SfuncMult.R')
source('LambMult1.R')
source('LambMult2.R')
source('LHMult1.R')
source('LHMult2.R')
source('List2VecMult.R')
source('Vec2ListMult.R')
source('PFDSMult.R')
source('PLFDMult.R')
source('PRESMult.R')
source('print.jmodelMult.R')
source('print.summary.jmodelMult.R')
source('summary.jmodelMult.R')

load("MultData21.RData")

control1 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8)
control2 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8, 
                 SE.method = "PFDS", delta = 10 ^ (- 2))
control3 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8, 
                 SE.method = "PLFD", delta = 10 ^ (- 2))

Data <- MultData21[[i]]
fitLME <- lme(Y ~ bs(ObsTime, knots = 2, Boundary.knots = c(0, 10)), random =~ 1 | ID, data = Data)
fitCOX <- coxph(Surv(start, stop, event) ~ Z, data = Data, x = TRUE)

time1 <- system.time(fit1 <- jmodelMult(fitLME, fitCOX, Data, model = 2, rho = 1, control = control1))[3]
time2 <- system.time(fit2 <- jmodelMult(fitLME, fitCOX, Data, model = 2, rho = 1, control = control2))[3]
time3 <- system.time(fit3 <- jmodelMult(fitLME, fitCOX, Data, model = 2, rho = 1, control = control3))[3]

theta.hat <- fit1$coefficients
step <- fit1$numIter
Var1 <- diag(fit1$Vcov)
Var2 <- diag(fit2$Vcov)
Var3 <- diag(fit3$Vcov)
time.SE1 <- fit1$time.SE
time.SE2 <- fit2$time.SE
time.SE3 <- fit3$time.SE

results = list(time1 = time1, time2 = time2, time3 = time3, time.SE1 = time.SE1, time.SE2 = time.SE2, 
               time.SE3 = time.SE3, theta.hat = theta.hat, Var1 = Var1, Var2 = Var2, Var3 = Var3, step = step)

save(results, file = paste("MultResult21", i, ".RData", sep = ""))
