
#=============== Simulation Study for Model I with rho = 1 ===============#

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

source('jmodelTM.R')
source('print.jmodelTM.R')
source('summary.jmodelTM.R')
source('print.summary.jmodelTM.R')
source('InitValTM1.R')
source('InitValTM2.R')
source('EMiterTM1.R')
source('EMiterTM2.R')
source('Lamb1.R')
source('Lamb2.R')
source('DQfunc1.R')
source('DQfunc2.R')
source('Sfunc.R')
source('PFDS.R')
source('PLFD.R')
source('PRES.R')
source('Indexing.R')
source('Vec2List.R')
source('List2Vec.R')
source('LH1.R')
source('LH2.R')
source('CheckEta.R')

load("TMData11.RData")

control1 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8)
control2 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8, 
                 SE.method = "PFDS", eta = 10 ^ (- 3))
control3 <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 200, nknot = 8, 
                 SE.method = "PLFD", eta = 10 ^ (- 3))

timeVarY <- "ObsTime"
funcY <- list()
funcY[[1]] <- function(x) x

Data <- TMData11[[i]]
fitLME <- lme(Y ~ X1 + X2 + ObsTime + X1:ObsTime + X2:ObsTime - 1, random = ~ ObsTime | ID, data = Data)
fitCOX <- coxph(Surv(start, stop, event) ~ X1 + X2, data = Data, x = TRUE)

time1 <- system.time(fit1 <- jmodelTM(fitLME, fitCOX, Data, rho = 1, timeVarY = timeVarY, funcY = funcY, control = control1))[3]
time2 <- system.time(fit2 <- jmodelTM(fitLME, fitCOX, Data, rho = 1, timeVarY = timeVarY, funcY = funcY, control = control2))[3]
time3 <- system.time(fit3 <- jmodelTM(fitLME, fitCOX, Data, rho = 1, timeVarY = timeVarY, funcY = funcY, control = control3))[3]

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

save(results, file = paste("TMResult11", i, ".RData", sep = ""))
