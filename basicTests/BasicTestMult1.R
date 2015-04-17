
#========== Initial test of the code for the NMRE model ==========#

library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(splines)
library(Rcpp)
library(RcppEigen)
library(lineprof)
library(microbenchmark)

source('logLik.jmodelMult.R')
source('vcov.jmodelMult.R')
source('AIC.jmodelMult.R')
source('BIC.jmodelMult.R')
source('jmodelMult.R')
source('InitValMultGeneric.R') 
source('EMiterMultGeneric.R') 
source('DQfuncMultGeneric.R')
source('CheckDeltaMult.R')
source('SfuncMult.R')
source('LambMultGeneric.R') 
source('LHMultGeneric.R') 
source('List2VecMult.R')
source('Vec2ListMult.R')
source('PFDSMult.R')
source('PLFDMult.R')
source('PRESMult.R')
source('print.jmodelMult.R')
source('print.summary.jmodelMult.R')
source('summary.jmodelMult.R')

load("aids.rda")

model <- 1
control <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 100, nknot = 4)

# source('HelperRcppEigenFunc.cxx');
cppToCompile <- list.files('src/')
for (i in 1:length(cppToCompile)){
  print( paste0('Compiling: ',cppToCompile[i]) )
  sourceCpp( paste0('src/',cppToCompile[i]) )
}



fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
time <- system.time(fitJT <- jmodelMult(fitLME, fitCOX, aids, model, control = control))[3]


load("MultData1.RData")
Data <- MultData[[1]]
fitLME <- lme(Y ~ bs(ObsTime, knots = 2, Boundary.knots = c(0, 5)), random =~ 1 | ID, data = Data)
fitCOX <- coxph(Surv(start, stop, event) ~ Z, data = Data, x = TRUE)
time <- system.time(fitJT <- jmodelMult(fitLME, fitCOX, Data, model, control = control))[3]


model <- 2
load("MultData2.RData")
Data <- MultData[[1]]
fitLME <- lme(Y ~ bs(ObsTime, knots = 2, Boundary.knots = c(0, 5)), random =~ 1 | ID, data = Data)
fitCOX <- coxph(Surv(start, stop, event) ~ Z, data = Data, x = TRUE)
time <- system.time(fitJT <- jmodelMult(fitLME, fitCOX, Data, model, control = control))[3]

