
library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(microbenchmark)
library(lineprof)
library(Rcpp)
library(RcppEigen)

source('logLik.jmodelTM.R')
source('vcov.jmodelTM.R')
source('AIC.jmodelTM.R')
source('BIC.jmodelTM.R')
source('jmodelTM.R')
source('print.jmodelTM.R')
source('summary.jmodelTM.R')
source('print.summary.jmodelTM.R')
source('InitValTMGeneric.R')
source('EMiterTMGeneric.R')
source('LambGeneric.R')
source('DQfuncGeneric.R')
source('Sfunc.R')
source('PFDS.R')
source('PLFD.R')
source('PRES.R')
source('Indexing.R')
source('Vec2List.R')
source('List2Vec.R')
source('LHGeneric.R')
source('CheckDelta.R')

# source('HelperRcppEigenFunc.cxx');
cppToCompile <- list.files('src/')
for (i in 1:length(cppToCompile)){
  print( paste0('Compiling: ',cppToCompile[i]) )
  sourceCpp( paste0('src/',cppToCompile[i]) )
}


load('aids.rda')
fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)

fitJTdam2r1n <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime') 
fitJTdam1r1n <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime') 
fitJTdam2r0n <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime') 
fitJTdam1r0n <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime') 

NewTimingsGeneric <- microbenchmark( jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime'),  jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime'), times=66 )

save.image('Generic_modelwide_benchmarks.RData.Pantelis')

