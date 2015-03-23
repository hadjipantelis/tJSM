library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(microbenchmark)

source('logLik.jmodelTM.R')
source('vcov.jmodelTM.R')
source('AIC.jmodelTM.R')
source('BIC.jmodelTM.R')
source('jmodelTM.R')
source('print.jmodelTM.R')
source('summary.jmodelTM.R')
source('print.summary.jmodelTM.R')
source('InitValTM1.R')
source('InitValTM2.R')
source('EMiterTM2.R')
source('Lamb2.R')
source('DQfunc2.R')
source('Sfunc.R')
source('PFDS.R')
source('PLFD.R')
source('PRES.R')
source('Indexing.R')
source('Vec2List.R')
source('List2Vec.R')
source('LH1.R')
source('Lamb1.R')
source('DQfunc1.R')
source('EMiterTM1.R')
source('LH2.R')
source('CheckDelta.R')

source('HelperRcppEigenFunc.cxx');

if(1==2){

load('aids.rda')
fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)

fitJTdam2r1n <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime') 
fitJTdam1r1n <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime') 
fitJTdam2r0n <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime') 
fitJTdam1r0n <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime') 

NewTimings <- microbenchmark( jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime'),  jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime'), times=66 )

source('../Final_Version_IV-test/EMiterTM1.R')
source('../Final_Version_IV-test/EMiterTM2.R')
source('../Final_Version_IV-test/Lamb1.R')
source('../Final_Version_IV-test/Lamb2.R')
source('../Final_Version_IV-test/DQfunc1.R')
source('../Final_Version_IV-test/DQfunc2.R')
source('../Final_Version_IV-test/LH1.R')
source('../Final_Version_IV-test/LH2.R')

fitJTdam2r1o <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime') 
fitJTdam1r1o <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime') 
fitJTdam2r0o <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime') 
fitJTdam1r0o <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime') 

OldTimings <- microbenchmark( jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime'),  jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime') , times=66)

save.image('Second_modelwide_benchmarks.RData.Pantelis')
}

if (1==1){ # This is to be done

load("liver.rda")
liver.unq <- liver[!duplicated(liver$ID), ]
liver.pns <- liver[liver$Trt == "prednisone", ]
liver.pcb <- liver[liver$Trt == "placebo", ]
times.pns <- sort(unique(liver.pns$obstime))
times.pcb <- sort(unique(liver.pcb$obstime))
means.pns <- tapply(liver.pns$proth, liver.pns$obstime, mean)
means.pcb <- tapply(liver.pcb$proth, liver.pcb$obstime, mean)

fitLME <- lme(proth ~ Trt * obstime, random = ~ obstime | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)

fitJTdlm2r1n <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1,timeVarY = 'obstime') 
fitJTdlm1r1n <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1,timeVarY = 'obstime') 
fitJTdlm2r0n <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=0,timeVarY = 'obstime') 
fitJTdlm1r0n <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0,timeVarY = 'obstime')  

NewTimings_Liver <- microbenchmark( jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, liver, model = 2, rho=0,timeVarY = 'obstime'),  jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0,timeVarY = 'obstime'), times=25 )

source('../Final_Version_IV-test/EMiterTM1.R')
source('../Final_Version_IV-test/EMiterTM2.R')
source('../Final_Version_IV-test/Lamb1.R')
source('../Final_Version_IV-test/Lamb2.R')
source('../Final_Version_IV-test/DQfunc1.R')
source('../Final_Version_IV-test/DQfunc2.R')
source('../Final_Version_IV-test/LH1.R')
source('../Final_Version_IV-test/LH2.R')

fitJTdlm2r1o <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1,timeVarY = 'obstime') 
fitJTdlm1r1o <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1,timeVarY = 'obstime') 
fitJTdlm2r0o <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=0,timeVarY = 'obstime') 
fitJTdlm1r0o <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0,timeVarY = 'obstime')  

OldTimings_Liver <- microbenchmark( jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1,timeVarY = 'obstime'), jmodelTM(fitLME, fitCOX, liver, model = 2, rho=0,timeVarY = 'obstime'),  jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0,timeVarY = 'obstime'), times=25 )

save.image('Third_modelwide_benchmarks.RData.Pantelis')

}

