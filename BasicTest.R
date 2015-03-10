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
source('CheckDelta.R')
source('HelperRcppEigenFunc.cxx');

load('aids.rda')

fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)

library(lineprof)


print('Running a lineprof with fitJT.ph2');
ll <- lineprof(fitJT.ph2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, timeVarY = 'obstime'))

