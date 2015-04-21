
#========== Initial test of the code for the NMRE model ==========#

library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(splines)
library(lineprof)
library(microbenchmark)
library(RcppEigen)
library(Rcpp)

source('vcov.jmodelMult.R');    source('vcov.jmodelTM.R')
source('logLik.jmodelMult.R');  source('logLik.jmodelTM.R')
source('AIC.jmodelMult.R');     source('BIC.jmodelMult.R')
source('AIC.jmodelTM.R')  ;     source('BIC.jmodelTM.R')
source('jmodelMult.R')    ;     source('jmodelTM.R')
source('InitValMultGeneric.R'); source('InitValTMGeneric.R')
source('EMiterMultGeneric.R');  source('EMiterTMGeneric.R')
source('DQfuncMultGeneric.R');  source('DQfuncGeneric.R')
source('CheckDeltaMult.R');     source('CheckDelta.R')
source('SfuncMult.R');          source('Sfunc.R')
source('LambMultGeneric.R');    source('LambGeneric.R')
source('LHMultGeneric.R');      source('LHGeneric.R')  
source('List2VecMult.R');       source('List2Vec.R')
source('Vec2ListMult.R');       source('Vec2List.R')
source('PFDSMult.R');           source('PFDS.R')
source('PLFDMult.R');           source('PLFD.R')
source('PRESMult.R');           source('PRES.R')
source('print.jmodelMult.R');   source('print.jmodelTM.R')
source('summary.jmodelMult.R'); source('summary.jmodelTM.R')

source('print.summary.jmodelMult.R')
source('print.summary.jmodelTM.R')
source('Indexing.R')


cppToCompile <- list.files('src/')
for (i in 1:length(cppToCompile)){
  print( paste0('Compiling: ',cppToCompile[i]) )
  sourceCpp( paste0('src/',cppToCompile[i]) )
}

load("aids.rda")
fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list( max.iter = 100, nknot = 5)
 
AIDS_MULT_RESULT_NEW <- microbenchmark( times = 100,  fitJTm1r0_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho = 0, control = control),  
fitJTm1r1_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho = 1, control = control), 
fitJTm2r0_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho = 0, control = control),
fitJTm2r1_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho = 1, control = control))

save(file='AIDS_MULT_RESULT_NEW.RData',AIDS_MULT_RESULT_NEW)


load("pbc.rda")
fitLME <- lme(log(serBilir) ~ bs(obstime, df = 6, degree = 2), random = ~ 1 | ID, data = pbc)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = pbc, x = TRUE) 
control <- list(tol.P = 10 ^ (- 3), max.iter = 100, nknot = 10)

PBC_MULT_RESULT_NEW <- microbenchmark( times = 70, fitJTm1r0ph_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 1, rho = 0, control = control),
fitJTm1r1ph_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 1, rho = 1, control = control), 
fitJTm2r0ph_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho = 0, control = control),  
fitJTm2r1ph_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho = 1, control = control))  

save(file='PBC_MULT_RESULT_NEW.RData',PBC_MULT_RESULT_NEW)


load('aids.rda')
fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list(nknot = 15)

AIDS_TM_RESULT_NEW <- microbenchmark( times = 100, fitJTdam1r0n_TM <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime', control = control),
fitJTdam1r1n_TM <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime', control = control),
fitJTdam2r0n_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime', control = control),
fitJTdam2r1n_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime', control = control))
 
save(file='AIDS_TM_RESULT_NEW.RData', AIDS_TM_RESULT_NEW)


load("liver.rda")
fitLME <- lme(proth ~ Trt * obstime, random = ~ obstime | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE) 
control <- list(tol.P = 10 ^ (- 3))
 
LIVER_TM_RESULT_NEW <- microbenchmark( times = 70, fitJTm1r0po_TM <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0, timeVarY = 'obstime', control = control),
fitJTm1r1po_TM <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1, timeVarY = 'obstime', control = control),
fitJTm2r0po_TM <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=0, timeVarY = 'obstime', control = control),
fitJTm2r0po_TM <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1, timeVarY = 'obstime', control = control))
 
save(file='LIVER_TM_RESULT_NEW.RData',LIVER_TM_RESULT_NEW)


