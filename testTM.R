
#========== Initial test of the code for the TM model ==========#

library(statmod)
library(mvtnorm)
library(nlme)
library(survival)

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



load('aids.rda')
fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), 
              random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)

fitJT.ph <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime')
summary(fitJT.ph)
fitJT.ph2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, timeVarY = 'obstime')
summary(fitJT.ph2)
fitJT.ph3 <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime', control = list(SE.method = 'PFDS'))
summary(fitJT.ph3)

fitJT.po <- jmodelTM(fitLME, fitCOX, aids, rho = 1, timeVarY = 'obstime')
summary(fitJT.po)
fitJT.po2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho = 1, timeVarY = 'obstime')
summary(fitJT.po2)

fitLME2 <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), 
               random = ~ obstime | ID, data = aids)
fitCOX2 <- coxph(Surv(start, stop, event) ~ drug + as.numeric(drug):obstime, data = aids, x = TRUE)
fitJT.ph4 <- jmodelTM(fitLME2, fitCOX2, aids, timeVarY = 'obstime', timeVarT = 'obstime')

load("liver.rda")
liver.unq <- liver[!duplicated(liver$ID), ]
liver.pns <- liver[liver$Trt == "prednisone", ]
liver.pcb <- liver[liver$Trt == "placebo", ]
times.pns <- sort(unique(liver.pns$obstime))
times.pcb <- sort(unique(liver.pcb$obstime))
means.pns <- tapply(liver.pns$proth, liver.pns$obstime, mean)
means.pcb <- tapply(liver.pcb$proth, liver.pcb$obstime, mean)
par(mfrow = c(1, 2))
plot(lowess(times.pns, means.pns, f = .3), type = "l", col = "blue", lwd = 2,
     ylab = "Prothrombin index", xlab = "Time (years)", ylim = c(50, 150),
     main = "Smoothed mean prothrombin index", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(lowess(times.pcb, means.pcb, f = .3), col = "red", lty = 2, lwd = 2)
legend(8, 150, c("Prednisone", "Placebo"), col=c("blue", "red"), lty = 1:2, lwd = 2, 
       bty = "n", cex = 1.5)
plot(survfit(Surv(Time, death) ~ Trt, data = liver.unq), mark.time = FALSE, 
     col=c("red", "blue"), lty = c(2, 1), ylab = "Survival", xlab = "Time (years)",lwd = 2,
     main = "Kaplan-Meier survival curves", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
legend(4, 1, legend = c("Prednisone", "Placebo"), col=c("blue", "red"), lty = 1:2, lwd = 2, 
       bty = "n", cex = 1.5)

fitLME <- lme(proth ~ Trt * obstime, random = ~ obstime | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)
fitJT.ph <- jmodelTM(fitLME, fitCOX, liver, timeVarY = 'obstime')
summary(fitJT.ph)
beta = fitJT.ph$coefficients$beta
plot(lowess(times.pns, means.pns, f = .3), type = "l", col = "blue", lwd = 2,
     ylab = "Prothrombin index", xlab = "Time (years)", ylim = c(50, 150))
lines(lowess(times.pcb, means.pcb, f = .3), col = "red", lty = 2, lwd = 2)
abline(a = beta[1], b = beta[3], col = "red", lty = 2)
abline(a = (beta[1] + beta[2]), b = (beta[3] + beta[4]), col = "blue")
legend(8, 150, c("Prednisone", "Placebo"), col=c("blue", "red"), lty = 1:2, lwd = 2, 
       bty = "n")

fitJT.po <- jmodelTM(fitLME, fitCOX, liver, rho = 1, timeVarY = 'obstime')
summary(fitJT.po)

load("liverResults.RData")
theta.ph <- liverResults$fitJT.ph$coefficients
theta.po <- liverResults$fitJT.po$coefficients
jump <- theta.ph$lamb[, 1]
beta.ph <- theta.ph$beta
beta.po <- theta.po$beta
phi.ph <- theta.ph$phi
phi.po <- theta.po$phi
alpha.ph <- theta.ph$alpha
alpha.po <- theta.po$alpha
lamb.ph <- theta.ph$lamb[, 2]
lamb.po <- theta.po$lamb[, 2]
mu.ph.pns <- beta.ph[1] + beta.ph[2] + (beta.ph[3] + beta.ph[4]) * jump
mu.ph.pcb <- beta.ph[1] + beta.ph[3] * jump
mu.po.pns <- beta.po[1] + beta.po[2] + (beta.po[3] + beta.po[4]) * jump
mu.po.pcb <- beta.po[1] + beta.po[3] * jump
Lamb.ph.pns <- cumsum(lamb.ph * exp(phi.ph + alpha.ph * mu.ph.pns))
Lamb.ph.pcb <- cumsum(lamb.ph * exp(alpha.ph * mu.ph.pcb))
Lamb.po.pns <- log(1 + cumsum(lamb.po * exp(phi.po + alpha.po * mu.po.pns)))
Lamb.po.pcb <- log(1 + cumsum(lamb.po * exp(alpha.po * mu.po.pcb)))
Surv.ph.pns <- exp(-Lamb.ph.pns); Surv.ph.pcb <- exp(-Lamb.ph.pcb)
Surv.po.pns <- exp(-Lamb.po.pns); Surv.po.pcb <- exp(-Lamb.po.pcb)
par(mfrow = c(1, 2))
plot(survfit(Surv(Time, death) ~ Trt, data = liver.unq), mark.time = FALSE,
     col=c("red", "blue"), lty = c(2, 1), ylab = "Survival", xlab = "Time (years)",lwd = 2,
     main = "Kaplan-Meier survival curves", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(jump, Surv.ph.pns, type = "s", col = "blue")
lines(jump, Surv.ph.pcb, type = "s", col = "red", lty = 2)

plot(survfit(Surv(Time, death) ~ Trt, data = liver.unq), mark.time = FALSE,
     col=c("red", "blue"), lty = c(2, 1), ylab = "Survival", xlab = "Time (years)",lwd = 2,
     main = "Kaplan-Meier survival curves", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(jump, Surv.po.pns, type = "s", col = "blue")
lines(jump, Surv.po.pcb, type = "s", col = "red", lty = 2)