
#========== Initial test of the code for the NMRE model ==========#

library(statmod)
library(mvtnorm)
library(nlme)
library(survival)
library(splines)

source('logLik.jmodelMult.R')
source('vcov.jmodelMult.R')
source('AIC.jmodelMult.R')
source('BIC.jmodelMult.R')
source('jmodelMult.R')
source('InitValMultGeneric.R') 
source('EMiterMultGeneric.R')Generic
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
control <- list(tol.P = 10 ^ (- 4), tol.L = 10 ^ (- 8), max.iter = 100, nknot = 5)


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

library(lattice)
load("pbc.rda")
pbc.unq <- pbc[!duplicated(pbc$ID), ]
dim(pbc.id[pbc.id$drug == "placebo", ])
[1] 154  16
dim(pbc.id[pbc.id$drug != "placebo", ])
[1] 158  16
xyplot(log(serBilir) ~ obstime | drug, group = ID, data = pbc,
       xlab = "Years", ylab = "log(serum bilirubin)", col = "darkgrey", type = "l")
pbc.Dpca <- pbc[pbc$drug == "D-penicil", ]
pbc.pcb <- pbc[pbc$drug == "placebo", ]
times.Dpca <- sort(unique(pbc.Dpca$obstime))
times.pcb <- sort(unique(pbc.pcb$obstime))
means.Dpca <- tapply(log(pbc.Dpca$serBilir), pbc.Dpca$obstime, mean)
means.pcb <- tapply(log(pbc.pcb$serBilir), pbc.pcb$obstime, mean)
plot(lowess(times.Dpca, means.Dpca, f = .5), type = "l", col = "blue", lwd = 2,
     ylab = "Prothrombin index", xlab = "Time (years)",
     main = "Smoothed mean prothrombin index", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
lines(lowess(times.pcb, means.pcb, f = .5), col = "red", lty = 2, lwd = 2)

plot(survfit(Surv(Time, death) ~ drug, data = pbc.unq), mark.time = FALSE, 
     col=c("red", "blue"), lty = c(2, 1), ylab = "Survival", xlab = "Time (years)", lwd = 2,
     main = "Kaplan-Meier survival curves", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

fitLME <- lme(log(serBilir) ~ bs(obstime, df = 6, degree = 2), random = ~ 1 | ID, data = pbc)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = pbc, x = TRUE)
fitJT.ph <- jmodelMult(fitLME, fitCOX, pbc, control = list(tol.P = 1e-3))
fitLME2 <- lme(log(serBilir) ~ bs(obstime, df = 4, degree = 2), random = ~ 1 | ID, data = pbc)
fitJT.ph2 <- jmodelMult(fitLME2, fitCOX, pbc, control = list(tol.P = 1e-3))
fitLME3 <- lme(log(serBilir) ~ bs(obstime, df = 2, degree = 2), random = ~ 1 | ID, data = pbc)
fitJT.ph3 <- jmodelMult(fitLME3, fitCOX, pbc, control = list(max.iter = 300, tol.P = 1e-3))
fitLME4 <- lme(log(serBilir) ~ bs(obstime, df = 3, degree = 1), random = ~ 1 | ID, data = pbc)
fitJT.ph4 <- jmodelMult(fitLME4, fitCOX, pbc, control = list(max.iter = 300, tol.P = 1e-3))

load("pbcResults.RData")
fit <- pbcResults$fitJT.ph3
bi <- fit$est.bi
gamma <- fit$coefficients$gamma
times <- sort(pbc$obstime)
# set.seed(12321)
# patients <- sort(sample(1:312, 6))
patients  <- c(24, 73, 138, 186, 218, 261)
B <- bs(times, df = 2, degree = 2)
mu <- gamma[1] + gamma[2] * B[, 1] + gamma[3] * B[, 2]
# plot(times, mu, type = "l", ylim = c(- 1, 4))
par(mfrow = c(2, 3))
for (i in 1:6) {
  data <- pbc[pbc$ID == patients[i], ]
  mui <- bi[patients[i]] * mu
  # ylim <- range(c(mui, log(data$serBilir))) + c(- 1, 1)
  plot(times, mui, type = "l", xlab = "Time (years)", ylab = "log(serBilir)", ylim = c(-2, 8), cex.main = 1.5,
       cex.lab = 1.5, cex.axis = 1.5, lwd = 2, col = "chocolate", main = paste("Patient", patients[i], seq = " "))
  points(data$obstime, log(data$serBilir), pch = 19, col = "darkgreen")
}

plot(survfit(Surv(Time, death) ~ drug, data = pbc.unq), mark.time = FALSE, fun = "cumhaz",
     col=c("red", "blue"), lty = c(2, 1), ylab = "Survival", xlab = "Time (years)", lwd = 2,
     main = "Cumulative hazard functions", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# Lamb <- cumsum(fit$coefficients$lamb[, 2])
lamb <- fit$coefficients$lamb[, 2]
jump <- fit$coefficients$lamb[, 1]
# plot(jump, Lamb, type = "s")
phi <- fit$coefficients$phi
alpha <- fit$coefficients$alpha
Bjump <- bs(jump, df = 2, degree = 2)
mujump <- gamma[1] + gamma[2] * Bjump[, 1] + gamma[3] * Bjump[, 2]
Lamb.Dpca <- cumsum(lamb * exp(phi + alpha * mujump))
Lamb.pcb <- cumsum(lamb * exp(alpha * mujump))
lines(jump, Lamb.Dpca, type = "s", col = "blue")
lines(jump, Lamb.pcb, type = "s", col = "red", lty = 2)



fitJT.po <- jmodelMult(fitLME, fitCOX, pbc, rho = 1)
fitJT.po2 <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho = 1)

fitCOX2 <- coxph(Surv(start, stop, event) ~ drug + as.numeric(drug):obstime, data = pbc, x = TRUE)
fitJT22 <- jmodelMult(fitLME, fitCOX2, pbc, timeVarT = 'obstime')
