
library(MASS)
library(nlme)

#===================== Function to Simulate Datasets For the Transformation Model II ====================#

SimuDataTM2 <- function (n, rho, beta, phi, alpha, Ysigma, BSigma, C0) {
  
  SurvTime <- rep(0, n)
  X1 <- rbinom(n, 1, 0.5) # Fixed covariate, Bernoulli r.v. with success prob. 0.5 #
  X2 <- runif(n, 0, 1) # Fixed covariate 2, uniform r.v. #
  CensorTime <- - log(runif(n)) * C0 + 0.20 # random sample from Exp(C0) #
  
  # ========== Generate Survival Time =========== #
  if (is.matrix(BSigma)) ncb <- ncol(BSigma) else ncb <- 1
  prob <- rep(0, n)
  k <- 1
  b <- matrix(0, n, ncb) # record the random effects #
  while (k <= n) {
    b[k, ] <- mvrnorm(1, rep(0, ncb), BSigma) # random effect #
    prob[k] <- runif(1)
    c1 <- phi[1] * X1[k] + phi[2] * X2[k] + alpha * b[k, 1]
    c2 <- alpha * b[k, 2]
    temp <- if(rho == 0) log(1 - prob[k]) else (1 - (1 - prob[k]) ^ (- rho)) / rho
    seed <- log(1 - c2 * temp / exp(c1)) / c2
    
    if (is.na(seed)) {
      print(c(X1[k], X2[k], c1, c2, temp, seed))
    }
    
    if (!is.na(seed)) {
      SurvTime[k] <- seed
      k <- k+1
    }              
  }
  
  # =========== Construct Censoring indicator and Observed Time ======== #
  ObservedTime <- SurvTime
  CensorIndex <- as.numeric(CensorTime > SurvTime)
  ObservedTime[CensorIndex == 0] <- CensorTime[CensorIndex == 0]
  
  # ========= Generate Longitudinal Measurements ========= #
  Rcov <- list() # random covariates #
  MeasureTime <- list()
  Ttemp <- 0.20 * ((1:1000) - 1)
  ni <- rep(0, n)
  for (i in 1:n) {
    MeasureTime[[i]] <- Ttemp[Ttemp <= ObservedTime[i]]
    ni[i] <- length(MeasureTime[[i]]) # the number of longitudinal observations #
    Rcov[[i]] <- beta[1] * X1[i] + beta[2] * X2[i] + beta[3] * MeasureTime[[i]] + beta[4] * X1[i] * MeasureTime[[i]] +
      beta[5] * X2[i] * MeasureTime[[i]] + b[i, 1] + b[i,2] * MeasureTime[[i]] + Ysigma * rnorm(ni[i])
    # N(0, Ysigma^2), measurement error #
  }
  
  # ========= Turn the Data into a Dataframe ========== #
  ID <- rep(1:n, ni)
  Time <- rep(ObservedTime, ni)
  death <- rep(CensorIndex, ni)
  Y <- unlist(Rcov)
  X1 <- rep(X1, ni)
  X2 <- rep(X2, ni)
  ObsTime <- unlist(MeasureTime)
  start <- unlist(MeasureTime)
  stop <- sapply(1:n, function(x) if(ni[x] >= 2) c(MeasureTime[[x]][2:ni[x]], ObservedTime[x]) 
                 else ObservedTime[x])
  stop <- unlist(stop)
  event <- sapply(1:n, function(x) if(ni[x] >= 2) c(rep(0,(ni[x]-1)), CensorIndex[x])
                  else CensorIndex[x])
  event <- unlist(event)
  
  result <- data.frame(ID = ID, Time = Time, death = death, Y = Y, X1 = X1, X2 = X2, ObsTime = ObsTime, 
                       start = start, stop = stop, event = event)
  return(result)
}


#=================================== TMData20: for Model II with rho = 0 ===================================#

beta <- c(-1, -1.5, 1, -0.5, 0.5); Ysigma <- sqrt(0.1)
phi <- c(-0.5, 1.5); alpha <- 0.5; C0 <- 2; n <- 200
BSigma <- matrix(c(0.5, -0.1, -0.1, 0.1), nrow = 2)
rho <- 0

Nsimu <- 100
TMData20 <- list()
set.seed(650921)
for (i in 1:Nsimu) {
  print(i)
  TMData20[[i]] <- SimuDataTM2(n, rho, beta, phi, alpha, Ysigma, BSigma, C0)
}
save(TMData20, file = "TMData20.RData")


#=============== Calculate Censoring Proportion ===============#
CenProp20 <- sapply(1:Nsimu, function(x) sum(TMData20[[x]]$event))
mean(CenProp20) * 100 / n # observed proportion is 81.72% #

Ni20 <- sapply(1:Nsimu, function(x) nrow(TMData20[[x]]))
mean(Ni20) / n # 3.23205 #

Ki20 <- sapply(1:Nsimu, function(x) {Temp <- TMData20[[x]]$Time[!duplicated(TMData20[[x]]$ID)] 
                                     return(length(Temp[Temp > 0.20]))})
mean(Ki20) * 100 / n # 68.955% #


#============== Test whether the data can be run by lme function ==============#
betas20 <- matrix(0, Nsimu, 5)
BSigmas20 <- matrix(0, nrow = Nsimu, ncol = 4)
Ysigmas20 <- rep(0, Nsimu)
for (i in 1:Nsimu) {
  print(i)
  Data <- TMData20[[i]]
  fitLME <- lme(Y ~ X1 + X2 + ObsTime + X1:ObsTime + X2:ObsTime - 1, random = ~ObsTime | ID, data = Data)
  betas20[i, ] <- as.vector(fixef(fitLME))
  BSigmas20[i, ] <- c(lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                             function(x) x*fitLME$sigma^2)[[1]])
  Ysigmas20[i] <- fitLME$sigma
}

colMeans(betas20)
[1] -0.9914439 -1.5024340  0.9572475 -0.4799075  0.4624599
apply(betas20, 2, sd)
[1] 0.09441617 0.11117556 0.10267751 0.10685368 0.21469223

colMeans(BSigmas20)
[1]  0.48694423 -0.10493611 -0.10493611  0.09609124
apply(BSigmas20, 2, sd)
[1] 0.05914280 0.04037592 0.04037592 0.03876701

mean(Ysigmas20); sd(Ysigmas20)
[1] 0.3170999
[1] 0.01210049




#=================================== TMData21: for Model II with rho = 1 ===================================#

beta <- c(-1, -1.5, 1, -0.5, 0.5); Ysigma <- sqrt(0.1)
phi <- c(-0.5, 1.5); alpha <- 0.5; C0 <- 2; n <- 200
BSigma <- matrix(c(0.5, -0.1, -0.1, 0.1), nrow = 2)
rho <- 1

Nsimu <- 100
TMData21 <- list()
set.seed(630726)
for (i in 1:Nsimu) {
  print(i)
  TMData21[[i]] <- SimuDataTM2(n, rho, beta, phi, alpha, Ysigma, BSigma, C0)
}
save(TMData21, file = "TMData21.RData")


#=============== Calculate Censoring Proportion ===============#
CenProp21 <- sapply(1:Nsimu, function(x) sum(TMData21[[x]]$event))
mean(CenProp21) * 100 / n # observed proportion is 69.82% #

Ni21 <- sapply(1:Nsimu, function(x) nrow(TMData21[[x]]))
mean(Ni21) / n # 4.4234 #

Ki21 <- sapply(1:Nsimu, function(x) {Temp <- TMData21[[x]]$Time[!duplicated(TMData21[[x]]$ID)] 
                                     return(length(Temp[Temp > 0.20]))})
mean(Ki21) * 100 / n # 73.05% #


#============== Test whether the data can be run by lme function ==============#
betas21 <- matrix(0, Nsimu, 5)
BSigmas21 <- matrix(0, nrow = Nsimu, ncol = 4)
Ysigmas21 <- rep(0, Nsimu)
for (i in 1:Nsimu) {
  print(i)
  Data <- TMData21[[i]]
  fitLME <- lme(Y ~ X1 + X2 + ObsTime + X1:ObsTime + X2:ObsTime - 1, random = ~ObsTime | ID, data = Data)
  betas21[i, ] <- as.vector(fixef(fitLME))
  BSigmas21[i, ] <- c(lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                             function(x) x*fitLME$sigma^2)[[1]])
  Ysigmas21[i] <- fitLME$sigma
}

colMeans(betas21)
[1] -1.0188476 -1.4970827  1.0276706 -0.4846753  0.4577615
apply(betas21, 2, sd)
[1] 0.10241660 0.10870653 0.08995010 0.08363135 0.15759618

colMeans(BSigmas21)
[1]  0.49566937 -0.09341272 -0.09341272  0.08808072
apply(BSigmas21, 2, sd)
[1] 0.06079638 0.03849692 0.03849692 0.02851020

mean(Ysigmas21); sd(Ysigmas21)
[1] 0.3162195
[1] 0.009152595
