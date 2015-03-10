
library(MASS)
library(nlme)

#===================== Function to Simulate Datasets For the Transformation Model I ====================#

SimuDataTM1 <- function (n, rho, beta, phi, alpha, Ysigma, BSigma, C0) {
  
  SurvTime <- rep(0, n)
  X1 <- rbinom(n, 1, 0.5) # Fixed covariate, Bernoulli r.v. with success prob. 0.5 #
  X2 <- runif(n, 0, 1) # Fixed covariate 2, uniform r.v. #
  CensorTime <- - log(runif(n)) * C0 + 0.25 # random sample from Exp(C0) #
  
  # ========== Generate Survival Time =========== #
  if (is.matrix(BSigma)) ncb <- ncol(BSigma) else ncb <- 1
  prob <- rep(0, n)
  k <- 1
  b <- matrix(0, n, ncb) # record the random effects #
  while (k <= n) {
    b[k, ] <- mvrnorm(1, rep(0, ncb), BSigma) # random effect #
    prob[k] <- runif(1)
    c1 <- phi[1] * X1[k] + phi[2] * X2[k] + alpha * beta[1] * X1[k] + alpha * beta[2] * X2[k] + alpha * b[k, 1]
    c2 <- alpha * (beta[3] + beta[4] * X1[k] + beta[5] * X2[k] + b[k, 2])
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
  Ttemp <- 0.25 * ((1:1000) - 1)
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


#=================================== TMData10: for Model I with rho = 0 ===================================#

beta <- c(-1, -1.5, 1, -0.5, 0.5); Ysigma <- sqrt(0.1)
phi <- c(-0.5, 1.5); alpha <- 0.5; C0 <- 2; n <- 200
BSigma <- matrix(c(0.5, -0.1, -0.1, 0.1), nrow = 2)
rho <- 0

Nsimu <- 100
TMData10 <- list()
set.seed(12321)
for (i in 1:Nsimu) {
  print(i)
  TMData10[[i]] <- SimuDataTM1(n, rho, beta, phi, alpha, Ysigma, BSigma, C0)
}
save(TMData10, file = "TMData10.RData")


#=============== Calculate Censoring Proportion ===============#
CenProp10 <- sapply(1:Nsimu, function(x) sum(TMData10[[x]]$event))
mean(CenProp10) * 100 / n # observed proportion is 76.92% #

Ni10 <- sapply(1:Nsimu, function(x) nrow(TMData10[[x]]))
mean(Ni10) / n # 3.3086 #

Ki10 <- sapply(1:Nsimu, function(x) {Temp <- TMData10[[x]]$Time[!duplicated(TMData10[[x]]$ID)] 
                                     return(length(Temp[Temp > 0.25]))})
mean(Ki10) * 100 / n # 75.69% #


#============== Test whether the data can be run by lme function ==============#
betas10 <- matrix(0, Nsimu, 5)
BSigmas10 <- matrix(0, nrow = Nsimu, ncol = 4)
Ysigmas10 <- rep(0, Nsimu)
for (i in 1:Nsimu) {
  print(i)
  Data <- TMData10[[i]]
  fitLME <- lme(Y ~ X1 + X2 + ObsTime + X1:ObsTime + X2:ObsTime - 1, random = ~ObsTime | ID, data = Data)
  betas10[i, ] <- as.vector(fixef(fitLME))
  BSigmas10[i, ] <- c(lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                             function(x) x*fitLME$sigma^2)[[1]])
  Ysigmas10[i] <- fitLME$sigma
}

colMeans(betas10)
[1] -1.0053577 -1.4935080  0.9783125 -0.4802163  0.4643478
apply(betas10, 2, sd)
[1] 0.08753473 0.10925874 0.10840754 0.10901864 0.14596493

colMeans(BSigmas10)
[1]  0.5022298 -0.1103429 -0.1103429  0.1084500
apply(BSigmas10, 2, sd)
[1] 0.06102760 0.04601910 0.04601910 0.03751946

mean(Ysigmas10); sd(Ysigmas10)
[1] 0.3176541
[1] 0.01154051




#=================================== TMData11: for Model I with rho = 1 ===================================#

beta <- c(-1, -1.5, 1, -0.5, 0.5); Ysigma <- sqrt(0.1)
phi <- c(-0.5, 1.5); alpha <- 0.5; C0 <- 2; n <- 200
BSigma <- matrix(c(0.5, -0.1, -0.1, 0.1), nrow = 2)
rho <- 1

Nsimu <- 100
TMData11 <- list()
set.seed(860416)
for (i in 1:Nsimu) {
  print(i)
  TMData11[[i]] <- SimuDataTM1(n, rho, beta, phi, alpha, Ysigma, BSigma, C0)
}
save(TMData11, file = "TMData11.RData")


#=============== Calculate Censoring Proportion ===============#
CenProp11 <- sapply(1:Nsimu, function(x) sum(TMData11[[x]]$event))
mean(CenProp11) * 100 / n # observed proportion is 65.04% #

Ni11 <- sapply(1:Nsimu, function(x) nrow(TMData11[[x]]))
mean(Ni11) / n # 4.2691 #

Ki11 <- sapply(1:Nsimu, function(x) {Temp <- TMData11[[x]]$Time[!duplicated(TMData11[[x]]$ID)] 
                                     return(length(Temp[Temp > 0.25]))})
mean(Ki11) * 100 / n # 79.135% #


#============== Test whether the data can be run by lme function ==============#
betas11 <- matrix(0, Nsimu, 5)
BSigmas11 <- matrix(0, nrow = Nsimu, ncol = 4)
Ysigmas11 <- rep(0, Nsimu)
for (i in 1:Nsimu) {
  print(i)
  Data <- TMData11[[i]]
  fitLME <- lme(Y ~ X1 + X2 + ObsTime + X1:ObsTime + X2:ObsTime - 1, random = ~ObsTime | ID, data = Data)
  betas11[i, ] <- as.vector(fixef(fitLME))
  BSigmas11[i, ] <- c(lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                             function(x) x*fitLME$sigma^2)[[1]])
  Ysigmas11[i] <- fitLME$sigma
}

colMeans(betas11)
[1] -1.0140828 -1.5028921  0.9775939 -0.4963007  0.5019190
apply(betas11, 2, sd)
[1] 0.08645914 0.12944817 0.08730133 0.08966833 0.13761006

colMeans(BSigmas11)
[1]  0.50126564 -0.10083440 -0.10083440  0.09865637
apply(BSigmas11, 2, sd)
[1] 0.05700463 0.03259590 0.03259590 0.02458211

mean(Ysigmas11); sd(Ysigmas11)
[1] 0.3154579
[1] 0.009307832

