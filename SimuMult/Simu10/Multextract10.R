
#=============== Extract the simulation results ===============#

Nsimu <- 100
phis <- alphas <- Ysigmas <- Bsigmas <- timeMaxs <- steps <- rep(0, Nsimu)
gammas <- matrix(0, Nsimu, 5)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 9)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("MultResult10", i, ".RData", sep = ""))
  gammas[i, ] <- results$theta.hat$gamma
  phis[i] <- results$theta.hat$phi
  alphas[i] <- results$theta.hat$alpha
  Ysigmas[i] <- results$theta.hat$Ysigma
  Bsigmas[i] <- results$theta.hat$Bsigma
  Vars1[i, ] <- results$Var1
  Vars2[i, ] <- results$Var2
  Vars3[i, ] <- results$Var3
  
  Lambs[[i]] <- cumsum(results$theta.hat$lamb[, 2])
  times[[i]] <- results$theta.hat$lamb[, 1]
  timeMaxs[i] <- max(times[[i]])
  steps[i] <- results$step
}

summary(timeMaxs)
Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
4.818   4.935   4.967   4.958   4.987   5.000 

pooltimes <- seq(from = 0, to = 5, length.out = 5000)
len <- length(pooltimes)
EstLamb <- matrix(0, Nsimu, len)
for (i in 1:Nsimu) {
  Time <- times[[i]]
  Lamb <- Lambs[[i]]
  Ni <- length(Time)
  EstLamb[i, ] <- sapply(pooltimes, function(x) {Ind <- (1 : Ni)[Time <= x]; if(length(Ind) == 0) return(0)
                                                 else return(Lamb[max(Ind)])})
}


colMeans(gammas)
[1]  3.0253253 -0.4780988 -1.6341274 -2.8343538 -3.5218093
apply(gammas, 2, sd)
[1] 0.13149810 0.09342734 0.17192764 0.18924858 0.21318600

mean(phis); sd(phis)
[1] 0.002465888
[1] 0.1682731

mean(alphas); sd(alphas)
[1] -1.027267
[1] 0.1177385

mean(Ysigmas); sd(Ysigmas)
[1] 0.5002198
[1] 0.01318249

mean(Bsigmas); sd(Bsigmas)
[1] 0.6410368
[1] 0.04239053

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.14057731 0.09080120 0.15035244 0.18612462 0.21190653 0.16328881 0.12178699 0.01453001 0.04462668

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.15358748 0.09169891 0.15458785 0.19445588 0.22640361 0.16348526 0.12151612 0.01456705 0.04678146

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.13813268 0.09049986 0.14939771 0.18465976 0.20813993 0.16218281 0.12120590 0.01460229 0.04435795


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 6), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
