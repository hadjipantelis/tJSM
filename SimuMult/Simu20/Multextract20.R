
#=============== Extract the simulation results ===============#

Nsimu <- 100
phis <- alphas <- Ysigmas <- Bsigmas <- timeMaxs <- steps <- rep(0, Nsimu)
gammas <- matrix(0, Nsimu, 5)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 9)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("MultResult20", i, ".RData", sep = ""))
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
4.502   4.816   4.901   4.864   4.957   5.000 

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
[1]  3.0084370 -0.4615769 -1.6452152 -2.7964479 -3.5223477
apply(gammas, 2, sd)
[1] 0.1432680 0.1025350 0.1540303 0.2060878 0.1862048

mean(phis); sd(phis)
[1] 0.5054297
[1] 0.160251

mean(alphas); sd(alphas)
[1] -1.031992
[1] 0.171752

mean(Ysigmas); sd(Ysigmas)
[1] 0.497743
[1] 0.01743509

mean(Bsigmas); sd(Bsigmas)
[1] 0.649532
[1] 0.04485174

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.14161453 0.09858144 0.16811192 0.19493626 0.20064099 0.16455486 0.15379699 0.01757273 0.04557020

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.14634833 0.09886720 0.16939782 0.19752063 0.20480519 0.16459145 0.15425816 0.01801497 0.04719446

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.14806091 0.09897674 0.16984311 0.19815283 0.20563731 0.16459895 0.15405329 0.01846840 0.04837077


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 8), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
