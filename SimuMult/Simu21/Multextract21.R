
#=============== Extract the simulation results ===============#

Nsimu <- 100
phis <- alphas <- Ysigmas <- Bsigmas <- timeMaxs <- steps <- rep(0, Nsimu)
gammas <- matrix(0, Nsimu, 5)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 9)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("MultResult21", i, ".RData", sep = ""))
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
8.587   9.411   9.709   9.628   9.895   9.993  

pooltimes <- seq(from = 0, to = 10, length.out = 5000)
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
[1]  2.9918251 -0.4572239 -2.8214974 -5.0847178 -7.0110318
apply(gammas, 2, sd)
[1] 0.15372991 0.09048032 0.21356740 0.29069464 0.39637065

mean(phis); sd(phis)
[1] 0.02656205
[1] 0.2498397

mean(alphas); sd(alphas)
[1] -0.9435104
[1] 0.2374088

mean(Ysigmas); sd(Ysigmas)
[1] 0.4991247
[1] 0.01595526

mean(Bsigmas); sd(Bsigmas)
[1] 0.6546695
[1] 0.04914546

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.14199568 0.08831818 0.20177794 0.29277102 0.35319053 0.25354918 0.21515602 0.01700819 0.04615686

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.15541398 0.08900613 0.20990747 0.31212792 0.38296003 0.25356373 0.21514364 0.01745344 0.04927876

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.16595764 0.08957804 0.21632363 0.32778698 0.40604048 0.25344670 0.21616668 0.01796036 0.05211467


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 15), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
