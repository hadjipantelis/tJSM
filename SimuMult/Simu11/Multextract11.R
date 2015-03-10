
#=============== Extract the simulation results ===============#

Nsimu <- 100
phis <- alphas <- Ysigmas <- Bsigmas <- timeMaxs <- steps <- rep(0, Nsimu)
gammas <- matrix(0, Nsimu, 5)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 9)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("MultResult11", i, ".RData", sep = ""))
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
6.593   6.865   6.926   6.893   6.968   6.997

pooltimes <- seq(from = 0, to = 7, length.out = 5000)
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
[1]  3.004164 -0.458974 -2.126240 -3.723975 -4.923485
apply(gammas, 2, sd)
[1] 0.06238178 0.09433020 0.17034236 0.21845182 0.20982069

mean(phis); sd(phis)
[1] 0.02536715
[1] 0.2644459

mean(alphas); sd(alphas)
[1] -0.9919962
[1] 0.3131931

mean(Ysigmas); sd(Ysigmas)
[1] 0.4962618
[1] 0.0152755

mean(Bsigmas); sd(Bsigmas)
[1] 0.249276
[1] 0.0170767

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.06241422 0.09674504 0.17426888 0.21623555 0.22961613 0.25723585 0.33986073 0.01472396 0.01649666

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.06288365 0.09688195 0.17441113 0.21631946 0.22969076 0.25721872 0.33954194 0.01508856 0.01709611

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.06335052 0.09701854 0.17455194 0.21640404 0.22976060 0.25721785 0.33939911 0.01546049 0.01771994


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(-1, 11), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
