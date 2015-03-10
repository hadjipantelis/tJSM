
#=============== Extract the simulation results ===============#

Nsimu <- 100
betas <- matrix(0, Nsimu, 5)
phis <- matrix(0, Nsimu, 2)
alphas <- Ysigmas <- timeMaxs <- steps <- rep(0, Nsimu)
BSigmas <- matrix(0, Nsimu, 4)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 12)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("TMResult11", i, ".RData", sep = ""))
  betas[i, ] <- results$theta.hat$beta
  phis[i, ] <- results$theta.hat$phi
  alphas[i] <- results$theta.hat$alpha
  Ysigmas[i] <- results$theta.hat$Ysigma
  BSigmas[i, ] <- c(results$theta.hat$BSigma)
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
2.949   4.085   4.843   5.134   5.781  12.610

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


colMeans(betas)
[1] -1.0151254 -1.5007982  0.9959338 -0.5068154  0.5094712
apply(betas, 2, sd)
[1] 0.08668025 0.12957422 0.09013072 0.09062991 0.13679709

colMeans(phis)
[1] -0.4946654  1.5053229
apply(phis, 2, sd)
[1] 0.3465345 0.5656073

mean(alphas); sd(alphas)
[1] 0.4922058
[1] 0.1937456

mean(Ysigmas); sd(Ysigmas)
[1] 0.315384
[1] 0.009313617

colMeans(BSigmas)
[1]  0.49642739 -0.09740614 -0.09740614  0.09327818
apply(BSigmas, 2, sd)
[1] 0.05645623 0.03215604 0.03215604 0.02376459

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.094476992 0.116907985 0.088533861 0.084040914 0.138087296 0.360975608 0.527169730 0.198801765 0.009356205 0.055837729
[11] 0.032503354 0.023767278

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.094479217 0.116994535 0.088592166 0.084053679 0.138093738 0.363327426 0.538404864 0.203823125 0.009393864 0.055903451
[11] 0.032633876 0.023987521

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.094476859 0.116896570 0.088547263 0.084046208 0.138067110 0.353776819 0.513489283 0.194179474 0.009431055 0.055967177
[11] 0.032747401 0.024202173


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(-1, 10), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)

# Dataset 42 provides an extreme estimate #
