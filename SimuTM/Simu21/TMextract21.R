
#=============== Extract the simulation results ===============#

Nsimu <- 100
betas <- matrix(0, Nsimu, 5)
phis <- matrix(0, Nsimu, 2)
alphas <- Ysigmas <- timeMaxs <- steps <- rep(0, Nsimu)
BSigmas <- matrix(0, Nsimu, 4)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 12)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("TMResult21", i, ".RData", sep = ""))
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
2.549   3.977   4.784   5.319   6.463  10.920 

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
[1] -1.0190251 -1.4925219  1.0429095 -0.4898691  0.4707946
apply(betas, 2, sd)
[1] 0.10270225 0.10885773 0.09080180 0.08256463 0.15719484

colMeans(phis)
[1] -0.4847724  1.4113613
apply(phis, 2, sd)
[1] 0.2587794 0.4135502

mean(alphas); sd(alphas)
[1] 0.5317907
[1] 0.1886242

mean(Ysigmas); sd(Ysigmas)
[1] 0.3162848
[1] 0.009139265

colMeans(BSigmas)
[1]  0.49076563 -0.08852344 -0.08852344  0.08067286
apply(BSigmas, 2, sd)
[1] 0.06004650 0.03750105 0.03750105 0.02657859

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.094333886 0.116631318 0.086146785 0.084881293 0.146695023 0.263168944 0.460949239 0.205898260 0.009031031 0.055232905
[11] 0.033743336 0.023628708

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.094340083 0.116738755 0.086095442 0.084895321 0.146659248 0.262924673 0.463439700 0.206723618 0.009067517 0.055301530
[11] 0.033898516 0.023889579

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.094332085 0.116611498 0.086155220 0.084892753 0.146668771 0.262694119 0.457662636 0.205555385 0.009103878 0.055370384
[11] 0.034039385 0.024149028


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 7), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)

# Actually, the cumulative hazard function estimate is not quite good under this case #
