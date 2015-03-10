
#=============== Extract the simulation results ===============#

Nsimu <- 100
betas <- matrix(0, Nsimu, 5)
phis <- matrix(0, Nsimu, 2)
alphas <- Ysigmas <- timeMaxs <- steps <- rep(0, Nsimu)
BSigmas <- matrix(0, Nsimu, 4)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 12)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("TMResult20", i, ".RData", sep = ""))
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
1.791   2.797   3.322   3.540   4.108   6.844

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
[1] -0.9916441 -1.4990298  0.9861832 -0.5000138  0.5153802
apply(betas, 2, sd)
1] 0.09411244 0.11048474 0.10216959 0.10762447 0.21471846

colMeans(phis)
[1] -0.5175303  1.5172244
apply(phis, 2, sd)
[1] 0.1778527 0.2998368

mean(alphas); sd(alphas)
[1] 0.5085623
[1] 0.1490432

mean(Ysigmas); sd(Ysigmas)
[1] 0.3168506
[1] 0.01203579

colMeans(BSigmas)
[1]  0.48252977 -0.09588877 -0.09588877  0.08546845
apply(BSigmas, 2, sd)
[1] 0.05848844 0.03945420 0.03945420 0.03510579

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.09372837 0.11599653 0.10850033 0.11051619 0.19867104 0.16929537 0.30026743 0.13676963 0.01132468 0.05527839 0.04304475
[12] 0.03440780

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.09373505 0.11605290 0.10845529 0.11056304 0.19864071 0.16930902 0.30044475 0.13687013 0.01137144 0.05536588 0.04330149
[12] 0.03476347

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.09372998 0.11599466 0.10858240 0.11058884 0.19872060 0.16926323 0.30004864 0.13675446 0.01141831 0.05544912 0.04350651
[12] 0.03511165


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 7), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
