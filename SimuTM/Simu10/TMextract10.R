
#=============== Extract the simulation results ===============#

Nsimu <- 100
betas <- matrix(0, Nsimu, 5)
phis <- matrix(0, Nsimu, 2)
alphas <- Ysigmas <- timeMaxs <- steps <- rep(0, Nsimu)
BSigmas <- matrix(0, Nsimu, 4)
Vars1 <- Vars2 <- Vars3 <- matrix(0, Nsimu, 12)
Lambs <- times <- list()

for (i in 1:Nsimu) {
  load(paste("TMResult10", i, ".RData", sep = ""))
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
2.269   3.101   3.493   3.585   3.927   8.395

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
[1] -1.0069711 -1.4925265  1.0212041 -0.5118970  0.4885787
apply(betas, 2, sd)
[1] 0.08739826 0.10936002 0.10941022 0.10943929 0.14728831

colMeans(phis)
[1] -0.4414383  1.5578296
apply(phis, 2, sd)
[1] 0.2148380 0.3379356

mean(alphas); sd(alphas)
[1] 0.5228767
[1] 0.1337732

mean(Ysigmas); sd(Ysigmas)
[1] 0.3174513
[1] 0.01150793

colMeans(BSigmas)
[1]  0.4974456 -0.1025253 -0.1025253  0.0995590
apply(BSigmas, 2, sd)
[1] 0.06027939 0.04612795 0.04612795 0.03540319

all(Vars1 > 0)
[1] TRUE
colMeans(sqrt(Vars1))
[1] 0.09483863 0.11784986 0.11233875 0.10947511 0.16972399 0.23837434 0.34829509 0.13780801 0.01139517 0.05676161
[11] 0.03983865 0.03301571

all(Vars2 > 0)
[1] TRUE
colMeans(sqrt(Vars2))
[1] 0.09484422 0.11791959 0.11236653 0.10948147 0.16968982 0.23877390 0.34887139 0.13821388 0.01144134 0.05683767
[11] 0.03999425 0.03328023

all(Vars3 > 0)
[1] TRUE
colMeans(sqrt(Vars3))
[1] 0.09483964 0.11784863 0.11237023 0.10948630 0.16972776 0.23767315 0.34708103 0.13734543 0.01148727 0.05691545
[11] 0.04013100 0.03354228


AverageLamb <- colMeans(EstLamb)
SELamb <- apply(EstLamb, 2, sd)

plot(pooltimes, AverageLamb, type = "l", ylim = c(0, 7), main = "Cumulative Hazard Function", 
     xlab = "Time", ylab = "Cumulative hazard", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
abline(a = 0, b = 1, lty = 2, col = 'blue')
upper <- AverageLamb + 1.96 * SELamb
lower <- AverageLamb - 1.96 * SELamb
points(pooltimes, upper, type = "l", lty = 4)
points(pooltimes, lower, type = "l", lty = 4)
