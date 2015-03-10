
#========== Function to print basic summary of joint model ==========#

print.jmodelTM <- function (result, digits = max(4, getOption("digits") - 4), ...) 
{
  cat("Joint modeling fitting of longitudinal and survival data by ML\n")
    
  cat("\nCoefficients:\n")
  cat("\tLongitudinal Process\n")
  print(round(result$coefficients$beta, digits))
  cat("\tSurvival Process\n")
  phis <- c(result$coefficients$phi, result$coefficients$alpha)
  print(round(phis, digits))
  
  cat("\nVariance Components:\n")
  BSigma <- result$coefficients$BSigma
  ncz <- if(length(BSigma) == 1) 1 else nrow(BSigma)
  SD <- if (ncz == 1) sqrt(BSigma) else sqrt(diag(BSigma))
  SD <- c(SD, result$coefficients$Ysigma)
  if (ncz > 1) {
    corr <- cov2cor(BSigma)
    corr[upper.tri(corr, TRUE)] <- 0
    corr <- rbind(corr, rep(0, ncz))
    tempMat <- round(cbind(SD, corr[ , - ncz]), digits)
    tempMat <- apply(tempMat, 2, sprintf, fmt = "% .4f")
    tempMat[tempMat == tempMat[1, 2]] <- ""
    tempMat[1, -1] <- abbreviate(colnames(tempMat)[- 1], 6)
    colnames(tempMat) <- c(colnames(tempMat)[1], rep("", ncz - 1))
    Mat <- data.frame(tempMat, check.rows = FALSE, check.names = FALSE)
    colnames(Mat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
    rownames(Mat) <- c(dimnames(BSigma)[[1]], "Residual")
  } else {
    Mat <- data.frame("StdDev" = SD, row.names = c(names(BSigma), "Residual"), 
                      check.rows = FALSE, check.names = FALSE)
  }
  print(if(!is.numeric(Mat)) Mat else round(Mat, digits))
  
  cat("\nLog-likelihood:", result$logLik, "\n")
  cat("Number of Observations:", result$N, "\n")
  cat("Number of Subjects:", result$n, "\n")
  invisible(result)
}
