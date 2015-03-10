
#=============== Joint Modeling Main Function with NMRE (nonparametric Multiplicative random effects) ===============#

jmodelMult <- function (fitLME, fitCOX, data, model = 1, rho = 0, timeVarT = NULL, 
                        control = list(), ...) 
{
  call <- match.call()
  if(!inherits(fitLME, "lme"))
    stop("\n'fitLME'must be a fit returned by lme().")
  if(length(fitLME$group) > 1)
    stop("\n nested random-effects are not allowed in lme().")
  if(!is.null(fitLME$modelStruct$corStruct))
    warning("\n correlation structure in 'fitLME' is ignored.")
  if(!is.null(fitLME$modelStruct$varStruct))
    warning("\n heteroscedasticity structure in 'fitLME' is ignored.")
  
  if(!inherits(fitCOX, "coxph"))
    stop("\n'fitCOX' must be a fit returned by coxph().")
  if(is.null(fitCOX$x))
    stop("\n must specify argument 'x=TRUE' when using coxph().")
  
  if (rho < 0) {
    rho <- 0 # fit Cox model if not specified #
    warning("\n rho<0 is not valid, Cox model is fitted instead!")
  }
  
  ID <- as.vector(unclass(fitLME$groups[[1]])) 
  ni <- as.vector(tapply(ID, ID, length))           
  nLong <- length(ni)
  
  if(ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")
  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  Time <- stop[cumsum(ni)]
  d<- event[cumsum(ni)] 
  nSurv <- length(Time) 
  if(sum(d) < 5)
    warning("\n more than 5 events are required.")
  if(nLong != nSurv)
    stop("\n sample sizes in the longitudinal and event processes differ.")
  
  Z <- as.matrix(fitCOX$x)
  phi.names <- colnames(Z)
  formSurv <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, data)[cumsum(ni), ]
  if(!is.null(timeVarT)) {
    if(!all(timeVarT %in% names(mfSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  Ztime <- as.matrix(model.matrix(formSurv, mfSurv))
  if(attr(TermsSurv, 'intercept')) Ztime <- as.matrix(Ztime[, - 1])
  # design matrix of the covariates in survival part, one row for each subject, excluding intercept #
  
  TermsLongX <- fitLME$terms
  mydata <- fitLME$data[all.vars(TermsLongX)] 
  formLongX <- formula(fitLME) 
  mfLongX <- model.frame(TermsLongX, data = mydata) 
  B <- as.matrix(model.matrix(formLongX, mfLongX)) # may include a column of intercept #
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed, of length N #
  alpha.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]
  
  tempForm <- strsplit(toString(splitFormula(formLongX)[[1]]), ",")[[1]]
  tempForm <- tempForm[- 1]
  tempForm[1] <- " bs(Time"
  tempForm <- paste(tempForm, collapse = ",")
  Btime <- as.matrix(eval(parse(text = tempForm)))
  if(attr(TermsLongX, 'intercept')) Btime <- cbind(rep(1, nLong), Btime)
  
  U <- sort(unique(Time[d == 1])) # ordered uncensored observed event time #
  tempU <- lapply(Time, function(t) U[t >= U]) 
  times <- unlist(tempU) # vector of length M #
  nk <- sapply(tempU, length)  # length of each element in times, vector of length n #
  M <- sum(nk)
  Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Index0 <- match(Time, U)
  # vector of length n, for the uncensored subjects, return the value as nk, otherwise return NA # 
  Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu #
  
  tempForm2 <- strsplit(toString(splitFormula(formLongX)[[1]]), ",")[[1]]
  tempForm2 <- tempForm2[- 1]
  tempForm2[1] <- " bs(times"
  tempForm2 <- paste(tempForm2, collapse = ",")
  Btime2 <- as.matrix(eval(parse(text = tempForm2)))
  if(attr(TermsLongX, 'intercept')) Btime2 <- cbind(rep(1, M), Btime2)
  
  mfSurv2 <- mfSurv[Index, ]
  if(!is.null(timeVarT)) {
    mfSurv2[timeVarT] <- times
  }
  Ztime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
  if(attr(TermsSurv, 'intercept')) Ztime2 <- as.matrix(Ztime2[, - 1]) # excluding intercept #
  
  n <- nLong
  N <- length(Y)
  nu <- length(U)
  ncz <- ncol(Z)
  ncb <- ncol(B)
  
  controlvals <- list(tol.P = 10 ^ (-4), tol.L = 10 ^ (-8), max.iter = 200, SE.method = 'PRES', delta = 10 ^ (- 5), 
                      nknot = 12)
  control <- c(control, list(...))
  namec <- names(control)
  if (length(uname <- namec[!namec %in% names(controlvals)]) > 0) 
    warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
  controlvals[namec] <- control
  if (controlvals$SE.method == 'PLFD' | controlvals$SE.method == 'PFDS') controlvals$delta <- 10 ^ (- 3)
  controlvals[namec] <- control
  
  tol.P <- controlvals$tol.P
  tol.L <- controlvals$tol.L
  iter <- controlvals$max.iter
  nknot <- controlvals$nknot
  
  GHQ <- gauss.quad(nknot, kind = "hermite")
  b <- GHQ$nodes
  wGQ <- GHQ$weights
  
  Y.st <- split(Y, ID)
  B.st <- lapply(split(B, ID), function(x) matrix(x, ncol = ncb))
  Ztime22 <- if(ncz > 1) t(apply(Ztime2, 1, function(x) x %o% x)) else Ztime2 ^ 2
  Btime22 <- if(ncb > 1) t(apply(Btime2, 1, function(x) x %o% x)) else Btime2 ^ 2
  B2 <- if(ncb > 1) t(apply(B, 1, function(x) x %o% x)) else B ^ 2
  
  if (model == 1) {
    environment(InitValMult1) <- environment(EMiterMult1) <- environment()
  } else {
    environment(InitValMult2) <- environment(EMiterMult2) <- environment()
  }
  
  tempResp <- strsplit(toString(formLongX), ", ")[[1]][c(2, 1)]
  tempResp <- paste(tempResp, collapse = "")
  tempForm3 <- strsplit(toString(splitFormula(formLongX)[[1]]), ", ")[[1]]
  tempForm3 <- tempForm3[-1]
  tempForm3[length(tempForm3) + 1] <- "data = fitLME$data"
  tempForm3 <- paste(tempForm3, collapse = ",")
  tempForm3 <- paste("lm(", tempResp, tempForm3, ")", sep = "")
  fitLM <- eval(parse(text = tempForm3))
  gamma <- as.vector(fitLM$coefficients)
  
  surv.init <- if (model == 1) InitValMult1(gamma) else InitValMult2(gamma)
  phi <- surv.init$phi
  alpha <- surv.init$alpha
  lamb <- surv.init$lamb
  Ysigma <- surv.init$Ysigma
  Bsigma <- surv.init$Bsigma
  
  theta.old <- list(gamma = gamma, phi = phi,  alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma, 
                   lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1
  
  while (step <= iter) {
    
    if (err.P < tol.P | err.L < tol.L) break
    
    theta.new <- if (model == 1) EMiterMult1(theta.old) else EMiterMult2(theta.old)
    
    new.P <- c(theta.new$gamma, theta.new$phi, theta.new$alpha, theta.new$Ysigma, theta.new$Bsigma)
    old.P <- c(theta.old$gamma, theta.old$phi, theta.old$alpha, theta.old$Ysigma, theta.old$Bsigma)
    err.P <- max(abs(new.P - old.P) / (abs(old.P) + tol.P))
    # add tol.P to avoid zero value of the estimated parameters #
    
    new.L <- theta.new$lgLik
    old.L <- theta.old$lgLik
    err.L <- abs(new.L - old.L) / (abs(old.L) + tol.P)
    
    step <- step + 1
    theta.old <- theta.new
  }
  converge <- as.numeric(err.P < tol.P | err.L < tol.L)
  
  delta <- controlvals$delta
  environment(SfuncMult) <- environment()  
  if (model == 1) {
    environment(LambMult1) <- environment(DQfuncMult1) <- environment(LHMult1) <- environment()
  } else {
    environment(LambMult2) <- environment(DQfuncMult2) <- environment(LHMult2) <- environment()
  }
  if (controlvals$SE.method == 'PFDS') {
    environment(PFDSMult) <- environment()
    time.SE <- system.time(Vcov <- PFDSMult(model, theta.new, min(tol.P, delta) / 100, iter, delta))[3]
    if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
      warning("NA's present in StdErr estimation due to numerical error!\n")
  } else if (controlvals$SE.method == 'PRES') {
    environment(PRESMult) <- environment()
    if (CheckDeltaMult(theta.new, delta)) {
      time.SE <- system.time(Vcov <- PRESMult(model, theta.new, min(tol.P, delta) / 100, iter, delta))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (controlvals$SE.method == 'PLFD') {
    environment(PLFDMult) <- environment()
    time.SE <- system.time(Vcov <- PLFDMult(model, theta.new, min(tol.P, delta) / 100, iter, delta))[3]
    if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
      warning("NA's present in StdErr estimation due to numerical error!\n")
  } else {
    Vcov <- time.SE <- NA
    warning("\n Standard error estimation method should be either 'PFDS', 'PRES' or 'PLFD'.")
  }
  
  theta.new$lamb <- cbind("time" = U, "bashaz" = theta.new$lamb)
  names(theta.new$gamma) <- paste("gamma.", 1:ncb, sep = "")
  names(theta.new$phi) <- phi.names
  names(theta.new$alpha) <- if (model == 1) alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"
  names(theta.new$Bsigma) <- "sigma.b"
  
  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if (converge == 1) "success" else "failure"
  result$control <- controlvals
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  class(result) <- "jmodelMult"
  result
}
