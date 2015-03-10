

jmodelTM <- function (fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, timeVarT = NULL, 
                      control = list(), ...)
{
  call <- match.call()
  if (!inherits(fitLME, "lme"))
    stop("\n'fitLME'must be a fit returned by lme().")
  if (length(fitLME$group) > 1)
    stop("\n nested random-effects are not allowed in lme().")
  if (!is.null(fitLME$modelStruct$corStruct))
    warning("\n correlation structure in 'fitLME' is ignored.")
  if (!is.null(fitLME$modelStruct$varStruct))
    warning("\n heteroscedasticity structure in 'fitLME' is ignored.")
  
  if (!inherits(fitCOX, "coxph"))
    stop("\n'fitCOX' must be a fit returned by coxph().")
  if (is.null(fitCOX$x))
    stop("\n must specify argument 'x=TRUE' when using coxph().")

  if (rho < 0) {
    rho <- 0 # fit Cox model if not specified #
    warning("\n rho<0 is not valid, Cox model is fitted instead!")
  }
  
  ID <- as.vector(unclass(fitLME$groups[[1]])) 
  ni <- as.vector(tapply(ID, ID, length))           
  bBLUP <- data.matrix(ranef(fitLME)) 
  dimnames(bBLUP) <- NULL
  nLong <- nrow(bBLUP)
  
  if (ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")
  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  Time <- stop[cumsum(ni)]
  d <- event[cumsum(ni)] 
  nSurv <- length(Time) 
  if (sum(d) < 5)
    warning("\n more than 5 events are required.")
  if (nLong != nSurv)
    stop("\n sample sizes in the longitudinal and event processes differ.")
  
  W <- as.matrix(fitCOX$x)
  phi.names <- colnames(W)
  formSurv <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, data)[cumsum(ni), ]
  if (!is.null(timeVarT)) {
    if (!all(timeVarT %in% all.vars(TermsSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  Wtime <- as.matrix(model.matrix(formSurv, mfSurv))
  if (attr(TermsSurv, 'intercept')) Wtime <- as.matrix(Wtime[, - 1])
  # design matrix in survival part, one row for each subject, excluding intercept #
  
  TermsLongX <- fitLME$terms
  mydata <- fitLME$data[all.vars(TermsLongX)] 
  formLongX <- formula(fitLME) 
  mfLongX <- model.frame(TermsLongX, data = mydata) 
  X <- as.matrix(model.matrix(formLongX, mfLongX))
  alpha.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]
  
  formLongZ <- formula(fitLME$modelStruct$reStruct[[1]]) 
  mfLongZ <- model.frame(terms(formLongZ), data = mydata)
  TermsLongZ <- attr(mfLongZ, "terms") 
  Z <- as.matrix(model.matrix(formLongZ, mfLongZ)) 
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed #
  
  data.id <- mydata[!duplicated(ID), ] # pick the first row of each subject in mydata, nrow=n #
  
  if (!is.null(timeVarY)) {
    if (!all(timeVarY %in% names(mydata)))
      stop("\n'timeVarY' does not correspond to columns in the fixed-effect design matrix of 'fitLME'.")
    data.id[timeVarY] <- Time
  }
  
  mfLongX.id <- model.frame(TermsLongX, data = data.id)
  Xtime <- as.matrix(model.matrix(formLongX, mfLongX.id)) # same structure with X, but with only n rows #
  
  mfLongZ.id <- model.frame(TermsLongZ, data = data.id)
  Ztime <- as.matrix(model.matrix(formLongZ, mfLongZ.id)) # same structure with Z, but with only n rows #
  
  U <- sort(unique(Time[d == 1])) # ordered uncensored observed event time #
  tempU <- lapply(Time, function(t) U[t >= U]) 
  times <- unlist(tempU) # vector of length M #
  nk <- sapply(tempU, length)  # length of each element in times, vector of length n #
  M <- sum(nk)
  Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Index0 <- match(Time, U)
  Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu #
  
  data.id2 <- data.id[Index, ]
  if (!is.null(timeVarY)) {
    data.id2[timeVarY] <- times
  }
  mfLongX2 <- model.frame(TermsLongX, data = data.id2)
  Xtime2 <- as.matrix(model.matrix(formLongX, mfLongX2))
  mfLongZ2 <- model.frame(TermsLongZ, data = data.id2)
  Ztime2 <- as.matrix(model.matrix(formLongZ, mfLongZ2))
  
  mfSurv2 <- mfSurv[Index, ]
  if (!is.null(timeVarT)){
    mfSurv2[timeVarT] <- times
  }
  Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
  if (attr(TermsSurv, 'intercept')) Wtime2 <- as.matrix(Wtime2[, - 1]) # excluding intercept #
  
  n <- nLong
  N <- length(Y)
  ni <- as.vector(tapply(ID, ID, length))
  nu <- length(U)
  ncz <- ncol(Z)
  ncx <- ncol(X)
  ncw <- ncol(W)
  ncz2 <- ncz ^ 2
  p <- ncz * (ncz + 1) / 2
  
  controlvals <- list(tol.P = 10 ^ (-4), tol.L = 10 ^ (-8), max.iter = 200, SE.method = 'PRES', delta = 10 ^ (- 5), 
                      nknot = if (ncz == 1) 12 else if (ncz == 2) 10 else 8)
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
  b <- as.matrix(expand.grid(rep(list(GHQ$nodes), ncz)))
  wGQ <- as.matrix(expand.grid(rep(list(GHQ$weights), ncz)))
  wGQ <- apply(wGQ, 1, prod)
  GQ <- nrow(b)
  
  Z.st <- lapply(split(Z, ID), function(x) matrix(x, ncol = ncz))
  Y.st <- split(Y, ID)
  X.st <- lapply(split(X, ID), function(x) matrix(x, ncol = ncx))
  Ztime2.st <- vector('list', n)
  for (i in (1:n)[nk != 0]) { Ztime2.st[[i]] <- matrix(Ztime2[Index == i, ], ncol = ncz) }
  Wtime22 <- if (ncw > 1) t(apply(Wtime2, 1, function(x) x %o% x)) else Wtime2 ^ 2
  Xtime22 <- if (ncx > 1) t(apply(Xtime2, 1, function(x) x %o% x)) else Xtime2 ^ 2
  X2 <- if(ncx > 1) t(apply(X, 1, function(x) x %o% x)) else X ^ 2
  X2.sum <- matrix(colSums(X2), nrow = ncx)  
  
  if (model == 1) {
    environment(InitValTM1) <- environment(EMiterTM1) <- environment()
  } else {
    environment(InitValTM2) <- environment(EMiterTM2) <- environment()
  }
    
  BSigma <- lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                   function(x) x * fitLME$sigma ^ 2)[[1]]
  # the estimated variance-covariance matrix for the random effects  #
  beta <- as.vector(fixef(fitLME))
  beta.names <- names(fixef(fitLME))
  Ysigma <- fitLME$sigma
  
  surv.init <- if (model == 1) InitValTM1(beta) else InitValTM2(beta)
  phi <- surv.init$phi
  alpha <- surv.init$alpha
  lamb <- surv.init$lamb
  
  theta.old <- list(beta = beta, phi = phi, alpha = alpha, Ysigma = Ysigma, BSigma = BSigma,  
                    lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1
  
  while (step <= iter) {
    
    if (err.P < tol.P | err.L < tol.L) break
    
    theta.new <- if (model == 1) EMiterTM1(theta.old) else EMiterTM2(theta.old)
    
    new.P <- c(theta.new$beta, theta.new$phi, theta.new$alpha, theta.new$Ysigma, theta.new$BSigma)
    old.P <- c(theta.old$beta, theta.old$phi, theta.old$alpha, theta.old$Ysigma, theta.old$BSigma)
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
  environment(Sfunc) <- environment()
  if (model == 1) {
    environment(Lamb1) <- environment(DQfunc1) <- environment(LH1) <- environment()
  } else {
    environment(Lamb2) <- environment(DQfunc2) <- environment(LH2) <- environment()
  }
  if (controlvals$SE.method == 'PFDS') {
    environment(PFDS) <- environment()
    if (CheckDeltaFD(theta.new, ncz, delta)) {
      time.SE <- system.time(Vcov <- PFDS(model, theta.new, min(tol.P, delta)/100, iter, delta))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (controlvals$SE.method == 'PRES') {
    environment(PRES) <- environment()
    if (CheckDeltaRE(theta.new, ncz, delta)) {
      time.SE <- system.time(Vcov <- PRES(model, theta.new, min(tol.P, delta)/100, iter, delta))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (controlvals$SE.method == 'PLFD') {
    environment(PLFD) <- environment()
    if (CheckDeltaFD(theta.new, ncz, delta)) {
      time.SE <- system.time(Vcov <- PLFD(model, theta.new, min(tol.P, delta)/100, iter, delta))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else {
    Vcov <- time.SE <- NA
    warning("\n Standard error estimation method should be either 'PFDS', 'PRES' or 'PLFD'.")
  }
  
  
  theta.new$lamb <- cbind("time" = U, "bashaz" = theta.new$lamb)
  names(theta.new$beta) <- beta.names
  names(theta.new$phi) <- phi.names
  names(theta.new$alpha) <- if (model == 1) alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"
  if (ncz > 1) dimnames(theta.new$BSigma) <- dimnames(BSigma) 
  else names(theta.new$BSigma) <- "sigma.b"
  
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
  class(result) <- "jmodelTM"
  result
}
