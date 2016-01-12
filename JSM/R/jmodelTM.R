#' Semiparametric Joint Models for Survival and Longitudinal Data
#' 
#' This function applies a maximum likelihood approach to fit the semiparametric joint models of survival and normal longitudinal data. The survival model is assumed to come from a class of transformation models, including the Cox proportional hazards model and the proportional odds model as special cases. The longitudinal process is modeled by liner mixed-effects models.
#' 
#'  @param fitLME An object inheriting from class \code{lme} representing a fitted linear mixed-effects model. See \bold{Note}.
#'  @param fitCOX An object inheriting from class \code{coxph} representing a fitted Cox proportional hazards regression model. Specifying \code{x = TRUE} is required in the call to \code{coxph()} to include the design matrix in the object fit. See \bold{Note}.
#'  @param data A data.frame containing all the variables included in the joint modeling. See \bold{Note}.
#'  @param model An indicator specifying the dependency between the survival and longitudinal outcomes. Default is 1. See\bold{Details}.
#'  @param rho A nonnegative real number specifying the transformation model you would like to fit. Default is 0, i.e. the Cox proportional hazards model. See \bold{Details}.
#'  @param timeVarY A character string indicating the time variable in the linear mixed-effects model. See \emph{Examples}.
#'  @param timeVarT A character string indicating the time variable in the \code{coxph} object. Normally it is \code{NULL}. See \bold{Note} and \emph{Examples}.
#'  @param control A list of control values for the estimation algorithm with components:
#'    \describe{
#'      \item{tol.P}{tolerance value for convergence in the parameters with default value 1e-04. See \bold{Details}.}
#'      \item{tol.L}{tolerance value for convergence in the log-likelihood with default value 1e-08. See \bold{Details}.}
#'      \item{max.iter}{the maximum number of EM iterations with default value 200.}
#'      \item{SE.method}{a character string specifying the standard error estimation method. Default is \code{"PRES"}. See \bold{Details} and \bold{Note}.}
#'      \item{delta}{a positive value used for numerical differentiation in the \code{SE.method}. Default is 1e-05 if \code{"PRES"} is used and 1e-03 otherwise. See \bold{Details}.}
#'      \item{nknot}{the number of Gauss-Hermite quadrature knots used to approximate the integrals over the random effects. Default is 12 and 10 for one- and two-dimensional integration, respectively, and 8 for those with higher dimensions.}
#'    }
#' @param \dots Additional options to be passed to the \code{control} argument.
#' 
#' @details The \code{jmodelTM} function fits joint models for survival and longitudinal data. Linear mixed-effects models are assumed for the longitudinal processes. With the Cox proportional hazards model and the proportional odds model as special cases, a general class of transformation models are assumed for the survival processes. The baseline hazard functions are left unspecified, i.e. no parametric forms are assumed, thus leading to semiparametric models. For detailed model formulation, please refer to Xu, Baines and Wang (2014). 
#' 
#' The longitudinal model is written as \deqn{Y_i(t)=\mu_i(t) + \varepsilon_i(t) = \mathbf{X}_i^\top(t)\boldsymbol\beta + \mathbf{Z}_i^\top(t)\mathbf{b}_i + \varepsilon_i(t).} If \code{model = 1}, then the linear predictor for the survival model is expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha\mu_i(t),} indicating that the entire longitudinal process (free of error) enters the survival model as a covariate. If other values are assigned to the \code{model} argument, the linear predictor for the surival model is then expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha\mathbf{Z}_i^\top(t)\mathbf{b}_i,} suggesting that the survival and longitudinal models only share the same random effects.
#' 
#' The survival model is written as \deqn{\Lambda(t|\eta(t)) = G\left[\int_0^t\exp{\eta(s)}d\Lambda_0(s)\right],} where \eqn{G(x) = \log(1 + \rho x) / \rho} with \eqn{\rho \geq 0} is the class of logarithmic transfomrations. If \code{rho = 0}, then \eqn{G(x) = x}, yielding the Cox proportional hazards model. If \code{rho = 1}, then \eqn{G(x) = \log(1 + x)}, yielding the proportional odds model. Users could assign any nonnegative real value to \code{rho}. 
#' 
#' An expectation-maximization (EM) algorithm is implemented to obtain parameter estimates. The convergence criterion is either of (i) \eqn{\max \{ | \boldsymbol\theta^{(t)} - \boldsymbol\theta^{(t - 1)} | / ( | \boldsymbol\theta^{(t - 1)} | + tol.P) \} < tol.P}, or (ii) \eqn{| L(\boldsymbol\theta^{(t)}) - L(\boldsymbol\theta^{(t - 1)})| / ( | L(\theta^{(t - 1)}) | + tol.P ) < tol.L}, is satisfied. Here \eqn{\boldsymbol\theta^{(t)}} and \eqn{\boldsymbol\theta^{(t-1)}} are the vector of parameter estimates at the \eqn{t}-th and \eqn{(t-1)}-th EM iterations, respectively; \eqn{L(\boldsymbol\theta)} is the value of the log-likelihood function evaluated at \eqn{\boldsymbol\theta}. Users could specify the tolerance values \code{tol.P} and \code{tol.L} through the \code{control} argument.
#' 
#' For standard error estimation for the parameter estimates, three methods are provided, namely \code{"PRES"}, \code{"PFDS"} and \code{"PLFD"} (detailed information are referred to Xu, Baines and Wang (2014)). In the \code{control} argument, if \code{SE.method = "PRES"}, numerically differentiating the profile Fisher score vector with Richardson extrapolation is applied; if \code{SE.method = "PFDS"}, numerically differentiating the profile Fisher score vector with forward difference is applied; if \code{SE.method = "PLFD"}, numerially (second) differentiating the profile likelihood with forward difference is applied. Generally, numerically differentiating a function \eqn{f(x)} (an arbitrary function) with forward difference is expressed as \deqn{f^\prime(x) = \frac{f(x + \delta) - f(x)}{\delta},} and that with Richardson extrapolation is expressed as \deqn{f^\prime(x) = \frac{f(x - 2\delta) - 8f(x - \delta) + 8f(x + \delta) - f(x + 2\delta)}{12\delta}.} Users could specify the value of \eqn{\delta} through the \code{delta} item in the \code{control} argument.
#' 
#' 
#' @return See \code{\link{jmodelTMObject}} for the components of the fit.
#' 
#' 
#' @examples
#' # linear mixed-effects model fit with random intercept
#' fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), 
#'              random = ~ 1 | ID, data = aids)
#' # Cox proportional hazards model fit with a single time-independent covariate
#' fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
#' 
#' # joint model fit which assumes the Cox proportional hazards model for the survival process
#' fitJT.ph <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime')
#' summary(fitJT.ph)
#' # joint model fit where the survival and longitudinal processes only share the same random effect
#' fitJT.ph2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, timeVarY = 'obstime')
#' summary(fitJT.ph2)
#' # joint model fit with standard error estimates calculated by the 'PFDS' method
#' fitJT.ph3 <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime', control = list(SE.method = 'PFDS'))
#' summary(fitJT.ph3)
#' 
#' # joint model fit which assumes the proportional odds model for the survival process
#' fitJT.po <- jmodelTM(fitLME, fitCOX, aids, rho = 1, timeVarY = 'obstime')
#' summary(fitJT.po)
#' # joint model fit where the survival and longitudinal processes only share the same random effect
#' fitJT.po2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho = 1, timeVarY = 'obstime')
#' summary(fitJT.po2)
#' 
#' # linear mixed-effects model fit with random intercept and random slope
#' fitLME2 <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), 
#'                random = ~ obstime | ID, data = aids)
#' # Cox proportional hazards model fit with a time-dependent covariate
#' fitCOX2 <- coxph(Surv(start, stop, event) ~ drug + as.numeric(drug) : obstime, data = aids, x = TRUE)
#' # joint model fit in which \code{timeVarT} must be specified
#' fitJT.ph4 <- jmodelTM(fitLME2, fitCOX2, aids, timeVarY = 'obstime', timeVarT = 'obstime')
#' summary(fitJT.ph4)
#' 
#' @references
#' \cite{Dabrowska, D. M. and Doksun K. A. (1988) Partial Likelihood in Transformation Models with Censored Data.  Scandinavian Journal of Statistics, 15, 1--23.} 
#' 
#' \cite{Tsiatis, A. A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. Statistica Sinica, 14, 809--834. } 
#' 
#' \cite{Wulfsohn, M. S. and Tsiatis, A. A. (1997) A joint model for survival and longitudinal data measured with error. Biometrics, 53, 330--339.}
#' 
#' \cite{Xu, C., Baines, P. D. and Wang, J. L. (2014) Standard error estimation using the EM algorithm for the joint modeling of survival and longitudinal data. Biostatistics. doi: 10.1093/biostatistics/kxu015}\#' 
#' 
#' \cite{Zeng, D. and Lin, D. (2007) Maximum likelihood estimation in semiparametric regression models with censored data.  Journal of the Royal Statistical Society: Series B, 69, 507--564}
#' 
#' @seealso \code{\link{jmodelTMObject}}, \code{\link{lme}}, \code{\link{coxph}},\code{\link{Surv}}
#' 
#' @note 
#'    1. Currently, \code{jmodelTM()} could only handle the \code{fitLME} object with a simple random-effects structure (only the \code{pdDiag()} class). Moreover, the within-group correlation and heteroscedasticity structures in the \code{fitLME} object (i.e. the \code{correlation} and \code{weights} argument of \code{lme()}) are ignored.
#'   
#'    2. The \code{data} argument in \code{jmodelTM()}, \code{lme()} and \code{coxph()} should be the same data frame.
#'  
#'    3. For the \code{fitCOX} object, only the \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) should be involved in the \code{formula} argument of \code{coxph{}}. Since \code{coxph()} uses the same data frame as \code{lme()} does, a time-dependent Cox model must be fitted by \code{coxph()} although \eqn{W_i(t)} may only contain time-independent covariates. See \emph{Examples}.
#'    
#'    4. If \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) does involve time-dependent covariate, then \code{timeVarT} must specify the name of the time variable involved. See \emph{Examples}.
#'    
#'    5. The standard error estimates are obtained by numerical approximations which is naturally subject to numerical errors. Therefore, in extreme cases, there may be \code{NA} values for part of the standard error estimates.
#' 
#' @keywords survival
#' @keywords longitudinal
#' 
#' @export


jmodelTM <- function (fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, timeVarT = NULL, 
                      control = list(), ...)
{
  call <- match.call()
  
  CheckInputs(fitLME, fitCOX, rho)
  
  cntrlLst <- GenerateControlList(control)  
  
  ID <- as.vector(unclass(fitLME$groups[[1]]))     # grouping factors as a vector
  ni <- as.vector(tapply(ID, ID, length))          # number of each factor as vector
  bBLUP <- data.matrix(ranef(fitLME))              # BLUP estimate of the random effects
  dimnames(bBLUP) <- NULL
  nLong <- nrow(bBLUP)                             # number of groupings in the longitudinal process
  if(ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")
  start <- as.vector(fitCOX$y[, 1])                # starting time of the interval which contains the time of measurements
  stop <- as.vector(fitCOX$y[, 2])                 # ending time of the interval which contains the time of measurements
  event <- as.vector(fitCOX$y[, 3])                # event indicator suggesting whether the event-of-interest happens in the interval given
  Time <- stop[cumsum(ni)]                         # survival time, i.e. time to death or censoring for each subject
  d <- event[cumsum(ni)]                           # event indicator for each subject
  nSurv <- length(Time)                            # number of event measurements
  if(sum(d) < 5)
    warning("\n more than 5 events are required.")
  if(nLong != nSurv)
    stop("\n sample sizes in the longitudinal and event processes differ.")
  
  W <- as.matrix(fitCOX$x)                         # observed covariates including baseline covariates as well as longitudinal
  
  varNames <- list()
  varNames$phi.names <- colnames(W)                # names of covariates
  formSurv <- formula(fitCOX)                      # formula of survival model
  TermsSurv <- fitCOX$terms                        # the ‘terms’ object used by the 'fitCOX'
  mfSurv <- model.frame(TermsSurv, data)[cumsum(ni), ]  # fixed-effect design matrix of 'fitCOX'
  if(!is.null(timeVarT)) {
    if(!all(timeVarT %in% all.vars(TermsSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  Wtime <- as.matrix(model.matrix(formSurv, mfSurv)) # design matrix for survival part
  if(attr(TermsSurv, 'intercept')) Wtime <- as.matri(x(Wtime[, - 1])
  # design matrix in survival part, one row for each subject, excluding intercept #
  
  TermsLongX <- fitLME$terms                       # the ‘terms’ object used by the 'fitLME'
  mydata <- fitLME$data[all.vars(TermsLongX)]      # get all the data in the LME object
  formLongX <- formula(fitLME)                     # formula of linear mixed model
  mfLongX <- model.frame(TermsLongX, data = mydata)# fixed-effect design matrix of 'fitLME'
  X <- as.matrix(model.matrix(formLongX, mfLongX)) # fixed-effect design matrix for both models
  varNames$alpha.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]
  
  formLongZ <- formula(fitLME$modelStruct$reStruct[[1]]) 
  mfLongZ <- model.frame(terms(formLongZ), data = mydata)
  TermsLongZ <- attr(mfLongZ, "terms") 
  Z <- as.matrix(model.matrix(formLongZ, mfLongZ)) 
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed #
  
  data.id <- mydata[!duplicated(ID), ] # pick the first row of each subject in mydata, nrow=n #
  
  if(!is.null(timeVarY)) {
    if(!all(timeVarY %in% names(mydata)))
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
  
  Indcs <- list(); 
  
  Indcs$Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Indcs$Index0 <- match(Time, U)
  Indcs$Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Indcs$Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu #
  
  data.id2 <- data.id[Indcs$Index, ]
  if(!is.null(timeVarY)) {
    data.id2[timeVarY] <- times
  }
  mfLongX2 <- model.frame(TermsLongX, data = data.id2)
  Xtime2 <- as.matrix(model.matrix(formLongX, mfLongX2))
  mfLongZ2 <- model.frame(TermsLongZ, data = data.id2)
  Ztime2 <- as.matrix(model.matrix(formLongZ, mfLongZ2))
  
  mfSurv2 <- mfSurv[Indcs$Index, ]
  if(!is.null(timeVarT)){
    mfSurv2[timeVarT] <- times
  }
  Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
  if(attr(TermsSurv, 'intercept')) Wtime2 <- as.matrix(Wtime2[, - 1]) # excluding intercept #
  
  n <- nLong
  N <- length(Y)
  ni <- as.vector(tapply(ID, ID, length))
  nu <- length(U)
  ncz <- ncol(Z)
  ncx <- ncol(X)
  ncw <- ncol(W)
  ncz2 <- ncz ^ 2
  p <- ncz * (ncz + 1) / 2
  
  GHQ <- gauss.quad(cntrlLst$nknot, kind = "hermite")
  b <- as.matrix(expand.grid(rep(list(GHQ$nodes), ncz)))
  wGQ <- as.matrix(expand.grid(rep(list(GHQ$weights), ncz)))
  wGQ <- apply(wGQ, 1, prod)
  GQ <- nrow(b)
  
  Z.st <- lapply(split(Z, ID), function(x) matrix(x, ncol = ncz))
  Y.st <- split(Y, ID)
  X.st <- lapply(split(X, ID), function(x) matrix(x, ncol = ncx))
  Ztime2.st <- vector('list', n)
  for (i in (1:n)[nk != 0]) { Ztime2.st[[i]] <- matrix(Ztime2[Indcs$Index == i, ], ncol = ncz) }
  Wtime22 <- if(ncw > 1) t(apply(Wtime2, 1, function(x) tcrossprod(x))) else Wtime2 ^ 2
  Xtime22 <- if(ncx > 1) t(apply(Xtime2, 1, function(x) tcrossprod(x))) else Xtime2 ^ 2
  X2 <- if(ncx > 1) t(apply(X, 1, function(x) tcrossprod(x))) else X ^ 2
  X2.sum <- matrix(colSums(X2), nrow = ncx)  
  
  Bsigma <- lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                   function(x) x * fitLME$sigma ^ 2)[[1]]
  
  # the estimated variance-covariance matrix for the random effects  #
  beta <- as.vector(fixef(fitLME))
  varNames$beta.names <- names(fixef(fitLME))
  Ysigma <- fitLME$sigma
  
  surv.init <- InitValTMGeneric(beta, model = model, n = n, X = X, Z = Z, bBLUP = bBLUP, ID = ID, Xtime = Xtime, Ztime = Ztime, Xtime2 = Xtime2, Ztime2 = Ztime2, Indcs = Indcs, start = start, event = event, stop = stop, W = W, ncw = ncw, Wtime2 = Wtime2, rho= rho, nk = nk, Wtime22 = Wtime22, d = d,  Wtime = Wtime, cvals = cntrlLst)
  phi <- surv.init$phi
  alpha <- surv.init$alpha
  lamb <- surv.init$lamb
  
  theta.old <- list(beta = beta, phi = phi, alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma,  
                    lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1
  
  while (step <= cntrlLst$max.iter) {
    
    if(err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L) break
    
    theta.new  <-   EMiterTMGeneric(theta.old, n = n, Z.st = Z.st, Ztime = Ztime, Ztime2.st = Ztime2.st, nk = nk, Indcs = Indcs, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ, rho = rho, wGQ = wGQ, d = d, Y.st = Y.st, X.st = X.st, ncz = ncz, ncz2 = ncz2, b = b, model =  model, Wtime = Wtime, Xtime = Xtime, X = X, Y = Y, ID = ID, N = N, ncw = ncw, Wtime22 = Wtime22, ncx = ncx, Xtime22 = Xtime22, Z = Z, X2.sum = X2.sum)
    new.P <- c(theta.new$beta, theta.new$phi, theta.new$alpha, theta.new$Ysigma, theta.new$Bsigma)
    old.P <- c(theta.old$beta, theta.old$phi, theta.old$alpha, theta.old$Ysigma, theta.old$Bsigma)
    err.P <- max(abs(new.P - old.P) / (abs(old.P) + cntrlLst$tol.P))
    # add tol.P to avoid zero value of the estimated parameters #
    
    new.L <- theta.new$lgLik
    old.L <- theta.old$lgLik
    err.L <- abs(new.L - old.L) / (abs(old.L) + cntrlLst$tol.P)
    
    step <- step + 1
    theta.old <- theta.new
  }
  converge <- as.numeric(err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L)
  
  if(cntrlLst$SE.method == 'PFDS') {
    if(CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PFDS(model, theta.new, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Indcs = Indcs, n = n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, rho = rho, d = d, wGQ = wGQ, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, X = X, Y = Y, ID = ID, N = N, Z = Z))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if(cntrlLst$SE.method == 'PRES') {
    if(CheckDeltaRE(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PRES(model, theta.new, ncz = ncz, ncx = ncx, ncw = ncw, n = n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, b = b, Ztime = Ztime, Ztime2.st = Ztime2.st, nk = nk, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, rho = rho, Indcs = Indcs, wGQ =wGQ, GQ = GQ, d = d, p = p, ncz2 = ncz2, X = X, Y = Y, Z = Z, ID = ID, N = N, cvals = cntrlLst, varNames = varNames))[3]
      if(any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if(cntrlLst$SE.method == 'PLFD') {
    if(CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PLFD(model, theta.new, n= n, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Z.st = Z.st, Y.st =Y.st, X.st = X.st, b = b, Ztime = Ztime, nk = nk, Indcs = Indcs, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ, rho = rho, d = d, wGQ = wGQ, Ztime2.st = Ztime2.st ))[3]
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
  names(theta.new$beta) <- varNames$beta.names
  names(theta.new$phi) <- varNames$phi.names
  names(theta.new$alpha) <- if(model == 1) varNames$alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"
  if(ncz > 1) dimnames(theta.new$Bsigma) <- dimnames(Bsigma) 
  else names(theta.new$Bsigma) <- "sigma.b"
  
  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if(converge == 1) "success" else "failure"
  result$control <- cntrlLst
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  class(result) <-  unlist(strsplit(deparse(sys.call()), split = '\\('))[1]
  
  return(result)
}
