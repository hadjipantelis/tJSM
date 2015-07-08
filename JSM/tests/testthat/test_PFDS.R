 cat("\nTests for 'PFDS (jmodelMult)'")

 myEps <- .Machine$double.eps
# browser()
cat(getwd())

#load("JSM/data/aids.rda")
fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list( max.iter = 100, nknot = 5, SE.method ='PFDS')
 
test_that(" basic PFDS jmodelMult test with for aids data model = 1, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=0,   control = control)
  expect_equal( mean (m_MULT$est.bi), 0.99912665200230888, tolerance = (10^1)*myEps, scale = 1)
  expect_equal( mean (m_MULT$Vcov),  0.00051055462959423, tolerance = (10^2)*myEps, scale = 1)
})

test_that(" basic PFDS jmodelMult test with for aids data model = 2, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=1, control = control)
  expect_equal( mean (m_MULT$est.bi), 0.99979161336117106, tolerance = (10^3)*myEps, scale = 1)
  expect_equal( mean (m_MULT$Vcov), 0.00202229509590529, tolerance = (10^2)*myEps, scale = 1)
})

 cat("\nTests for 'PFDS (jmodelTM)'")

fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime^2) + drug:obstime + drug:I(obstime^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list(nknot = 6, SE.method = 'PFDS')

test_that(" basic PLFD jmodelTM test with for aids data model = 1, rho = 0 ", {
  m_TM <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime', control= control, model=1, rho=0)
  expect_equal( mean (m_TM$est.bi), -0.000211429242162, tolerance = (10^2)*myEps, scale = 1)
  expect_equal( mean (m_TM$Vcov),  0.000407428109278, tolerance = (10^5)*myEps, scale = 1)
})

test_that(" basic PLFD jmodelTM test with for aids data model = 2, rho = 1 ", {
  m_TM <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime', control= control, model=2, rho=1)
  expect_equal( mean (m_TM$est.bi), -0.000700248613528, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$Vcov), 0.000677869446869, tolerance = (10^4)*myEps, scale = 1)
})

 
