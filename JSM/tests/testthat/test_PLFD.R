 cat("\nTests for 'PLFD (jmodelMult)'")

 myEps <- .Machine$double.eps

#load("data/aids.rda")
fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list( max.iter = 100, nknot = 5, SE.method ='PLFD')

# We use increased tolerances as this procedure is very susceptible to numeric noise. 
 
test_that(" basic PLFD jmodelMult test with for aids data model = 1, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=0,   control = control)
  expect_equal( mean (m_MULT$est.bi), 0.99912665200230888, tolerance = (10^2)*myEps, scale = 1)
  expect_equal( mean (m_MULT$Vcov), 0.00050059556321931, tolerance = (10^5)*myEps, scale = 1)
})


test_that(" basic PLFD jmodelMult test with for aids data model = 2, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=1, control = control)
  expect_equal( mean (m_MULT$est.bi), 0.99979161336117106, tolerance = (10^2)*myEps, scale = 1)
  expect_equal( mean (m_MULT$Vcov), 0.00198699228755987, tolerance = (10^6)*myEps, scale = 1)
})
 
