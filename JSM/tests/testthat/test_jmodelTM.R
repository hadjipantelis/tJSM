 cat("\nTests for 'jmodelTM'")

 myEps <- .Machine$double.eps

fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list(nknot = 15, tol.L = 1e-04, tol.P = 1e-02)

test_that(" basic jmodelTM test with for aids data model = 1, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2523.329037070506274, tolerance = 2, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.379207628026059, tolerance = 0.02, scale = 1)
})

test_that(" basic jmodelTM test with for aids data model = 2, rho = 0 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2521.566933974915173, tolerance = 0.07, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.340359673028605, tolerance = (10^11)*myEps, scale = 1)
})

test_that(" basic jmodelTM test with for aids data model = 2, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2522.506475889341800, tolerance = (10^15)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.341256927560805, tolerance = (10^11)*myEps, scale = 1)
})

