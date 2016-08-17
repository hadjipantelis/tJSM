 cat("\nTests for 'jmodelMult'")

 myEps <- .Machine$double.eps
 
fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list( max.iter = 50, nknot = 5, tol.L = 1e-04, tol.P = 1e-02)
 
#test_that(" basic jmodelMult test with for aids data model = 1, rho = 0 ", { 
#  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=0,   control = control)
#  expect_equal(      m_MULT$coefficients$lgLik,-2548.893797585262746, tolerance = (10^4)*myEps, scale = 1)
#  expect_equal( mean (m_MULT$coefficients$lamb), 4.353943942192298, tolerance = (10^1)*myEps, scale = 1)
#})

test_that(" basic jmodelMult test with for aids data model = 1, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2551.601832947596904, tolerance = 3, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 4.374275722571750, tolerance = (10^14)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for aids data model = 2, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=0, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2542.155094254993401, tolerance = 3, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 4.358013652751516, tolerance = (10^14)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for aids data model = 2, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2542.768126923503132, tolerance = 4, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb),4.393738379111060, tolerance = 0.03, scale = 1)
})
