 cat("\nTests for 'jmodelMult'")

 myEps <- .Machine$double.eps

#load("../../../data/aids.rda")
fitLME <- lme(sqrt(CD4) ~ bs(obstime, 4), random =~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list( max.iter = 100, nknot = 5)
 
test_that(" basic jmodelMult test with for aids data model = 1, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=0,   control = control)
  expect_equal(      m_MULT$coefficients$lgLik,-2548.893797585262746, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 4.353943942192298, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for aids data model = 1, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 1, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2551.601832947596904, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 4.374275722571750, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for aids data model = 2, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=0, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2542.155094254993401, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 4.358013652751516, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for aids data model = 2, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, aids, model = 2, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2542.768126923503132, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb),4.393738379111060, tolerance = (10^1)*myEps, scale = 1)
})
 
#load("../../../data/pbc.rda")
fitLME <- lme(log(serBilir) ~ bs(obstime, df = 6, degree = 2), random = ~ 1 | ID, data = pbc)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = pbc, x = TRUE) 
control <- list(tol.P = 10 ^ (- 3), max.iter = 100, nknot = 10)

test_that(" basic jmodelMult test with for pbc data model = 1, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 1, rho=0, control = control)
  expect_equal(      m_MULT$coefficients$lgLik, -2602.822699596125403, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 2.207232400382697, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for pbc data model = 1, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 1, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2611.630755259258876, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb),2.208544472332525, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for pbc data model = 2, rho = 0 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho=0, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2536.769504475529629, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb),  2.206643985199276, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelMult test with for pbc data model = 2, rho = 1 ", { 
  m_MULT <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho=1, control = control)
  expect_equal(  m_MULT$coefficients$lgLik, -2541.720928353124691, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_MULT$coefficients$lamb), 2.207320505861047, tolerance = (10^1)*myEps, scale = 1)
})

