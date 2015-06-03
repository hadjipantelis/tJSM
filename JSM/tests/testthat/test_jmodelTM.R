 cat("\nTests for 'jmodelTM'")

 myEps <- .Machine$double.eps

#load('../../../data/aids.rda')
fitLME <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^2), random = ~ 1 | ID, data = aids)
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)
control <- list(nknot = 15)

#test_that(" basic jmodelTM test with for aids data model = 1, rho = 0 ", { 
#  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=0,timeVarY = 'obstime', control = control)
#  expect_equal(      m_TM$coefficients$lgLik, -2522.306954085800498, tolerance = (10^4)*myEps, scale = 1)
#  expect_equal( mean (m_TM$coefficients$lamb), 4.354565736461517, tolerance = (10^1)*myEps, scale = 1)
#})

test_that(" basic jmodelTM test with for aids data model = 1, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 1, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2523.329037070506274, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.379207628026059, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelTM test with for aids data model = 2, rho = 0 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=0,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2521.566933974915173, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.340359673028605, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelTM test with for aids data model = 2, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik,-2522.506475889341800, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 4.341256927560805, tolerance = (10^1)*myEps, scale = 1)
})

if(1==3){

#load("../../../data/liver.rda")
fitLME <- lme(proth ~ Trt * obstime, random = ~ obstime | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE) 
control <- list(tol.P = 10 ^ (- 3))

test_that(" basic jmodelTM test with for liver data model = 1, rho = 0 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=0,timeVarY = 'obstime', control = control)
  expect_equal(      m_TM$coefficients$lgLik, -15054.946130276501208 , tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 1.584051297495303, tolerance = (10^1)*myEps, scale = 1)
})
test_that(" basic jmodelTM test with for liver data model = 1, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, liver, model = 1, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik, -15053.795691753390201, tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 2.790275789183873, tolerance = (10^1)*myEps, scale = 1)
})

test_that(" basic jmodelTM test with for liver data model = 2, rho = 1 ", { 
  m_TM <- jmodelTM(fitLME, fitCOX, liver, model = 2, rho=1,timeVarY = 'obstime', control = control)
  expect_equal(  m_TM$coefficients$lgLik, -15053.524758074829151 , tolerance = (10^4)*myEps, scale = 1)
  expect_equal( mean (m_TM$coefficients$lamb), 1.513523123979168 , tolerance = (10^1)*myEps, scale = 1)
})

}
