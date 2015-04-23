 cat("\nTests for 'PFDS (jmodelMult)'")

 myEps <- .Machine$double.eps

load("aids.rda")
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
  expect_equal( mean (m_MULT$est.bi), 0.99979161336117106, tolerance = (10^1)*myEps, scale = 1)
  expect_equal( mean (m_MULT$Vcov), 0.00202229509590529, tolerance = (10^2)*myEps, scale = 1)
})
 
