context("Test estimators for univariate parameters")
test_that("rob_estimate_poisson returns a non negative real number", {
  x <- rpois(100,1)
  estimate <- rob_estimate_poisson(x)
  expect_is(estimate, "numeric")
  expect_equal(estimate>=0,TRUE)
  })
