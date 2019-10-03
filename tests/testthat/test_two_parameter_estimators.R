context("Test estimators for bivariate parameters")
test_that("rob_estimate_negbin returns a vector of two non negative real number", {
  x <- rnbinom(100,mu=1,size=2)
  estimate <- rob_estimate_negbin(x)
  expect_is(estimate, "numeric")
  expect_equal(estimate[2]>=0,TRUE)
  expect_equal(estimate[1]>=0,TRUE)
})
