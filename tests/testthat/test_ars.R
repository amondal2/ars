context("ars")

tol = 1e-5

test_that("test inverse cdf", {
  prob = .8
  result = eval_inverse_cdf(prob, pnorm)
  expected = qnorm(prob)
  expect_equal(result,expected, tolerance=tol)
  
  prob = .3
  pois_func = function(x) ppois(x,5)
  result = eval_inverse_cdf(prob, pois_func)
  expected = qpois(prob,5)
  expect_equal(result,expected, tolerance=tol)
})

