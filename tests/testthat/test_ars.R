context("ars")

tol = 1e-1

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

test_that("test deriv", {
  # the deriv of the log normal dist d/dx log(e^(-x^2/2)/sqrt(2 π))= -x
  x=3
  result = deriv(dnorm, x)
  expected = -x
  expect_equal(result,expected, tolerance=1e-1)
  
  #unif should be 0
  x=.3
  result = deriv(dunif, x)
  expected = 0
  expect_equal(result,expected, tolerance=1e-1)
})

test_that("test cdf", {
  x = 3
  cdf = get_cdf(dnorm)
  result = cdf(x)
  expected = pnorm(x)
  expect_equal(result,expected, tolerance=tol)
  
  x = 3
  n = 100
  p=.2
  binom_func = function(x) dbinom(x, n, p)
  cdf = get_cdf(binom_func)
  result = suppressWarnings(cdf(x))
  expected = pbinom(x,n, p)
  expect_equal(result,expected, tolerance=tol)
  
})

test_that("test calc tangents", {
  # d/dx of log(std normal) = -x
  abscis = c(-1,0)
  
  #tangent to 0 is  y = -1/2 log(2 π)
  #and tangent line to -1 is  = x + 1 + 1/2 (-1 - log(2 π)) 
  # and they intersect at x=-.5
  expected = -.5
  result = calculate_tangents(abscis, dnorm)
  expect_equal(result,expected, tolerance=tol)
  
})